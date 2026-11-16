#!/usr/bin/env python3

"""Reference-anchored gene/oligo comparison on an existing phylogeny. See README.md."""

import argparse

import csv

import gzip

import hashlib

import json

import math

import os

from pathlib import Path

import re

import shutil

import subprocess

import sys

from datetime import datetime, timezone
from concurrent.futures import ThreadPoolExecutor



from Bio import Phylo, SeqIO

from Bio.Seq import Seq



VERSION = '1.4.0'

REF = '__GENE_REFERENCE__'

DNA = set('ACGTRYSWKMBDHVN')

IUPAC = dict(zip('ACGTRYSWKMBDHVN', ['A','C','G','T','AG','CT','GC','AT','GT','AC','CGT','AGT','ACT','ACG','ACGT']))

COMPONENTS = [('fwd','fwdPrimer','Forward RPA','fwdRPA_aligned.fasta'),

              ('rev','revPrimer','Reverse RPA','revRPA_aligned.fasta'),

              ('gRNA','gRNA','gRNA / supplied sequence','gRNA_aligned.fasta')]



def log(msg):

    print(msg, flush=True)



def rc(s):

    return str(Seq(s).reverse_complement())



def read_fasta(path, aligned=False):

    """Validate all records; IDs are the first whitespace-delimited header token."""

    records = {}

    opener = gzip.open if str(path).lower().endswith('.gz') else open

    with opener(path, 'rt') as handle:

        for r in SeqIO.parse(handle, 'fasta'):

            if r.id in records:

                raise ValueError(f'Duplicate FASTA identifier {r.id!r} in {path}')

            s = str(r.seq).upper()

            if not s:

                raise ValueError(f'Empty sequence: {r.id} in {path}')

            bad = set(s) - DNA - ({'-'} if aligned else set())

            if bad:

                raise ValueError(f'Invalid DNA characters {sorted(bad)}: {r.id} in {path}; use DNA (T), not RNA (U).')

            records[r.id] = s

    if not records:

        raise ValueError(f'No FASTA records in {path}')

    if aligned and len({len(s) for s in records.values()}) != 1:

        raise ValueError(f'Inconsistent alignment lengths in {path}')

    return records



def fasta(path, records):

    with open(path, 'w') as out:

        for name, seq in records.items():

            if not seq:

                raise ValueError(f'Refusing to write empty sequence: {name}')

            out.write(f'>{name}\n')

            for i in range(0, len(seq), 80):

                out.write(seq[i:i+80]+'\n')



def table(path, rows, fields):

    with open(path, 'w', newline='') as handle:

        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction='ignore')

        writer.writeheader()

        writer.writerows(rows)



def run(cmd, stdout_path, stderr_path):

    with open(stdout_path, 'w') as out, open(stderr_path, 'w') as err:

        r = subprocess.run([str(x) for x in cmd], stdout=out, stderr=err)

    if r.returncode:

        raise RuntimeError(f'Command failed ({r.returncode}): {cmd[0]}; see {stderr_path}')



def resolve_threads(value):

    if value is not None:

        if value < 1: raise ValueError('--threads must be a positive integer.')

        return value

    if not sys.stdin.isatty():

        raise ValueError('No interactive terminal available; supply --threads N to specify the CPU budget.')

    available = len(os.sched_getaffinity(0)) if hasattr(os, 'sched_getaffinity') else (os.cpu_count() or 1)

    while True:

        try: answer = input(f'How many CPU cores/threads can this run use? (available: {available}): ').strip()

        except EOFError as exc: raise ValueError('No thread count entered; supply --threads N.') from exc

        if answer.isdecimal() and int(answer) > 0: return int(answer)

        log('Please enter a positive whole number of CPU threads.')


def mismatch_positions(ref, seq):

    # A degenerate oligo accepts its concrete members. Ambiguous genome bases

    # are unresolved, never silently considered matches (checked by caller).

    return [i+1 for i,(a,b) in enumerate(zip(ref,seq)) if b not in IUPAC[a]]



def locate_oligo(gene, oligo, max_mismatches=0):

    candidates = []

    for strand in ['+','-']:

        for i in range(len(gene)-len(oligo)+1):

            site = gene[i:i+len(oligo)]

            if strand == '-': site = rc(site)

            if set(site)-set('ACGT'): continue

            mm = len(mismatch_positions(oligo,site))

            if mm <= max_mismatches: candidates.append((mm,i,i+len(oligo),strand))

    if not candidates:

        raise ValueError('Oligo does not map to the reference gene at the allowed reference mismatch threshold.')

    candidates.sort()

    best = [x for x in candidates if x[0] == candidates[0][0]]

    if len(best) != 1:

        raise ValueError(f'Oligo maps ambiguously to {len(best)} equally good reference locations. Supply a unique gene/oligo reference.')

    mm,start,end,strand = best[0]

    return dict(start=start,end=end,strand=strand,reference_mismatches=mm)



def inspect_tree(path):

    # Bio.Phylo can permissively reinterpret malformed branch-length tokens.

    text = Path(path).read_text()

    for token in re.findall(r':\s*([^,();\[\]\s]+)', text):

        if not re.fullmatch(r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?', token):

            raise ValueError(f'Invalid Newick branch length: {token!r}')

    tree = Phylo.read(path,'newick')

    tips = [x.name for x in tree.get_terminals()]

    if not tips or any(not x for x in tips) or len(tips) != len(set(tips)):

        raise ValueError('Tree tips must have nonempty, unique identifiers.')

    for node in tree.find_clades():

        if node.branch_length is not None and (not math.isfinite(node.branch_length) or node.branch_length < 0):

            raise ValueError('Tree branch lengths must be finite and nonnegative.')

    # Midpoint-root so the figure reflects a rooted topology rather than

    # whatever rooting happened to be present in the supplied Newick file.

    # A midpoint requires branch lengths on every edge. Fail before drawing

    # rather than presenting an unrooted tree as midpoint-rooted.

    edges = [node for node in tree.find_clades() if node is not tree.root]

    if len(tips) < 3 or any(node.branch_length is None for node in edges) or not any(node.branch_length > 0 for node in edges):

        raise ValueError('Midpoint rooting requires at least three tips and branch lengths on every edge (with at least one positive length).')

    tree.root_at_midpoint()

    return tree



def bootstrap_support(clade):

    """Read UFBoot from an IQ-TREE SH-aLRT/UFBoot pair, or a single support."""

    for raw in (clade.name, clade.confidence, clade.comment):

        if raw is None:

            continue

        parts = str(raw).strip().split('/')

        if parts and all(re.fullmatch(r'(?:\d+(?:\.\d*)?|\.\d+)', part) for part in parts):

            return float(parts[-1])

    return None


def identifier_parts(identifier):

    """Return an accession-aware key, keeping any supplied assembly version."""

    name = re.sub(r'\.(?:fna|fasta|fa)(?:\.gz)?$', '', identifier.strip(), flags=re.I)

    name = re.sub(r'\s+', '_', name)

    match = re.fullmatch(r'(?P<species>.+)_(?P<accession>GC[AF]_\d+)(?:\.(?P<version>\d+))?', name)

    if match:

        return (match['species'], match['accession'], match['version'])

    return (name, '', None)



def reconcile_identifiers(tipnames, genome_names):

    """Resolve only unique one-to-one matches; never guess between duplicates."""

    unmatched_tips, unmatched_genomes = set(tipnames), set(genome_names)

    matched = {}

    methods = {}

    ambiguous = []



    def stage(method, tip_key, genome_key, compatible=lambda t, g: True):

        for tip in tipnames:

            if tip not in unmatched_tips:

                continue

            key = tip_key(tip)

            candidates = [g for g in unmatched_genomes if genome_key(g) == key and compatible(tip, g)]

            competing = [t for t in unmatched_tips if tip_key(t) == key and any(genome_key(g) == key and compatible(t, g) for g in unmatched_genomes)]

            if len(candidates) == 1 and len(competing) == 1:

                genome = candidates[0]

                matched[tip] = genome

                methods[tip] = method

                unmatched_tips.remove(tip)

                unmatched_genomes.remove(genome)

            elif candidates:

                ambiguous.append((method, tip, tuple(sorted(candidates))))



    stage('exact', lambda s: s, lambda s: s)

    stage('normalized label / FASTA extension', identifier_parts, identifier_parts)

    # If one side omits an accession version, match the species and accession

    # only when they identify exactly one tree tip and one genome. Two explicit

    # but different version numbers never match.

    stage('unique accession with version missing on one side',

          lambda s: identifier_parts(s)[:2], lambda s: identifier_parts(s)[:2],

          lambda t, g: identifier_parts(t)[2] is None or identifier_parts(g)[2] is None)

    unresolved = [(m, t, c) for m, t, c in ambiguous if t in unmatched_tips]

    return matched, methods, sorted(unmatched_tips), sorted(unmatched_genomes), unresolved



def select_hits(lines, contigs, reference_length, args):

    """Chain nonoverlapping, collinear HSPs only within a bounded locus.



    Dynamic programming retains score/coverage Pareto states per start/end HSP.

    Every member independently passes identity; chain coverage never double-counts.

    """

    groups = {}

    for line in lines:

        if not line.strip(): continue

        fields = line.rstrip().split('\t')

        if len(fields) != 10: raise ValueError('Unexpected BLAST output columns.')

        sid,ident,alen,qs,qe,ss,se,bits,evalue,nident = fields

        ident,bits,evalue = float(ident),float(bits),float(evalue)

        alen,qs,qe,ss,se,nident = map(int,(alen,qs,qe,ss,se,nident))

        if sid not in contigs: raise ValueError(f'Unknown BLAST contig identifier: {sid}')

        if not (1 <= min(qs,qe) <= max(qs,qe) <= reference_length):

            raise ValueError('BLAST query coordinates are out of range.')

        low,high = sorted((ss,se))

        if not (1 <= low <= high <= len(contigs[sid])):

            raise ValueError(f'Invalid BLAST coordinates on {sid}: {ss}, {se}')

        if ident < args.min_identity: continue

        strand = '+' if (qe-qs)*(se-ss) >= 0 else '-'

        a,b = (low,high) if strand=='+' else (-high,-low)

        h=dict(qs=min(qs,qe),qe=max(qs,qe),a=a,b=b,score=bits,alen=alen,nident=nident,evalue=evalue)

        groups.setdefault((sid,strand),[]).append(h)

    candidates = []

    max_gap=getattr(args,'max_hsp_gap',1000)

    max_qgap=getattr(args,'max_query_gap',200)

    for (sid,strand),hsps in groups.items():

        hsps.sort(key=lambda h:(h['qs'],h['qe'],h['a']))

        states=[]

        for j,h in enumerate(hsps):

            options=[dict(first=j,path=[j],score=h['score'],covered=h['qe']-h['qs']+1)]

            for i,prev in enumerate(hsps[:j]):

                if not (0 <= h['qs']-prev['qe']-1 <= max_qgap and 0 <= h['a']-prev['b']-1 <= max_gap): continue

                for state in states[i]:

                    if h['b']-hsps[state['first']]['a']+1 > args.max_length_factor*reference_length: continue

                    options.append(dict(first=state['first'],path=state['path']+[j],score=state['score']+h['score'],covered=state['covered']+h['qe']-h['qs']+1))

            # Different starts cannot dominate each other because their spans differ.

            frontier=[]

            for option in sorted(options,key=lambda v:(-v['score'],-v['covered'])):

                if any(x['first']==option['first'] and x['score']>=option['score'] and x['covered']>=option['covered'] for x in frontier): continue

                frontier.append(option)

            states.append(frontier)

            for state in frontier:

                coverage=100*state['covered']/reference_length

                first=hsps[state['first']]

                low,high=(first['a'],h['b']) if strand=='+' else (-h['b'],-first['a'])

                span=high-low+1

                if coverage < args.min_coverage or not args.min_length_factor*reference_length <= span <= args.max_length_factor*reference_length: continue

                chain=[hsps[i] for i in state['path']]

                identity=100*sum(v['nident'] for v in chain)/sum(v['alen'] for v in chain)

                if identity < args.min_identity: continue

                seq=contigs[sid][low-1:high]

                if strand=='-': seq=rc(seq)

                if not seq or set(seq)<={'N'}: continue

                candidates.append(dict(identity=identity,coverage=coverage,contig=sid,start=low,end=high,strand=strand,sequence=seq,bitscore=state['score'],evalue=min(v['evalue'] for v in chain),query_start=first['qs'],query_end=h['qe'],extracted_length=span,hsp_count=len(chain)))

    candidates.sort(key=lambda h:(-h['bitscore'],-h['coverage'],-h['identity'],h['contig'],h['start']))

    # Alternative HSP chains can resolve to precisely the same locus. Keep
    # the highest-scoring chain per coordinate span, but retain other loci.
    distinct=[];seen=set()
    for hit in candidates:
        locus=(hit['contig'],hit['strand'],hit['start'],hit['end'])
        if locus not in seen:
            seen.add(locus);distinct.append(hit)
    return distinct



def search_genome(index, path, reference_path, reference_length, work, logs, args):

    """Search one assembly with isolated temporary files; return its output row."""

    contigs=read_fasta(path)

    idmap={f'contig_{i}':name for i,name in enumerate(contigs)}

    safe={sid:contigs[name] for sid,name in idmap.items()}

    subject=work/f'{index:04d}_genome.fasta'

    hits_path=work/f'{index:04d}_blast.tsv'

    try:

        fasta(subject,safe)

        cmd=[args.blastn,'-task','blastn','-query',reference_path,'-subject',subject,'-dust','no','-soft_masking','false','-word_size','7','-evalue',str(args.evalue),'-max_target_seqs',str(max(100,len(safe))),'-outfmt','6 sseqid pident length qstart qend sstart send bitscore evalue nident']

        run(cmd,hits_path,logs/f'{index:04d}_blast.log')

        with hits_path.open() as handle: candidates=select_hits(handle,safe,reference_length,args)

        row=dict(genome_id=path.name,gene_present=int(bool(candidates)),identity='',coverage='',contig='',start='',end='',strand='',bitscore='',evalue='',query_start='',query_end='',extracted_length='',hsp_count='',valid_candidate_count=len(candidates),selected_candidate_id='',assay_compatible_candidate_count=0)

        for rank,hit in enumerate(candidates,1):
            hit['candidate_id']=f'HIT_{index:04d}_{rank:05d}'
            hit['genome_id']=path.name
            hit['contig']=idmap[hit['contig']]

        return row,candidates

    finally:

        subject.unlink(missing_ok=True)


def site_result(aligned_ref, aligned_gene, oligo, anchor, threshold):

    """Project only reference columns, while rejecting real internal indels."""

    cols = [i for i,c in enumerate(aligned_ref) if c != '-']

    selected = cols[anchor['start']:anchor['end']]

    result = dict(present=0,status='not_detected',mismatches='',mismatch_positions='',

                  sequence='',candidate_sequence='',candidate_mismatches='',insertions=0,deletions=0,

                  start='',end='',strand=anchor['strand'],threshold=threshold)

    if len(selected) != len(oligo): raise ValueError('Reference anchor length is inconsistent.')

    first,last = selected[0],selected[-1]

    seq = ''.join(aligned_gene[i] for i in selected)

    result['insertions'] = sum(aligned_ref[i]=='-' and aligned_gene[i]!='-' for i in range(first,last+1))

    result['deletions'] = seq.count('-')

    non_gap = [i for i,b in enumerate(aligned_gene) if b!='-']

    if not non_gap or first < non_gap[0] or last > non_gap[-1]:

        result['status']='outside_extracted_region'; return result

    if result['insertions'] or result['deletions']:

        result['status']='indel'; return result

    result['start'] = sum(c!='-' for c in aligned_gene[:first])+1

    result['end'] = sum(c!='-' for c in aligned_gene[:last+1])

    if anchor['strand']=='-': seq=rc(seq)

    result['candidate_sequence']=seq

    if set(seq)-set('ACGT'):

        result.update(present='',status='unresolved_bases');return result

    mm=mismatch_positions(oligo,seq)

    result['candidate_mismatches']=len(mm)

    result['mismatch_positions']=';'.join(map(str,mm))

    if len(mm)>threshold:

        result['status']='above_mismatch_threshold';return result

    result.update(present=1,status='detected',mismatches=len(mm),sequence=seq)

    return result


def evaluate_candidate(hit, aligned_ref, aligned_gene, oligos, anchors, args):

    """Evaluate all three sites on one physical locus, never across hits."""

    sites={}
    for key,alias,_,_ in COMPONENTS:
        threshold=getattr(args,key.lower()+'_max_mismatches')
        if threshold is None: threshold=args.max_mismatches
        site=site_result(aligned_ref,aligned_gene,oligos[alias],anchors[key],threshold)
        site['genomic_start']='';site['genomic_end']='';site['genomic_strand']=''
        if site['start']!='':
            if hit['strand']=='+':
                site['genomic_start']=hit['start']+site['start']-1
                site['genomic_end']=hit['start']+site['end']-1
            else:
                site['genomic_start']=hit['end']-site['end']+1
                site['genomic_end']=hit['end']-site['start']+1
            site['genomic_strand']='+' if hit['strand']==site['strand'] else '-'
        sites[key]=site

    f,r,g=sites['fwd'],sites['rev'],sites['gRNA']
    amplicon_length='';geometry='binding_site_unresolved'
    if all(s['start']!='' for s in (f,r,g)):
        if f['strand']==r['strand']:
            geometry='primers_same_strand'
        elif f['strand']=='+' and f['end']<r['start']:
            amplicon_length=r['end']-f['start']+1
            geometry='valid' if f['end']<g['start'] and g['end']<r['start'] else 'guide_outside_primers'
        elif r['strand']=='+' and r['end']<f['start']:
            amplicon_length=f['end']-r['start']+1
            geometry='valid' if r['end']<g['start'] and g['end']<f['start'] else 'guide_outside_primers'
        else:
            geometry='primers_not_inward_facing'
        if geometry=='valid' and not args.min_amplicon_len<=amplicon_length<=args.max_amplicon_len:
            geometry='amplicon_out_of_range'

    detected=sum(s['present']==1 for s in sites.values())
    compatible=geometry=='valid' and detected==3
    assessed=sum(isinstance(s['candidate_mismatches'],int) for s in sites.values())
    mismatch_penalty=sum(s['candidate_mismatches'] if isinstance(s['candidate_mismatches'],int)
                         else len(oligos[alias])+1 for key,alias,_,_ in COMPONENTS
                         for s in [sites[key]])
    rank=(int(compatible),int(geometry=='valid'),detected,assessed,-mismatch_penalty,
          hit['bitscore'],hit['coverage'],hit['identity'])
    return sites,dict(assay_compatible=int(compatible),site_geometry=geometry,
                      amplicon_length=amplicon_length,sites_detected=detected,
                      total_mismatches=mismatch_penalty if assessed==3 else ''),rank



def plot_figure(tree, genes, sites, oligos, output, args, tip_to_genome=None):

    import matplotlib

    matplotlib.use('Agg')

    import matplotlib.pyplot as plt

    from matplotlib.patches import Rectangle

    plt.rcParams.update({'font.family':'DejaVu Sans','pdf.fonttype':42,'ps.fonttype':42})

    tips=tree.get_terminals(); n=len(tips)

    if tip_to_genome is None:

        tip_to_genome={t.name:t.name for t in tips if t.name in genes}

    ys={t:float(n-1-i) for i,t in enumerate(tips)}

    def set_y(node):

        if node not in ys: ys[node]=(set_y(node.clades[0])+set_y(node.clades[-1]))/2

        return ys[node]

    set_y(tree.root)

    depths=tree.depths(); mx=max(depths[t] for t in tips)

    cladogram=mx==0

    if cladogram: depths=tree.depths(unit_branch_lengths=True);mx=max(depths[t] for t in tips) or 1

    def label(name):

        s=re.sub(r'\.(fna|fasta|fa)(\.gz)?$','',name,flags=re.I)

        m=re.match(r'Phytophthora_(.+)_(GC[AF]_.+)$',s)

        return 'P. '+m.group(1).replace('_',' ')+'  '+m.group(2) if m else s

    labels={t:label(t.name) for t in tips}

    left=.35; tree_width=3.0; label_x=left+tree_width+.14

    label_width=max(3.7,max(len(s) for s in labels.values())*.056)

    gene_x=label_x+label_width+.25

    step=.092; panel_x=gene_x+.85; panels=[]

    for key,alias,title,filename in COMPONENTS:

        L=len(oligos[alias]);panels.append((key,alias,title,panel_x,L))

        panel_x+=L*step+.72

    width=panel_x+.18

    row_h=args.row_height; bottom=.75; top=.8

    height=max(4.5,(n-1)*row_h+bottom+top)

    fig=plt.figure(figsize=(width,height),facecolor='white')

    ax=fig.add_axes([0,0,1,1]);ax.set_xlim(0,width);ax.set_ylim(-bottom/row_h,n-1+top/row_h);ax.axis('off')

    # x uses inches; y uses shared leaf rows for branches, labels and all panels.

    ink='#111111';muted='#737373';teal='#087F8C';red='#C23B45'

    # Subtle alternating rows connect labels, gene symbols and assay panels.

    for i,t in enumerate(tips):

        if i % 2 == 0:

            ax.add_patch(Rectangle((label_x-.06,ys[t]-.5),width-label_x-.18,1.0,

                                   facecolor='#EEF0F2',edgecolor='none',zorder=-1))

    dx={c:left+depths[c]/mx*tree_width for c in tree.find_clades()}

    supports=[bootstrap_support(c) for c in tree.find_clades() if c.clades and c is not tree.root]

    supports=[s for s in supports if s is not None]

    support_scale=100 if supports and max(supports)<=1 else 1

    if any(not math.isfinite(s) or s<0 or s*support_scale>100 for s in supports):

        raise ValueError('Bootstrap support must be percentages (0–100) or proportions (0–1), not replicate counts.')

    supported_branches=[]

    def draw_node(node):

        if node.clades:

            ax.plot([dx[node]]*2,[ys[node.clades[0]],ys[node.clades[-1]]],color=ink,lw=args.branch_width)

            for child in node.clades:

                ax.plot([dx[node],dx[child]],[ys[child]]*2,color=ink,lw=args.branch_width)

                if child.clades and bootstrap_support(child) is not None:

                    supported_branches.append((node,child))

                draw_node(child)

    draw_node(tree.root)

    if args.show_support:

        # Put support above the branch it describes, after all tree lines have
        # been drawn. A small white backing keeps crossings from obscuring it.
        for parent,node in supported_branches:

            support=bootstrap_support(node)*support_scale

            x=(dx[parent]+dx[node])/2

            ax.annotate(f'{support:g}',(x,ys[node]),xytext=(0,2.5),
                        textcoords='offset points',ha='center',va='bottom',
                        fontsize=5.5,color=muted,
                        bbox=dict(facecolor='white',edgecolor='none',pad=.4),zorder=5)

    for t in tips:

        y=ys[t];genome_id=tip_to_genome.get(t.name);g=genes.get(genome_id)

        ax.plot([dx[t],label_x-.06],[y,y],color='#909090',lw=args.guide_width,ls=(0,(3,2)),zorder=0)

        label_text=labels[t]

        if label_text.startswith('P. ') and '  ' in label_text:

            species,accession=label_text.split('  ',1)

            label_text=r'$\mathit{'+species.replace(' ',r'\ ')+r'}$  '+accession

        ax.text(label_x,y,label_text,fontsize=7.3,va='center',color=ink)

        if g is None: ax.scatter([gene_x],[y],marker='s',s=14,facecolors='#CDD4D8',edgecolors='#84939A',linewidths=.45)

        elif g['gene_present']: ax.scatter([gene_x],[y],s=22,c=teal,linewidths=0)

        else: ax.scatter([gene_x],[y],s=22,facecolors='white',edgecolors=teal,linewidths=.9)

        for key,alias,title,x,L in panels:

            s=sites.get(genome_id,{}).get(key)

            if g is None: message='No data available'

            elif not g['gene_present']: message='Gene not detected'

            elif s['status']=='detected': message=None

            elif s['status']=='unresolved_bases': message='Ambiguous base(s)'

            elif s['status']=='indel': message='Indel in target region'

            elif s['status']=='outside_extracted_region': message='Target not present in gene extract'

            else: message='Greater than {} mismatch{}'.format(s['threshold'],'es' if s['threshold']!=1 else '')

            if message:

                ax.text(x,y,message,fontsize=7,va='center',color=muted)

            else:

                mismatch=set(mismatch_positions(oligos[alias],s['sequence']))

                for j,b in enumerate(s['sequence']):

                    if j+1 in mismatch:

                        ax.add_patch(Rectangle((x+j*step-step*.48,y-.36),step*.96,.72,facecolor='#FBE7E8',lw=0))

                    ax.text(x+j*step,y,b,fontfamily='DejaVu Sans Mono',fontsize=7.4,ha='center',va='center',color=red if j+1 in mismatch else ink,fontweight='bold' if j+1 in mismatch else 'normal')

                ax.text(x+L*step+.13,y,str(s['mismatches']),fontsize=7.5,fontweight='bold',ha='center',va='center',color=red if s['mismatches'] else teal)

    hy=n-1

    header_y=hy+.48/row_h

    heading=dict(fontsize=args.header_size,fontweight='bold',color='black',va='baseline')

    ax.text(left,header_y,'Phylogeny',**heading)

    ax.text(label_x,header_y,'Isolate Species and Accession Number',**heading)

    ax.text(gene_x,header_y,'Gene Presence',ha='center',**heading)

    for key,alias,title,x,L in panels:

        heading_text={'fwd':'Forward RPA','rev':'Reverse RPA','gRNA':'gRNA'}[key]

        ax.text(x,header_y,heading_text,**heading)

        for j,b in enumerate(oligos[alias]):

            ax.text(x+j*step,hy+.23/row_h,b,fontfamily='DejaVu Sans Mono',fontsize=7.4,ha='center',color=ink,va='center',fontweight='bold')

        ax.text(x+L*step+.13,hy+.23/row_h,'MM',fontsize=7,ha='center',color=muted,va='center')

        for pos in [1]+list(range(5,L+1,5)):

            ax.text(x+(pos-1)*step,hy+.075/row_h,str(pos),fontsize=5,ha='center',color=muted)

        ax.plot([x-.055,x+L*step+.33],[hy+.13/row_h]*2,color='#777777',lw=.6)

    # Scale bar uses the supplied branch-length units; the tree is midpoint-rooted but never re-ladderized.

    scale=10**math.floor(math.log10(mx/5))

    scale*=max(1,int((mx/5)/scale))

    bar=tree_width*scale/mx;y=-.28/row_h

    ax.plot([left,left+bar],[y,y],color=ink,lw=1)

    ax.text(left,y-.14/row_h,f'{scale:g} '+('edges (branch lengths unavailable)' if cladogram else 'branch-length units'),fontsize=7,ha='left',color=muted)

    y=-.62/row_h

    legend_items=[(left,'o',teal,teal,'Gene detected'),(left+1.85,'o','white',teal,'Gene not detected')]

    has_unmatched_tips=any(t.name not in tip_to_genome for t in tips)

    if has_unmatched_tips:

        legend_items.append((left+4.45,'s','#CDD4D8','#84939A','No data available'))

    for x,marker,fc,ec,text_ in legend_items:

        ax.scatter([x+.04],[y],s=22,marker=marker,facecolors=fc,edgecolors=ec,linewidths=.7)

        ax.text(x+.17,y,text_,fontsize=8,va='center',color=ink)

    mismatch_x=left+(7.2 if has_unmatched_tips else 4.45)

    ax.text(mismatch_x,y,'A',fontsize=9,fontweight='bold',color=red,va='center')

    ax.text(mismatch_x+.19,y,'Mismatch',fontsize=8,va='center',color=ink)

    ax.text(mismatch_x+1.6,y,'MM: mismatch count',fontsize=8,va='center',color=ink)

    fig.savefig(output/'assay_phylogeny_figure.pdf')

    fig.savefig(output/'assay_phylogeny_figure.png',dpi=args.dpi)

    plt.close(fig)



def parser():

    p=argparse.ArgumentParser(description=__doc__,formatter_class=argparse.ArgumentDefaultsHelpFormatter)

    for name in ['genomes','gene-reference','assay-oligos','tree','output-dir']:

        p.add_argument('--'+name,required=True,type=Path)

    p.add_argument('--min-identity',type=float,default=60,help='Minimum BLAST nucleotide identity, percent.')

    p.add_argument('--min-coverage',type=float,default=50,help='Minimum reference coverage across nonoverlapping chained HSPs, percent.')

    p.add_argument('--min-length-factor',type=float,default=.25,help='Minimum extracted length / reference length.')

    p.add_argument('--max-length-factor',type=float,default=1.5,help='Maximum extracted length / reference length.')

    p.add_argument('--max-hsp-gap',type=int,default=1000,help='Maximum unaligned genomic gap between collinear HSPs.')

    p.add_argument('--max-query-gap',type=int,default=200,help='Maximum unaligned reference gap between collinear HSPs.')

    p.add_argument('--max-mismatches',type=int,default=6,help='Default maximum substitutions per full-length ungapped site.')

    p.add_argument('--min-amplicon-len',type=int,default=100,help='Minimum length of a complete RPA product, including both primers.')

    p.add_argument('--max-amplicon-len',type=int,default=250,help='Maximum length of a complete RPA product, including both primers.')

    for key in ['fwd','rev','grna']: p.add_argument('--'+key+'-max-mismatches',type=int)

    p.add_argument('--reference-max-mismatches',type=int,default=0,help='Allowed oligo mismatches when anchoring to reference gene.')

    p.add_argument('--reverse-input',choices=['oligo','target'],default='oligo',help='oligo: preserve supplied reverse sequence; target: reverse-complement it before comparison.')

    p.add_argument('--threads',type=int,help='Total CPU thread budget. If omitted, prompt in an interactive terminal.')

    p.add_argument('--genome-jobs',type=int,help='Maximum concurrent genome searches (default: the full --threads budget). Reduce for large assemblies or limited RAM.')

    p.add_argument('--blastn',default='blastn')

    p.add_argument('--mafft',default='mafft')

    p.add_argument('--evalue',type=float,default=1e-5)

    p.add_argument('--dpi',type=int,default=300)

    p.add_argument('--branch-width',type=float,default=.9,help='Tree branch width in points.')

    p.add_argument('--guide-width',type=float,default=.65,help='Dashed tip-to-label guide width in points.')

    p.add_argument('--header-size',type=float,default=10,help='Uniform black bold column-heading size in points.')

    p.add_argument('--row-height',type=float,default=.20,help='Leaf row spacing in inches.')

    p.add_argument('--show-support',dest='show_support',action='store_true',default=True,help='Show available bootstrap support as percentages (default).')

    p.add_argument('--hide-support',dest='show_support',action='store_false',help='Hide branch support labels.')

    p.add_argument('--force',action='store_true',help='Replace the output directory only if it belongs to this workflow.')

    return p



def main(argv=None):

    args=parser().parse_args(argv)

    args.threads=resolve_threads(args.threads)

    if args.genome_jobs is not None and args.genome_jobs < 1: raise ValueError('--genome-jobs must be a positive integer.')

    if not 0<args.min_identity<=100 or not 0<args.min_coverage<=100: raise ValueError('Identity and coverage must be percentages in (0,100].')

    if not 0<args.min_length_factor<=args.max_length_factor: raise ValueError('Invalid length-factor limits.')

    for name in ['max_mismatches','reference_max_mismatches','fwd_max_mismatches','rev_max_mismatches','grna_max_mismatches']:

        value=getattr(args,name)

        if value is not None and value<0: raise ValueError('Mismatch limits must be nonnegative.')

    if args.max_hsp_gap<0 or args.max_query_gap<0: raise ValueError('HSP gap limits must be nonnegative.')

    if not 1<=args.min_amplicon_len<=args.max_amplicon_len: raise ValueError('Invalid amplicon length limits.')

    if min(args.branch_width,args.guide_width,args.header_size)<=0: raise ValueError('Plot line widths and heading size must be positive.')

    if args.threads<1 or args.dpi<72 or args.row_height<.13 or args.evalue<=0: raise ValueError('Invalid threads, dpi, row height or evalue.')

    log('[1/6] Validating inputs and mapping oligos to the reference gene')

    if not args.genomes.is_dir(): raise ValueError('--genomes must be a directory containing unpacked FASTA assemblies.')

    files=sorted(p for p in args.genomes.rglob('*') if p.is_file() and re.search(r'\.(fna|fa|fasta)(\.gz)?$',p.name,re.I))

    if not files: raise ValueError('No .fna, .fa, .fasta or corresponding .gz genomes found.')

    names=[p.name for p in files]

    if len(names)!=len(set(names)): raise ValueError('Duplicate genome basenames; identifiers must be unique across subdirectories.')

    if any(re.search(r'\s',n) for n in names): raise ValueError('Genome filenames must not contain whitespace.')

    for executable in [args.blastn,args.mafft]:

        if not shutil.which(executable): raise ValueError(f'Missing executable: {executable}. Install dependencies as described in README.md.')

    references=read_fasta(args.gene_reference)

    if len(references)!=1: raise ValueError('Reference gene FASTA must contain exactly one record.')

    reference=next(iter(references.values()))

    if set(reference)-set('ACGT'): raise ValueError('Reference gene must contain only A/C/G/T for unambiguous site anchoring.')

    raw=read_fasta(args.assay_oligos)

    aliases={'fwdprimer':'fwdPrimer','forwardrpaprimer':'fwdPrimer','revprimer':'revPrimer','reverserpaprimer':'revPrimer','grna':'gRNA'}

    oligos={}

    for name,seq in raw.items():

        canon=aliases.get(name.lower())

        if canon is None or canon in oligos: raise ValueError('Oligo IDs must uniquely identify fwdPrimer, revPrimer and gRNA (ForwardRPAprimer/ReverseRPAprimer also accepted).')

        oligos[canon]=seq

    if set(oligos)!={'fwdPrimer','revPrimer','gRNA'}: raise ValueError('Supply all three oligos.')

    original=oligos.copy()

    if args.reverse_input=='target': oligos['revPrimer']=rc(oligos['revPrimer'])

    anchors={}

    for key,alias,_,_ in COMPONENTS:

        try: anchors[key]=locate_oligo(reference,oligos[alias],args.reference_max_mismatches)

        except ValueError as e: raise ValueError(f'{alias}: {e}') from e

        a=anchors[key];log(f'  {alias}: {len(oligos[alias])} nt; reference {a["start"]+1}-{a["end"]} ({a["strand"]})')

    warnings=[]

    def warn(msg): warnings.append(msg);log('WARNING: '+msg)

    _,reference_assessment,_=evaluate_candidate(
        dict(start=1,end=len(reference),strand='+',bitscore=0,coverage=100,identity=100),
        reference,reference,oligos,anchors,args)
    if not reference_assessment['assay_compatible']:
        warn('The reference assay is not complete under the configured primer orientation, guide placement, '
             f'amplicon length and mismatch criteria ({reference_assessment["site_geometry"]}); '
             'check oligo orientation and --min-amplicon-len/--max-amplicon-len.')

    if anchors['fwd']['strand']==anchors['rev']['strand']:

        warn('Both primer inputs map to the same reference strand. Sequences are preserved as supplied. If the reverse entry is the forward-strand target, rerun with --reverse-input target. This run is not evidence of a valid primer pair.')

    tree=inspect_tree(args.tree)

    out=args.output_dir.resolve()

    # Do not let --force remove supplied data or a broader working directory.

    inputs=[args.genomes.resolve(),args.tree.resolve(),args.gene_reference.resolve(),args.assay_oligos.resolve()]

    if any(x==out or out in x.parents for x in inputs): raise ValueError('Output directory must not contain any input path.')

    if out.exists() and any(out.iterdir()):

        if not args.force or not (out/'.assay_phylogeny_output').is_file(): raise ValueError('Output directory is nonempty; choose a new directory or use --force for a previous workflow output.')

        shutil.rmtree(out)

    out.mkdir(parents=True,exist_ok=True);(out/'.assay_phylogeny_output').write_text(VERSION)

    work=out/'work';work.mkdir();logs=out/'logs';logs.mkdir()

    fasta(work/'reference.fasta',{REF:reference})

    fasta(out/'effective_assay_oligos.fasta',oligos)

    table(out/'oligo_reference_mapping.csv',[dict(component=key,input_sequence=original[alias],effective_sequence=oligos[alias],reference_start=anchors[key]['start']+1,reference_end=anchors[key]['end'],reference_strand=anchors[key]['strand'],reference_mismatches=anchors[key]['reference_mismatches']) for key,alias,_,_ in COMPONENTS],['component','input_sequence','effective_sequence','reference_start','reference_end','reference_strand','reference_mismatches'])

    genes={};candidate_hits={};candidates_by_genome={};extracted={}

    log(f'[2/6] Searching {len(files)} genome(s) with BLASTn; retaining all qualifying hits')

    # Each BLAST -subject search uses a single core, so by default every genome
    # runs concurrently up to your full --threads budget. Pass --genome-jobs to
    # cap this lower if peak memory becomes a problem on very large assemblies.
    workers=min(len(files),args.threads,args.genome_jobs or args.threads)
    log(f'  {workers} concurrent genome search(es); MAFFT: {args.threads} thread(s)')

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures=[pool.submit(search_genome,index,path,work/'reference.fasta',len(reference),work,logs,args)
                 for index,path in enumerate(files,1)]
        for index,(path,future) in enumerate(zip(files,futures),1):
            try: row,candidates=future.result()
            except Exception as exc: raise RuntimeError(f'Genome search failed for {path.name}: {exc}') from exc
            genes[path.name]=row
            candidates_by_genome[path.name]=candidates
            candidate_hits.update({hit['candidate_id']:hit for hit in candidates})
            log(f'  [{index}/{len(files)}] {path.name}: {len(candidates)} qualifying candidate(s)')

    log('[3/6] Aligning every candidate hit with the gene reference')

    alignment_input={REF:reference,**{cid:h['sequence'] for cid,h in candidate_hits.items()}}
    fasta(work/'alignment_input.fasta',alignment_input)
    fasta(out/'all_candidate_loci.fasta',{cid:h['sequence'] for cid,h in candidate_hits.items()})
    all_aligned_path=out/'all_candidate_loci_aligned.fasta'

    if candidate_hits:
        run([args.mafft,'--auto','--thread',str(args.threads),'--inputorder',work/'alignment_input.fasta'],
            all_aligned_path,logs/'mafft.log')
        aligned=read_fasta(all_aligned_path,aligned=True)
    else:
        aligned=alignment_input
        fasta(all_aligned_path,aligned)
        log('  No qualifying gene hits; continuing with empty gene and assay-site outputs.')

    if set(aligned)!=set(alignment_input): raise ValueError('MAFFT changed or lost candidate identifiers.')
    for name,seq in aligned.items():
        if seq.replace('-','')!=alignment_input[name]: raise ValueError(f'MAFFT sequence integrity check failed: {name}')

    log('[4/6] Comparing all candidate loci and selecting the best assay match per genome')

    evaluations={};candidate_rows=[];selected_ids={}
    for cid,hit in candidate_hits.items():
        sites,assessment,rank=evaluate_candidate(hit,aligned[REF],aligned[cid],oligos,anchors,args)
        evaluations[cid]=dict(sites=sites,assessment=assessment,rank=rank)

    for name,candidates in candidates_by_genome.items():
        if not candidates: continue
        best=max(candidates,key=lambda h:evaluations[h['candidate_id']]['rank'])
        selected_ids[name]=best['candidate_id']
        extracted[name]=best['sequence']
        genes[name].update({k:v for k,v in best.items() if k not in ('sequence','candidate_id','genome_id')})
        genes[name]['selected_candidate_id']=best['candidate_id']
        genes[name]['assay_compatible_candidate_count']=sum(
            evaluations[h['candidate_id']]['assessment']['assay_compatible'] for h in candidates)
        chosen=evaluations[best['candidate_id']]['assessment']
        log(f'  {name}: selected {best["candidate_id"]} on {best["contig"]}:{best["start"]}-{best["end"]}; '
            f'{chosen["sites_detected"]}/3 sites within mismatch limits; geometry={chosen["site_geometry"]}')

    for cid,hit in candidate_hits.items():
        evaluation=evaluations[cid];row={k:v for k,v in hit.items() if k!='sequence'}
        row.update(evaluation['assessment'])
        row['selected']=int(selected_ids.get(hit['genome_id'])==cid)
        for key,_,_,_ in COMPONENTS:
            site=evaluation['sites'][key]
            for field in ('status','present','candidate_mismatches','mismatch_positions','candidate_sequence',
                          'genomic_start','genomic_end','genomic_strand'):
                row[f'{key}_{field}']=site[field]
        candidate_rows.append(row)

    candidate_fields=['genome_id','candidate_id','selected','contig','start','end','strand','identity','coverage',
                      'bitscore','evalue','extracted_length','hsp_count','assay_compatible','site_geometry',
                      'amplicon_length','sites_detected','total_mismatches']+[f'{key}_{field}'
                      for key,_,_,_ in COMPONENTS for field in ('status','present','candidate_mismatches',
                      'mismatch_positions','candidate_sequence','genomic_start','genomic_end','genomic_strand')]
    table(out/'all_candidate_loci.csv',candidate_rows,candidate_fields)

    fasta(out/'extracted_genes.fasta',extracted)
    gene_fields=['genome_id','gene_present','identity','coverage','contig','start','end','strand','bitscore',
                 'evalue','query_start','query_end','extracted_length','hsp_count','valid_candidate_count',
                 'selected_candidate_id','assay_compatible_candidate_count']
    table(out/'gene_presence.csv',genes.values(),gene_fields)

    selected_alignment={REF:aligned[REF],**{name:aligned[cid] for name,cid in selected_ids.items()}}
    fasta(out/'reference_anchored_alignment.fasta',selected_alignment)
    fasta(out/'extracted_genes_aligned.fasta',{name:aligned[cid] for name,cid in selected_ids.items()})

    sites={};summary=[];detail=[]
    for name,gene in genes.items():
        sites[name]={};row=dict(genome_id=name,gene_present=gene['gene_present'],
                                selected_candidate_id=gene['selected_candidate_id'])
        for key,alias,_,_ in COMPONENTS:
            threshold=getattr(args,key.lower()+'_max_mismatches')
            if threshold is None: threshold=args.max_mismatches
            if name in selected_ids:
                site=evaluations[selected_ids[name]]['sites'][key]
            else:
                site=dict(present=0,status='gene_not_detected',mismatches='',mismatch_positions='',sequence='',
                          candidate_sequence='',candidate_mismatches='',insertions=0,deletions=0,start='',end='',
                          strand=anchors[key]['strand'],genomic_start='',genomic_end='',genomic_strand='',
                          threshold=threshold)
            sites[name][key]=site
            for field in ('present','mismatches','status','mismatch_positions'):
                row[f'{key}_{field}']=site[field]
            detail.append(dict(genome_id=name,component=key,max_mismatches=threshold,**site))
        summary.append(row)

    fields=['genome_id','gene_present','selected_candidate_id']+[f'{key}_{field}'
            for key,_,_,_ in COMPONENTS for field in ('present','mismatches','status','mismatch_positions')]
    table(out/'assay_site_summary.csv',summary,fields)
    table(out/'assay_site_details.csv',detail,['genome_id','component','max_mismatches','present','status',
          'mismatches','mismatch_positions','sequence','candidate_sequence','candidate_mismatches',
          'insertions','deletions','start','end','strand','genomic_start','genomic_end','genomic_strand'])

    for key,alias,_,filename in COMPONENTS:

        rec={f'REFERENCE_{alias}':oligos[alias]}

        rec.update({name:sites[name][key]['sequence'] for name in genes if sites[name][key]['present']==1})

        if len({len(s) for s in rec.values()})!=1: raise ValueError(f'Inconsistent short-site lengths for {key}')

        fasta(out/filename,rec)

    log('[5/6] Reconciling complete filenames with Newick tip labels')

    tipnames=[t.name for t in tree.get_terminals()]

    tip_to_genome,match_methods,missing,extra,ambiguous=reconcile_identifiers(tipnames,list(genes))

    sections=[('Tree tips with no matching gene-extraction record',missing),('Gene-extraction records with no matching tree tip',extra)]

    with (out/'identifier_mismatch_report.txt').open('w') as handle:

        handle.write('Matching: exact filename, then recognized FASTA extension and whitespace differences, then unique species/accession when only one side lacks a version. Two explicit different versions never match.\n\n')

        handle.write('Matched tree tip -> original genome filename (method)\n')

        for tip in tipnames:

            if tip in tip_to_genome:

                handle.write(f'{tip} -> {tip_to_genome[tip]} ({match_methods[tip]})\n')

        handle.write('\n')

        for title,items in sections:

            handle.write(title+f' ({len(items)})\n'+'\n'.join(items or ['(none)'])+'\n\n')

            if items:

                preview=', '.join(items[:5])

                suffix=f', ... ({len(items)-5} more in report)' if len(items)>5 else ''

                warn(f'{title} ({len(items)}): {preview}{suffix}')

        handle.write(f'Ambiguous unmatched tree tips ({len(ambiguous)})\n')

        for method,tip,candidates in ambiguous:

            handle.write(f'{tip} ({method}) -> {", ".join(candidates)}\n')

    if ambiguous: warn(f'{len(ambiguous)} ambiguous tip match(es) were left unmatched; see identifier_mismatch_report.txt')

    normalized=sum(method!='exact' for method in match_methods.values())

    if normalized: log(f'  Matched {normalized} tip(s) after conservative label normalization.')

    table(out/'plotted_tip_order.csv',[dict(row=i+1,tree_tip=t,genome_id=tip_to_genome.get(t,''),match_method=match_methods.get(t,''),has_genome_result=int(t in tip_to_genome)) for i,t in enumerate(tipnames)],['row','tree_tip','genome_id','match_method','has_genome_result'])

    log('[6/6] Drawing the shared-row phylogeny and assay panels; exporting PDF and PNG')

    plot_figure(tree,genes,sites,oligos,out,args,tip_to_genome)

    stats=dict(genomes_searched=len(genes),genes_detected=len(extracted),
               candidate_hits_evaluated=len(candidate_hits),
               assay_compatible_candidates=sum(e['assessment']['assay_compatible'] for e in evaluations.values()),
               genomes_with_complete_assay=sum(g['assay_compatible_candidate_count']>0 for g in genes.values()),
               tree_tips_plotted=len(tipnames),tree_tips_with_results=len(tip_to_genome),
               tree_tips_without_results=len(missing),records_without_tree_tip=len(extra),
               ambiguous_tip_matches=len(ambiguous),
               sites_detected={key:sum(s[key]['present']==1 for s in sites.values()) for key,_,_,_ in COMPONENTS})

    versions={}

    for key,cmd in [('blastn',[args.blastn,'-version']),('mafft',[args.mafft,'--version'])]:

        r=subprocess.run(cmd,capture_output=True,text=True);versions[key]=(r.stdout+r.stderr).strip()

    import Bio,matplotlib

    versions.update(python=sys.version,biopython=Bio.__version__,matplotlib=matplotlib.__version__)

    hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [args.gene_reference,args.assay_oligos,args.tree]}

    manifest=dict(workflow_version=VERSION,utc=datetime.now(timezone.utc).isoformat(),command=sys.argv,parameters={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()},versions=versions,input_sha256=hashes,genomes=[dict(path=str(p.resolve()),size_bytes=p.stat().st_size) for p in files],summary=stats,warnings=warnings)

    (out/'run_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')

    (out/'run_summary.json').write_text(json.dumps(stats,indent=2)+'\n')

    (out/'warnings.txt').write_text('\n'.join(warnings)+'\n')

    log('COMPLETE\n'+json.dumps(stats,indent=2)+'\nOutputs: '+str(out))

    return stats



if __name__=='__main__':

    try: main()

    except (ValueError,RuntimeError,OSError) as exc:

        print(f'ERROR: {exc}',file=sys.stderr);sys.exit(1)
