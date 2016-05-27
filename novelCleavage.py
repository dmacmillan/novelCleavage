import pysam
import argparse, os, sys
from bisect import bisect_left
from Kleat import *
from GTF import *

def sprint(text):
    sys.stdout.write(text)
    sys.stdout.flush()

def parseSummary(summary):
    results = {}
    with open(summary, 'r') as f:
        header = f.readline()
        for line in f:
            line = line.strip().split('\t')
            # cell_line_name, tissue_site, disease
            results[line[11]] = [line[0], line[2], line[7]]
    return results

def minDistanceBinarySearch(a, x, lo=0, hi=None):   # can't use a to specify default for hi
    hi = hi if hi is not None else len(a) # hi defaults to len(a)   
    pos = bisect_left(a,x,lo,hi)          # find insertion position
    #return (pos if pos != hi and a[pos] == x else -1) # don't walk off the end
    if pos != hi and a[pos] == x:
        return pos
    elif pos != 0:
        return pos-1
    else:
        return 0

def groupBy(group, att):
	for c in group:
		for g in group[c]:
			strand = group[c][g][0].strand
			transcripts = {}
			for x in group[c][g]:
				tid = x.attribute[att]
				if tid not in transcripts:
					transcripts[tid] = [x]
				else:
					transcripts[tid].append(x)
			group[c][g] = transcripts

def findNovelCleavageEvents(kleats_grouped, all_annotations, window=20, cutoff=5000):
    results = []
    kg = kleats_grouped
    aa = all_annotations
    for chrom in kg:
        for gene in kg[chrom]:
            sites = [x.cleavage_site for x in kg[chrom][gene]]
            strand = kg[chrom][gene][0].transcript_strand
            for i,site in enumerate(sites):
                potentials = []
                for ann in all_annotations:
                    try:
                        potentials += [g for g in ann.fetch(kg[chrom][gene][i].chromosome, site-window, site+window) if g.strand == strand and g.feature == 'exon']
                    except ValueError:
                        continue
                if not potentials:
#                    print '^'*20
#                    print 'Cleavage event:'
#                    print kg[chrom][gene][i]
#                    print 'i:'
#                    print i
#                    print 'Site:'
#                    print site
                    # find closest site
                    add = window
                    while not potentials and add <= cutoff:
                        try:
                            potentials += [g for g in ann.fetch(kg[chrom][gene][i].chromosome, site-window-add, site+window+add) if g.strand == strand and g.feature == 'exon']
                            add += window
                        except ValueError:
                            continue
                    if not potentials:
                        results.append([kg[chrom][gene][i], None, None])
                        continue
                    if strand == '+':
                        my_min, idx = min([abs(site - val.end), idx] for [idx, val] in enumerate(potentials))
                    else:
                        my_min, idx = min([abs(site - val.start), idx] for [idx, val] in enumerate(potentials))
                    closest_utr = potentials[idx]
#                    print 'Potentials:'
#                    for x in potentials:
#                        print x
#                    print 'Closest UTR:'
#                    print closest_utr
#                    print 'my_min:'
#                    print my_min
                    #if my_min > window:
                    results.append([kg[chrom][gene][i], closest_utr, my_min])
    return results

parser = argparse.ArgumentParser(description='')

parser.add_argument('kleat', nargs='+', help='Kleat files')
parser.add_argument('ensembl', help='GTF annotation file for ensembl. Must have "UTR" feature for every gene')
parser.add_argument('aceview', help='GTF annotation file for aceview. Must have "UTR" feature for every gene')
parser.add_argument('refseq', help='GTF annotation file for refseq')
parser.add_argument('ucsc', help='GTF annotation file for ucsc')
parser.add_argument('summary', help='Summary file for CCLE data')
parser.add_argument('-x', '--min_novel_distance', type=int, default=20, help='The minimum distance between a cleavage site and an annotated transcript for the site to be deemed novel. Default = 20')
parser.add_argument('-d', '--max_dist_ann', type=int, default=float('inf'), help='Maximum distance from annotated site. Default = infinity')
parser.add_argument('-l', '--min_len_tail_contig', type=int, default=0, help='Minimum length of tail in contig. Default = 0')
parser.add_argument('-t', '--min_num_tail_reads', type=int, default=0, help='Minimum number of tail reads. Default = 0')
parser.add_argument('-b', '--min_num_bridge_reads', type=int, default=0, help='Minimum number of bridge reads. Default = 0')
parser.add_argument('-bl', '--min_bridge_read_tail_len', type=int, default=0, help='Minimum length of each bridge read. Default = 0')
parser.add_argument('-p', '--has_polyadenylation_signal', action='store_true', help='Filter out cleavage events that do not have a polyadenylation signal')
#parser.add_argument('',help='')
parser.add_argument('-n', '--name', default='novel.sites', help='Name for the file output. Default is "novel.sites"')
parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Path to output to. Default is "{}"'.format(os.getcwd()))

args = parser.parse_args()

if not os.path.isdir(args.outdir):
    os.makedirs(args.outdir)

kleats = []

# Parse KLEAT data
N = len(args.kleat)
for i,kleat in enumerate(args.kleat):
    sprint('Loading KLEAT {}/{} ...\r'.format(i+1,N))
    kleats += Kleat.parseKleat(kleat,
                               args.max_dist_ann,
                               args.min_len_tail_contig,
                               args.min_num_tail_reads,
                               args.min_num_bridge_reads,
                               args.min_bridge_read_tail_len,
                               has_pas = args.has_polyadenylation_signal)
print 'Loading KLEAT {}/{} ...DONE\r'.format(N,N)

# Group KLEAT data
sprint('Grouping kleat data ...')
kleats = Kleat.groupKleat(kleats)
print 'DONE'

# Parse ensembl
sprint('Loading ensembl annotation ...')
ensembl = pysam.TabixFile(args.ensembl, parser=pysam.asGTF())
print 'DONE'

# Parse aceview
sprint('Loading aceview annotation ...')
aceview = pysam.TabixFile(args.aceview, parser=pysam.asGTF())
print 'DONE'

# Parse refseq
sprint('Loading refseq annotation ...')
refseq = pysam.TabixFile(args.refseq, parser=pysam.asGTF())
print 'DONE'

# Parse ucsc
sprint('Loading ucsc annotation ...')
ucsc = pysam.TabixFile(args.ucsc, parser=pysam.asGTF())
print 'DONE'

all_annotations = [ensembl, aceview, refseq, ucsc]

# Parse summary file
sprint('Parsing summary file ...')
summary = parseSummary(args.summary)
print 'DONE'

# Find novel sites
sprint('Finding novel sites ...')
sites = findNovelCleavageEvents(kleats, all_annotations, window=args.min_novel_distance)
print 'DONE'

outfile = os.path.join(args.outdir, args.name)

sprint('Writing to: {} ...'.format(outfile))
with open(outfile, 'w') as o:
    o.write(('\t').join(['ID','CELL_LINE','TISSUE','DISEASE','CHROM','GENE','CLEAVAGE_SITE','GENE_START','GENE_END','TRANSCRIPT_ID','STRAND','DISTANCE']) + '\n')
    for clv_evnt, utr, my_min in sites:
        if not utr:
            o.write(('\t').join([str(x) for x in [clv_evnt.ID, summary[clv_evnt.ID][0], summary[clv_evnt.ID][1], summary[clv_evnt.ID][2], 'NA', clv_evnt.gene, clv_evnt.cleavage_site, 'NA', 'NA', 'NA', 'NA', 'NA']]) + '\n')
        else:
            o.write(('\t').join([str(x) for x in [clv_evnt.ID, summary[clv_evnt.ID][0], summary[clv_evnt.ID][1], summary[clv_evnt.ID][2], utr.contig, clv_evnt.gene, clv_evnt.cleavage_site, utr.start, utr.end, utr.asDict()['transcript_id'], utr.strand, my_min]]) + '\n')
print 'DONE'
