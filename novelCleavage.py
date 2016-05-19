import argparse, os, sys
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

def parseSimpleKleat(kleat,
                     max_dist_ann = float('inf'), 
                     min_len_tail_contig = None,
                     min_num_tail_reads = None,
                     min_num_bridge_reads = None,
                     min_bridge_read_tail_len = None,
                     has_pas = None):
    results = []
    with open(kleat, 'r') as f:
        header = f.readline()
        for line in f:
            result = Kleat(*line.strip().split('\t'))
            if (result.length_of_tail_in_contig < min_len_tail_contig):
                if (result.number_of_bridge_reads < min_num_bridge_reads):
                    continue
            elif min_num_tail_reads and (0 <= result.number_of_tail_reads < min_num_tail_reads):
                continue
            elif min_bridge_read_tail_len and result.number_of_bridge_reads > 0 and (0 <= result.max_bridge_read_tail_length < min_bridge_read_tail_len):
                continue
            elif (result.distance_from_annotated_site > max_dist_ann):
                continue
            elif has_pas and not result.pas:
                continue
            result.ID = os.path.basename(kleat).split('.')[0]
            results.append(result)
    return results

def parseKleat(kleat, 
               max_dist_ann = float('inf'), 
               min_len_tail_contig = None,
               min_num_tail_reads = None,
               min_num_bridge_reads = None,
               min_bridge_read_tail_len = None,
               has_pas = None):
    results = []
    with open(kleat, 'r') as f:
        header = f.readline()
        for line in f:
            result = Kleat(*line.strip().split('\t'))
            if min_len_tail_contig and (0 <= result.length_of_tail_in_contig < min_len_tail_contig):
                continue
            elif min_num_tail_reads and (0 <= result.number_of_tail_reads < min_num_tail_reads):
                continue
            elif min_num_bridge_reads and (0 <= result.number_of_bridge_reads < min_num_bridge_reads):
                continue
            elif min_bridge_read_tail_len and (0 <= result.max_bridge_read_tail_length < min_bridge_read_tail_len):
                continue
            elif (result.distance_from_annotated_site > max_dist_ann):
                continue
            elif has_pas and not result.pas:
                continue
            result.ID = os.path.basename(kleat).split('.')[0]
            results.append(result)
    return results

def groupKleat(parsed):
    results = {}
    for r in parsed:
        if r.chromosome not in results:
            results[r.chromosome] = {r.gene: [r]}
        if r.gene not in results[r.chromosome]:
            results[r.chromosome][r.gene] = [r]
        else:
            results[r.chromosome][r.gene].append(r)
    return results

def findNovelCleavageEvents(kleats_grouped, gtf_grouped, distance=20):
    results = []
    kg = kleats_grouped
    gg = gtf_grouped
    for chrom in kg:
        if chrom not in gg:
            continue
        for gene in kg[chrom]:
            if gene not in gg[chrom]:
                continue
            strand = gg[chrom][gene][0].strand
            sites = [x.cleavage_site for x in kg[chrom][gene]]
            utrs = gg[chrom][gene]
            for i,site in enumerate(sites):
                my_min = float('inf')
                closest_utr = None
                for utr in utrs:
                    current_min = min([abs(site - utr.start), abs(utr.end - site)])
                    if current_min <= my_min:
                        my_min = current_min
                        closest_utr = utr
                if my_min > distance:
                    results.append([kg[chrom][gene][i], closest_utr, my_min])
    return results

def groupEnsemblGTF(gtf):
    result = {}
    for g in gtf:
        if g.seqname not in result:
            result[g.seqname] = {g.attribute['gene_name']: [g]}
        if g.attribute['gene_name'] not in result[g.seqname]:
            result[g.seqname][g.attribute['gene_name']] = [g]
        else:
            result[g.seqname][g.attribute['gene_name']].append(g)
    #for chrom in result:
    #    for gene in result[chrom]:
    #        result[chrom][gene].sort(key=lambda x: x.start)
    return result

def groupAceviewGTF(gtf, result=None):
    if not result:
        result = {}
    for g in gtf:
        if g.seqname not in result:
            result[g.seqname] = {g.attribute['gene_id']: [g]}
        if g.attribute['gene_id'] not in result[g.seqname]:
            result[g.seqname][g.attribute['gene_id']] = [g]
        else:
            result[g.seqname][g.attribute['gene_id']].append(g)
    for chrom in result:
        for gene in result[chrom]:
            result[chrom][gene].sort(key=lambda x: x.start)
    return result

parser = argparse.ArgumentParser(description='')

parser.add_argument('kleat', nargs='+', help='Kleat files')
parser.add_argument('ensembl', help='GTF annotation file for ensembl. Must have "UTR" feature for every gene')
parser.add_argument('aceview', help='GTF annotation file for aceview. Must have "UTR" feature for every gene')
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
annots = []

# Parse KLEAT data
N = len(args.kleat)
for i,kleat in enumerate(args.kleat):
    sprint('Loading KLEAT {}/{} ...\r'.format(i+1,N))
#    kleats += parseKleat(kleat,
#                         args.max_dist_ann,
#                         args.min_len_tail_contig,
#                         args.min_num_tail_reads,
#                         args.min_num_bridge_reads,
#                         has_pas = args.has_polyadenylation_signal)
    kleats += parseSimpleKleat(kleat,
                         args.max_dist_ann,
                         args.min_len_tail_contig,
                         args.min_num_tail_reads,
                         args.min_num_bridge_reads,
                         has_pas = args.has_polyadenylation_signal)
print 'Loading KLEAT {}/{} ...DONE\r'.format(N,N)

# Group KLEAT data
sprint('Grouping kleat data ...')
kleats = groupKleat(kleats)
print 'DONE'

# Parse ensembl
sprint('Parsing ensembl annotation ...')
ensembl = parseGTF(args.ensembl, seqnames=set(['chr' + str(x) for x in range(1,24)] + ['X', 'Y']), sources=set(['protein_coding','retained_intron']), features=set(['transcript']))
print 'DONE'

# Parse aceview
sprint('Parsing aceview annotation ...')
aceview = parseGTF(args.aceview, seqnames=set(['chr' + str(x) for x in range(1,24)] + ['X', 'Y']))
print 'DONE'

# Group ensembl
sprint('Grouping ensembl data ...')
gensembl = groupEnsemblGTF(ensembl)
print 'DONE'

# Group aceview
sprint('Grouping aceview data ...')
ggtf = groupAceviewGTF(aceview, gensembl)
print 'DONE'

# Parse summary file
sprint('Parsing summary file ...')
summary = parseSummary(args.summary)
print 'DONE'

# Find novel sites
sprint('Finding novel sites ...')
sites = findNovelCleavageEvents(kleats, ggtf, distance=args.min_novel_distance)
print 'DONE'

outfile = os.path.join(args.outdir, args.name)

sprint('Writing to: {} ...'.format(outfile))
with open(outfile, 'w') as o:
    o.write(('\t').join(['ID','CELL_LINE','TISSUE','DISEASE','CHROM','GENE','CLEAVAGE_SITE','GENE_START','GENE_END','TRANSCRIPT_ID','STRAND','DISTANCE']) + '\n')
    for clv_evnt, utr, my_min in sites:
        o.write(('\t').join([str(x) for x in [clv_evnt.ID, summary[clv_evnt.ID][0], summary[clv_evnt.ID][1], summary[clv_evnt.ID][2], utr.seqname, clv_evnt.gene, clv_evnt.cleavage_site, utr.start, utr.end, utr.attribute['transcript_id'] ,utr.strand, my_min]]) + '\n')
print 'DONE'
