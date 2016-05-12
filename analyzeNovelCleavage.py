import argparse, os, sys
from Kleat import *
from GTF import *

class Novel:
    
    @staticmethod
    def strToInt(string):
        try:
            return int(string)
        except ValueError:
            return None
    
    def __init__(self, sample_id=None, cell_line_name=None,
                 tissue_site=None, disease=None, chrom=None,
                 gene=None, cleavage_site=None, 
                 gene_start=None, gene_end=None, 
                 transcript_id=None, strand=None, 
                 distance=None):
        self.sample_id = sample_id
        self.cell_line_name = cell_line_name
        self.tissue_site = tissue_site
        self.disease = disease
        self.chrom = chrom
        self.gene = gene
        self.cleavage_site = Novel.strToInt(cleavage_site)
        self.gene_start = Novel.strToInt(gene_start)
        self.gene_end = Novel.strToInt(gene_end)
        self.transcript_id = transcript_id
        self.strand = strand
        self.distance = Novel.strToInt(distance)

    def __str__(self):
        return ('\t').join([str(x) for x in [self.sample_id, self.cell_line, self.tissue, self.disease, self.chrom, self.gene, self.cleavage_site, self.gene_start, self.gene_end, self.transcript_id, self.strand, self.distance]])

def sprint(text):
    sys.stdout.write(text)
    sys.stdout.flush()

def parseNovel(novel):
    results = {}
    with open(novel, 'r') as f:
        header = f.readline()
        for i,line in enumerate(f):
            novel = Novel(*line.strip().split('\t'))
            if novel.tissue_site not in results:
                results[novel.tissue_site] = {}
            if novel.chrom not in results[novel.tissue_site]:
                results[novel.tissue_site][novel.chrom] = {}
            if novel.gene not in results[novel.tissue_site][novel.chrom]:
                results[novel.tissue_site][novel.chrom][novel.gene] = []
            results[novel.tissue_site][novel.chrom][novel.gene].append(novel)
    return results

def analyzeNovel(novel, distance=20):
    tissues = novel.keys()
    n = len(tissues)
    for i in xrange(n-1):
        for j in xrange(i+1,n):
            #print 'Comparing {} to {}'.format(tissues[i], tissues[j])
            for chrom in novel[tissues[i]]:
                if chrom not in novel[tissues[j]]:
                    continue
                for gene in novel[tissues[i]][chrom]:
                    if gene not in novel[tissues[j]][chrom]:
                        continue
                    for site in novel[tissues[i]][chrom][gene]:
                        match = False
                        for site2 in novel[tissues[j]][chrom][gene]:
                            if abs(site2.cleavage_site - site.cleavage_site) <= distance:
                                match = True
                                break
                        if not match:
                            print '{}\t{}\t{}\t{}'.format(site.cell_line_name, site.sample_id, site.tissue_site, site.cleavage_site)

parser = argparse.ArgumentParser(description='')

parser.add_argument('novel_sites', help='novelCleavage.py output file')
parser.add_argument('-n', '--name', default='novel.sites.analysis', help='Name for the file output. Default is "novel.sites"')
parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Path to output to. Default is "{}"'.format(os.getcwd()))

args = parser.parse_args()

if not os.path.isdir(args.outdir):
    os.makedirs(args.outdir)

#with open(args.novel_sites) as f:
#    N = len(f.readlines()) - 1

sprint('Parsing novel sites ...')
novels = parseNovel(args.novel_sites)
print 'DONE'

analysis = analyzeNovel(novels)
