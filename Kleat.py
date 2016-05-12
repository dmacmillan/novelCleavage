class Kleat:
    
    @staticmethod
    def strToInt(string):
        try:
            return int(string)
        except ValueError:
            return 0

    def __init__(self, gene, transcript, transcript_strand, coding, contig, chromosome, cleavage_site, within_UTR, distance_from_annotated_site, ESTs, length_of_tail_in_contig, number_of_tail_reads, number_of_bridge_reads, max_bridge_read_tail_length, bridge_read_identities, tail_and_bridge_reads, number_of_link_pairs, max_link_pair_length, link_pair_identities, pas, utr3):
        self.gene = gene
        self.transcript = transcript
        self.transcript_strand = transcript_strand
        self.coding = coding
        self.contig = contig
        self.chromosome = chromosome
        self.cleavage_site = int(cleavage_site)
        if (within_UTR == 'no'):
            self.within_UTR = False
        else:
            self.within_UTR = True
        self.distance_from_annotated_site = Kleat.strToInt(distance_from_annotated_site)
        self.ESTs = ESTs
        self.length_of_tail_in_contig = Kleat.strToInt(length_of_tail_in_contig)
        self.number_of_tail_reads = Kleat.strToInt(number_of_tail_reads)
        self.number_of_bridge_reads = Kleat.strToInt(number_of_bridge_reads)
        self.max_bridge_read_tail_length = Kleat.strToInt(max_bridge_read_tail_length)
        self.bridge_read_identities = bridge_read_identities
        self.tail_and_bridge_reads = Kleat.strToInt(tail_and_bridge_reads)
        self.number_of_link_pairs = Kleat.strToInt(number_of_link_pairs)
        self.max_link_pair_length = Kleat.strToInt(max_link_pair_length)
        self.link_pair_identities = link_pair_identities
        self.pas = self.utr3 = None
        try:
            self.pas = [int(x) for x in pas.split(':')]
        except ValueError:
            pass
        try:
            self.utr3 = [int(x) for x in utr3.split('-')]
        except ValueError:
            pass

    def __str__(self):
        if self.pas:
            self.pas = (':').join([str(x) for x in self.pas])
        atts = [self.gene, self.transcript, self.transcript_strand, self.coding, self.contig, self.chromosome, self.cleavage_site, self.within_UTR, self.distance_from_annotated_site, self.ESTs, self.length_of_tail_in_contig, self.number_of_tail_reads, self.number_of_bridge_reads, self.max_bridge_read_tail_length, self.bridge_read_identities, self.tail_and_bridge_reads, self.number_of_link_pairs, self.max_link_pair_length, self.link_pair_identities, self.pas, self.utr3]
        atts = [str(x) if x is not None else '-' for x in atts]
        return ('\t').join(atts)
