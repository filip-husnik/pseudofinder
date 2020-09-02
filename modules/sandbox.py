
class RegionInfo_WIP:  # TODO: This is the start of changing data types into classes. Perhaps long term project, not urgent
    def __init__(self):
        self.contig = None
        self.genbank_locus_tag = None
        self.genbank_locus_tag_list = []
        self.pseudo_locus_tag = None
        self.start = None
        self.end = None
        self.strand = None
        self.hits = []
        self.source = None
        self.note = None
        self.region_type = None

    def __str__(self):
        return False #stub

    def nucleotide_length(self):
        #TODO: Make sure this is uniform across blastp/blastx and all that
        return self.end - self.start

    def ratio_gene_length_to_avg_hit_length(self):
        hit_lengths = [hit.nucleotide_length for hit in self.hits]
        avg_hit_length = sum(hit_lengths) / len(hit_lengths)
        return self.nucleotide_length() / avg_hit_length

    def lastItem(ls):
        x = ''
        for i in ls:
            if i != "":
                x = i
        return x

    def allButTheLast(iterable, delim):
        x = ''
        length = len(iterable.split(delim))
        for i in range(0, length - 1):
            x += iterable.split(delim)[i]
            x += delim
        return x[0:len(x) - 1]

    def pseudogene_reasoning(self):
        """
        Cases:
            1. Reason: Predicted fragmentation of a single gene.
            2. Reason: ORF is %s%% of the average length of hits to this gene.
            3. Reason: Intergenic region with %s blast hits.
        """
        if self.region_type == RegionType.fragmentedpseudo:
            return "Reason: Predicted fragmentation of a single gene."
        elif self.region_type == RegionType.shortpseudo:
            return "Reason: ORF is %s%% of the average length of hits to this gene." % self.ratio_gene_length_to_avg_hit_length()
        elif self.region_type == RegionType.intergenicpseudo:
            return "Reason: Intergenic region with %s blast hits." % len(self.hits)

    def write_gff_note(self):
        note = "Note=pseudogene candidate. %s" % self.pseudogene_reasoning()
        colour = "colour=229 204 255"
        locus_tag = "locus_tag=%s" % self.pseudo_locus_tag
        genbank_locus_tags = "gbk_locus_tags=%s" % ",".join(self.genbank_locus_tag_list)
        return ";".join([note, colour, locus_tag, genbank_locus_tags])

    def gff_entry(self):  # gff-version 3 compliant entry
        seqid = "gnl|Prokka|%s" % self.contig
        source = "pseudofinder"
        type = "gene"
        start = self.start
        end = self.end
        score = "."
        strand = self.strand
        phase = "."
        attributes = self.write_gff_note()

        return "\t".join([seqid, source, type, start, end, score, strand, phase, attributes])
