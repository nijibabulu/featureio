import itertools


def complement_char(c):
    if c in complement_char.compd:
        return complement_char.compd[c]
    else:
        return c


complement_char.chars = 'acgtumrwsykvhdbnACGTUMRWSYKVHDBN'
complement_char.compl = 'tgcaakywsrmbdhvnTGCAAKYWSRMBDHVN'
complement_char.compd = dict(zip(complement_char.chars, complement_char.compl))


def reverse_complement(seq):
    return ''.join(complement_char(c) for c in seq)[::-1]


class Gene(object):
    # TODO: make a separate "from_blocks" instantiator for bed-style formats
    # TODO: gene-transcript concept
    # separate from gff-style
    def __init__(self, chrom, start, end, name, score, strand, cds_start,
                 cds_end, item_rgb, block_count, block_sizes, block_starts,
                 attrs=None, *args, **kwargs):
        self.chrom = chrom
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.score = int(score)
        self.strand = strand
        self.cds_start = int(cds_start)
        self.cds_end = int(cds_end)
        self.item_rgb = item_rgb
        self.block_count = int(block_count)
        self.block_sizes = [int(s) for s in block_sizes.split(',') if len(s)]
        self.length = sum(self.block_sizes)
        self.block_starts = [int(s) for s in block_starts.split(',') if len(s)]
        self.attrs = {} if attrs is None else attrs

        self.__dict__.update(kwargs)
        self._aux_attrs = kwargs.keys()

        self.exons = []
        self.cds_exons = []

        for start, size in zip(self.block_starts, self.block_sizes):
            start += self.start
            end = start + size
            sc = sorted([start, end, self.cds_start, self.cds_end])
            self.exons.append((start, end))
            if (sc[0] == start and sc[1] == end) or (
                    sc[0] == self.cds_start and sc[1] == self.cds_end):
                continue
            self.cds_exons.append(
                (max(self.cds_start, start), min(self.cds_end, end)))

        self.length = sum([e[1] - e[0] + 1 for e in self.exons])
        self.cds_length = sum([e[1] - e[0] + 1 for e in self.cds_exons])

        if strand == '+':
            self.fivep = self.start
            self.cds_fivep = self.cds_start
        else:
            self.fivep = self.end
            self.cds_fivep = self.cds_end

    def copy(self):
        return Gene(self.chrom, self.start, self.end, self.name, self.score,
                    self.strand, self.cds_start, self.cds_end, self.item_rgb,
                    self.block_count, ','.join(map(str, self.block_sizes)),
                    ','.join(map(str, self.block_starts)),
                    {k: self.__dict__[k] for k in self._aux_attrs})

    def modified(self, **kwargs):
        copy = self.copy()
        for a, v in kwargs.items():
            setattr(copy, a, v)
        return copy

    def _get_seq(self, seq, exons):
        gseq = ''.join(seq[c[0] - 1:c[1]]
                       for c in sorted(exons, key=lambda a: a[0]))
        if self.strand == '-':
            gseq = reverse_complement(gseq)
        return gseq

    def get_cds(self, seq):
        return self._get_seq(seq, self.cds_exons)

    def get_exons(self, seq):
        return self._get_seq(seq, self.exons)

    def locus_overlap(self, other):
        if self.start > other.end or other.start > self.end:
            return False
        if self.chrom != other.chrom:
            return False
        if self.strand != other.strand:
            return False
        return True

    def is_isoform(self, other, comparison='exons'):  # , check=False):
        if not self.locus_overlap(other):
            return False
        for a, b in itertools.product(getattr(self, comparison),
                                      getattr(other, comparison)):
            if a[0] == b[0] and a[1] == b[1]:
                return True
        return False

    def identical(self, other, comparison='exons'):
        self_exons = sorted(getattr(self, comparison))
        other_exons = sorted(getattr(other, comparison))
        return len(self_exons) == len(other_exons) and \
            all(a == b for a, b in zip(self_exons, other_exons))

    def overlap(self, other, comparison='exons'):
        if not self.locus_overlap(other):
            return False
        for a, b in itertools.product(getattr(self, comparison),
                                      getattr(other, comparison)):
            if a[0] <= b[1] and b[0] <= a[1]:
                return True
        return False

    def overlap_length(self, other, comparison='exons'):
        if not self.locus_overlap(other):
            return 0
        length = 0
        for a, b in itertools.product(getattr(self, comparison),
                                      getattr(other, comparison)):
            length += max(0, min(a[1], b[1]) - max(a[0], b[0]))
        return length

    def __str__(self):
        return '\t'.join(str(item) for item in
                         [self.chrom, self.start, self.end, self.name,
                          self.score, self.strand,
                          self.cds_start, self.cds_end, self.item_rgb,
                          self.block_count,
                          ','.join(str(s) for s in self.block_sizes),
                          ','.join(str(s) for s in self.block_starts)])
