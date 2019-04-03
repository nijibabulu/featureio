#! /usr/bin/env python
import functools

from . import gene


def BedIterator(handle, cls=gene.Gene):
    for n, line in enumerate(handle):
        if line.startswith('#'):
            continue
        if not len(line.strip()):
            continue
        fields = line.strip().split('\t')
        if len(fields) != 12:
            raise ValueError(
                'Incorrect number of fields on line {}:\n{}'.format(
                    n, line))
        yield cls(*fields)


def parse_psl_line(matches, misMatches, repMatches, nCount, qNumInsert,
                   qBaseInsert, tNumInsert, tBaseInsert, strand,
                   name, qSize, qStart, qEnd, chrom, tSize, start, end,
                   block_count,
                   block_sizes, qStarts, block_starts, cls=gene.Gene):
    corr_block_starts = ','.join([str(int(s) - int(start))
                                  for s in block_starts.split(',') if len(s)])
    return cls(chrom, start, end, name, matches, strand, 0, 0, 0,
               block_count, block_sizes, corr_block_starts)


def PslIterator(handle, cls=gene.Gene):
    for line in handle:
        fields = line.strip().split()
        yield parse_psl_line(*fields, cls=cls)


def BlatPslIterator(handle, cls=gene.Gene):
    for _ in range(5):  # skip the header
        next(handle)
    for line in handle:
        yield parse_psl_line(line, cls=cls)


def AugustusGtfIterator(handle, cls=gene.Gene):
    while True:
        line = handle.readline()
        if not line:
            return
        if line.startswith("# start gene"):
            break
    while True:
        line = handle.readline()
        if not line:
            return  # stopiteration
        f = line.strip().split('\t')
        if len(f) > 1 and f[2] == 'transcript':
            chrom, start, end, strand, tid = (
                f[0], int(f[3]), int(f[4]), f[6], f[-1])
            bst, bsz = [[], []]
            cs, ce = [None, None]
            seq = ''
            gene_id = None
            while True:
                line = handle.readline()
                if not line:
                    return  # stopiteration
                # if we have already seen protein sequence or we are initiating
                # the protein sequence add to the current sequence
                if line.startswith('# protein sequence') or len(seq):
                    istart = (line.index('[') + 1) if line.count('[') else 2
                    iend = line.index(']') if line.count(']') else len(line)
                    seq += line[istart:iend].strip()
                else:
                    # TODO: change to a non-block initializer to avoid these
                    # calcluations
                    f = line.strip().split('\t')
                    if f[2] == 'exon':
                        gene_id = f[-1].split()[-1].strip('";')
                        bst.append(int(f[3]) - start)
                        bsz.append(int(f[4]) - int(f[3]))
                    elif f[2] == 'CDS':
                        gene_id = f[-1].split()[-1].strip('";')
                        cs = min(cs, int(f[3])) if cs else int(f[3])
                        ce = max(cs, int(f[4])) if ce else int(f[4])
                if seq is not None and ']' in line:
                    yield cls(
                        chrom, start, end, tid, 0, strand, cs, ce, 0, len(bst),
                        ','.join(str(x) for x in bsz),
                        ','.join(str(x) for x in bst),
                        seq=seq, gene_id=gene_id)
                    break


'''
def GFF3Iterator(handle, cls=Gene):
    genes = {}
    children = {}
    for lineno,line in enumerate(handle):
        if line.startswith('#'):
            continue
        try:
            chrom,source,ftype,start,end,score,strand,_,a = line.strip().split('\t')
            attrs = dict(tuple(a.strip() for a in attr.split('='))
                         for attr in a.split(';'))
        except ValueError:
            raise ValueError('Could not parse line {}: {}'.format(lineno,line))

        if ftype == 'gene':
            gene_name = attrs.get('ID', attrs.get('Name', attrs.get('gene_id')))
            if gene_name is None:
                raise ValueError('No ID, Name or gene_id found for gene on '
                                 'line {}:\n{}'.format(lineno, line))
            gene = children.get(gene_name, genes.setdefault(
                gene_name, Gene(chrom, start, end, gene_name, score, strand,
                                -1, -1, 0, 0, 0, [], attrs)))

        gene = genes.setdefault(gene_name, Gene(chrom, start, end, ))
'''

_readers = {"bed12": BedIterator, "psl": PslIterator,
            "blatpsl": BlatPslIterator,
            "augustusgtf": AugustusGtfIterator}


# "gff3": GFF3Iterator}


def parse(maybe_handle, format, mode='r', cls=gene.Gene, **kwargs):
    # type: (Union[TextIO, str], str, str, Callable[[...], gene.Gene], ...) -> List[gene.Gene]
    # this can be better handled with contextlib.contextmanager
    if isinstance(maybe_handle, str):
        fp = open(maybe_handle, mode, **kwargs)
    else:
        fp = maybe_handle

    if format in _readers:
        gen = _readers[format]
        i = gen(fp, cls=cls)

        for g in i:
            yield g
    else:
        raise ValueError('Unknown format {}. Should be one of {}'.format(
            format, ','.join(_readers.keys())))


def to_dict(gene_iterable):
    return {g.name: g for g in gene_iterable}


class GeneWriter(object):
    def __init__(self, handle, **kwargs):
        self.handle = handle

    def write_header(self):
        pass

    def write_genes(self, genes):
        for i, gene in enumerate(genes):
            self.count = i + 1
            self.write_gene(gene)

    def write_gene(self, _):
        raise NotImplementedError('{} does not implement write_gene'.format(
            self.__class__))

    def write_footer(self):
        pass

    def write_file(self, genes):
        self.write_header()
        for gene in genes:
            self.write_gene(gene)
        self.write_footer()


class AugustusExonHintWriter(GeneWriter):
    def __init__(self, handle, cds_exons=True, feature_type='exon',
                 source='featureio', augustus_source='E', priority=4):
        super(AugustusExonHintWriter, self).__init__(handle)
        self.exon_attr = 'cds_exons' if cds_exons else 'exons'
        self.feature_type = feature_type
        self.source = source
        self.priority = priority
        self.augustus_source = augustus_source

    def write_gene(self, gene):
        for exon in getattr(gene, self.exon_attr):
            attrs = 'grp={};pri={};src={}'.format(
                gene.name, self.priority, self.augustus_source)
            self.handle.write('\t'.join(str(x) for x in [
                gene.chrom,
                self.source,
                self.feature_type,
                exon[0],
                exon[1],
                '.',
                gene.strand,
                '.',
                attrs]))
            self.handle.write('\n')


class Bed12Writer(GeneWriter):
    def write_gene(self, gene):
        self.handle.write('\t'.join(str(item) for item in
                                    [gene.chrom, gene.start, gene.end,
                                     gene.name, gene.score, gene.strand,
                                     gene.cds_start, gene.cds_end,
                                     gene.item_rgb, gene.block_count,
                                     ','.join(str(s) for s in gene.block_sizes),
                                     ','.join(str(s) for s in
                                              gene.block_starts)]) + '\n')


class GFF3Writer(GeneWriter):
    def __init__(self, handle, source='GFF3Conv', **kwargs):
        super(GFF3Writer, self).__init__(handle, **kwargs)
        self.source = source

    def write_feature(self, gene, ftype, start, end, feature_number=1,
                      toplevel=False, phase='.', score='.',
                      name=None, attrs=None):
        attrs = dict() if attrs is None else attrs
        if toplevel:
            attrs.update(gene.attrs)
        self.handle.write('\t'.join(map(str, [
            gene.chrom, self.source, ftype, start, end, score, gene.strand,
            phase])))
        name_base = name or gene.name
        if not toplevel:
            attrs.setdefault('Parent', gene.name)
            attrs['ID'] = '{}.{}.{}'.format(gene.name, ftype, feature_number)
            if name is not None:
                attrs['Name'] = name
        else:
            attrs['Name'] = name_base
            attrs['ID'] = gene.name
        self.handle.write('\t' + ';'.join('{}={}'.format(k, v)
                                          for k, v in attrs.items()))
        self.handle.write('\n')

        return attrs['ID']

    def write_header(self):
        self.handle.write('##gff-version 3\n')

    @staticmethod
    def _sorted_cds(gene):
        return sorted(gene.cds_exons, key=lambda c: c[0],
                      reverse=gene.strand == '-')

    def cds_writer(self, gene, transcript_id):
        def write_cds(phase, enumerated_cds):
            feature_number, cds = enumerated_cds
            self.write_feature(gene, 'CDS', cds[0], cds[1], phase=phase,
                               feature_number=feature_number,
                               attrs={'Parent': transcript_id})
            return (3 - ((cds[1] - cds[0] - phase) % 3)) % 3

        return write_cds

    def write_gene(self, gene):
        self.write_feature(gene, 'gene', gene.start, gene.end, toplevel=True)
        transcript_id = self.write_feature(gene, 'mRNA', gene.start, gene.end)
        for n, e in enumerate(gene.exons):
            self.write_feature(gene, 'exon', e[0], e[1], n,
                               attrs={'Parent': transcript_id})
        functools.reduce(self.cds_writer(gene, transcript_id),
                         enumerate(GFF3Writer._sorted_cds(gene)), 0)


_writers = {"bed12": Bed12Writer,
            "augustus_exon_hints": AugustusExonHintWriter,
            "gff3": GFF3Writer}
valid_writers = _writers.keys()
valid_readers = _readers.keys()


def write(genes, maybe_handle, format, mode='w', **kwargs):
    if isinstance(maybe_handle, str):
        fp = open(maybe_handle, mode, **kwargs)
    else:
        fp = maybe_handle

    if format in _writers:
        writer = _writers[format](fp, **kwargs)
        writer.write_file(genes)
    else:
        raise ValueError('Unknown format {}. Should be one of {}'.format(
            format, ','.join(_writers.keys())))


def main():
    from Bio import SeqIO

    s = str(SeqIO.read('unmask_split_nv2i5/Chr1.fa', 'fasta').seq)

    for g in parse('bookends_rb/Chr1.BE0.1.gff', 'augustusgtf'):
        print(g.name)
        print(g.strand)
        cds = g.get_cds(s)
        print(len(cds))
        print(g.cds_length)
        print(cds)


if __name__ == '__main__':
    main()
