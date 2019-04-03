import os
import pathlib

import pytest
import featureio


@pytest.mark.fasta
def test_create(fasta_dir):
    indexed_fasta = featureio.IndexedFasta(os.path.join(fasta_dir, 'random.fa'))
    assert indexed_fasta is not None


@pytest.mark.fasta
def test_extract_string(fasta_dir, output_dir):
    indexed_fasta = featureio.IndexedFasta(os.path.join(
        fasta_dir, 'GCF_000744065.1_ASM74406v1_genomic.fna'))
    seq = indexed_fasta.get_sequence('NZ_BBIY01000160.1')[20:30]
    out_file = 'GCF_000744065.1_ASM74406v1_genomic.fna_NZ_BBIY01000160.1.20-30'
    with open(os.path.join(output_dir, out_file)) as f:
        expected_substr = f.read().strip()

    assert seq == expected_substr


@pytest.mark.fasta
def test_complete_index(fasta_dir):
    indexed_fasta = featureio.IndexedFasta(os.path.join(fasta_dir, 'random.fa'))
    with open(indexed_fasta.index_filename) as f:
        lines = f.read().strip().split('\n')
    assert len(indexed_fasta) == len(lines)


def test_read_empty_file(tmp_path):
    p = tmp_path / 'test'
    p.write_text(' ')
    with pytest.raises(StopIteration):
        featureio.parse_fasta_record(p.open('r'))


def test_read_bad_index(tmp_path):
    fa = tmp_path / 'test.fa'
    fai = tmp_path / 'test.fa.fai'
    fa.write_text('>hi\nAAAAAAA')
    fai.write_text('not a well formed file')
    with pytest.raises(ValueError):
        _ = featureio.IndexedFasta(str(fa))


@pytest.mark.fasta
def test_lengths_random_fasta(fasta_dir, output_dir):
    indexed_fasta = featureio.IndexedFasta(os.path.join(fasta_dir, 'random.fa'))
    length_dict = {k: len(indexed_fasta[k]) for k in indexed_fasta.sequences()}
    with open(os.path.join(output_dir, 'random.fa-lengths')) as f:
        def process_pair(p): return tuple([p[1], int(p[0])])

        expected_length_dict = dict(process_pair(l.split()) for l in f)
    assert length_dict == expected_length_dict


@pytest.mark.fasta
def test_make_fasta_collection(fasta_dir):
    indexed_fasta_group = featureio.IndexedFastaCollection(
        [os.path.join(fasta_dir, fn) for fn in
         ['GCF_000744065.1_ASM74406v1_genomic.fna', 'random.fa']]
    )
    assert indexed_fasta_group is not None


@pytest.mark.fasta
def test_nonexistent_fasta(tmp_path):
    p: pathlib.Path = tmp_path / 'nonexistent'
    p.write_text('>a\nA')
    p.rename(tmp_path / 'existent')
    with pytest.raises(ValueError):
        missing = tmp_path / 'nonexistent'
        _ = featureio.IndexedFasta(str(missing))


@pytest.mark.fasta
def test_nonexistent_fasta_index(tmp_path):
    p: pathlib.Path = tmp_path / 'nonexistent_index'
    p.write_text('>a\nA')
    with pytest.raises(ValueError):
        _ = featureio.IndexedFasta(str(p))


@pytest.mark.fasta
def test_redundant_record_names(tmp_path):
    fa: pathlib.Path = tmp_path / 'redundant'
    fa.write_text('>redundant_name\nACGT\n>redundant_name\nACGT')
    fai: pathlib.Path = tmp_path / 'redundant.fai'
    fai.write_text('redundant_name\t4\t16\t4\t5\nredundant_name\t4\t37\t4\t5')
    with pytest.raises(ValueError):
        _ = featureio.IndexedFasta(str(fa))


@pytest.mark.fasta
def test_getitem_non_string(fasta_dir):
    indexed_fasta = featureio.IndexedFasta(os.path.join(fasta_dir, 'random.fa'))
    with pytest.raises(ValueError):
        _ = indexed_fasta[1]


@pytest.mark.fasta
def test_getitem_nonexistent_sequence(fasta_dir):
    indexed_fasta = featureio.IndexedFasta(os.path.join(fasta_dir, 'random.fa'))
    with pytest.raises(KeyError):
        _ = indexed_fasta['idontexist']


@pytest.mark.fasta
def test_collisions_collection(fasta_dir):
    try:
        _ = featureio.IndexedFastaCollection(
            [os.path.join(fasta_dir, fn) for fn in ['random.fa', 'random.fa']])
    except ValueError as ve:
        assert str(ve).startswith('Key collision')
    else:
        assert False


@pytest.mark.fasta
def test_collections_contains(fasta_dir):
    inedexed_fasta_collection = featureio.IndexedFastaCollection(
        [os.path.join(fasta_dir, fn) for fn in
         ['GCF_000744065.1_ASM74406v1_genomic.fna', 'random.fa']]
    )
    assert 'seq1' in inedexed_fasta_collection


@pytest.mark.fasta
def test_getitem_collection(fasta_dir):
    indexed_fasta = featureio.IndexedFastaCollection(
        [os.path.join(fasta_dir, fn) for fn in
         ['random.fa', 'GCF_000744065.1_ASM74406v1_genomic.fna']])
    assert indexed_fasta['seq1'].name == 'seq1'
    with pytest.raises(ValueError):
        _ = indexed_fasta[1]
    with pytest.raises(KeyError):
        _ = indexed_fasta['idontexist']


@pytest.mark.fasta
def test_keys_collection(fasta_dir):
    indexed_fasta = featureio.IndexedFastaCollection(
        [os.path.join(fasta_dir, fn) for fn in
         ['random.fa', 'GCF_000744065.1_ASM74406v1_genomic.fna']])
    assert len(indexed_fasta.keys()) == 470
    assert indexed_fasta.keys() == indexed_fasta.sequences_names()


def test_index_contains(fasta_dir):
    indexed_fasta = featureio.IndexedFasta(os.path.join(fasta_dir, 'random.fa'))
    assert 'seq1' in indexed_fasta


@pytest.mark.fasta
def test_index_keys(fasta_dir):
    indexed_fasta = featureio.IndexedFasta(os.path.join(fasta_dir, 'random.fa'))
    with open(indexed_fasta.index_filename) as f:
        lines = f.read().strip().split('\n')
    assert len(indexed_fasta.sequences()) == len(lines)


@pytest.mark.fasta
def test_seq_must_be_slice():
    seq = featureio.Seq('testname', 'ACGT')
    assert seq[1:2] == 'C'
    with pytest.raises(ValueError):
        _ = seq['no']


@pytest.mark.fasta
def test_fasta_string():
    seq = featureio.Seq('testname', 'AAAAAAAA')
    assert featureio.fasta_string(seq, wrap=2) == '>testname\nAA\nAA\nAA\nAA\n'


@pytest.mark.fasta
def test_write_fasta(tmp_path):
    seq = featureio.Seq('testname', 'AAAAAAAA')
    p: pathlib.Path = tmp_path / 'test'
    featureio.write_fasta_record(seq, p.open('w'), wrap=2)
    assert p.read_text() == '>testname\nAA\nAA\nAA\nAA\n'
