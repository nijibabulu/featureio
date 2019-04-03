import os
from typing import TextIO, List, Dict

import attr


class Seq(object):
    """Placeholder biological sequence object."""

    def __init__(self, name, sequence, description=None):
        """Initialize a seuqence object

        :param name: Name of the sequence
        :param sequence: The actual sequence, without spaces or newlines
        :param description: an optional description. This equates to the
            information after the first space in a fasta file
        """
        self.name = name
        self.sequence = sequence
        self.description = description

    def __len__(self):
        return len(self.sequence)

    def __getitem__(self, key):
        if isinstance(key, int) or isinstance(key, slice):
            return self.sequence[key]
        else:
            raise ValueError("Cannot getitem from a Seq with {}".format(
                type(key)
            ))


def parse_fasta_record(file: TextIO) -> Seq:
    """Read a single fasta record from an opened file

    :param file: an opened file object that implements read() and seek().
        streams will not work on this function
    :return: a ``Seq`` object
    :raises StopIteration: when no header could be found before the end of file.
    """
    while True:
        c = file.read(1)
        if len(c) == 0:
            raise StopIteration
        if c == '>':
            break
    header = file.readline()
    fields = header.strip().split(maxsplit=1)
    id = fields[0]
    description = fields[1] if len(fields) > 1 else None
    sequence = ''
    end_of_file = False
    try:
        while True:
            c = file.read(1)
            if c.isspace():
                continue
            if len(c) == 0:
                raise StopIteration
            sequence += c
            if c == '>':
                break
    except (EOFError, StopIteration):
        end_of_file = True
    if not end_of_file:
        file.seek(file.tell() - 1)
        sequence = sequence[:-1]
    return Seq(id, sequence, description)


def fasta_string(seq: Seq, wrap: int = 100) -> str:
    """Convert a seq object to a fasta-formatted string.

    :param seq: a ``Seq`` object.
    :param wrap: number of lines to wrap the sequence when outputing
    :return:
    """
    description = " " + seq.description if seq.description is not None else ""
    wrapped_seq = '\n'.join(seq.sequence[i:i + wrap]
                            for i in range(0, len(seq.sequence), wrap))
    return f">{seq.name}{description}\n{wrapped_seq}\n"


def write_fasta_record(seq: Seq, file: TextIO, wrap=100) -> None:
    """Write a fasta record to a file

    :param seq: a ``Seq`` object.
    :param file: a file opened for writing.
    :param wrap: number of lines to wrap the sequence when outputing
    :return:
    """
    file.write(fasta_string(seq, wrap))


@attr.s
class FastaIndexRecord(object):
    """A fasta index record.

    This represents a single fasta index record. This should be used in the
    context of a ``IndexedFasta`` object and is internal.

    :param name: the name of the sequence
    :param length: length of the sequence
    :param offset: the position in the corresponding file
    :param line_bases: the number of bases in a given line
    :param line_width: the total number of bytes per line (usually
        1+``line_bases``)
    """
    name: str = attr.ib()
    length: int = attr.ib(converter=int)
    offset: int = attr.ib(converter=int)
    line_bases: int = attr.ib(converter=int)
    line_width: int = attr.ib(converter=int)

    @classmethod
    def from_string(cls, string):
        """Make a FastaIndexRecord from a string. Used for parsing.

        :param string: a line representing a fasta index record
        :return: a FastaIndexRecord object
        :rtype: FastaIndexRecord
        """
        fields = string.split('\t')
        if len(fields) != len(attr.fields(cls)):
            raise ValueError(f"Malformed fasta index record. Expected "
                             f"{len(attr.fields(cls))} but got "
                             f"{len(fields)}. Offending line was: \n{string}")
        return FastaIndexRecord(*fields)


class IndexedFasta(object):
    """An index to a fasta file and its file.

    This object can be used to retrieve sequences from large files without
    loading the fasta into memory.
    """

    def __init__(self, filename: str):
        """Initialize an IndexedFasta object.

        :param filename: a path to a fasta file which has an associated ``.fai``
            file in the same path.
        """
        self.filename = filename
        self.index_filename = filename + ".fai"
        if not os.path.exists(self.filename):
            raise ValueError(f"{self.filename} does not exist")
        if not os.path.exists(self.index_filename):
            raise ValueError(f"No {self.index_filename} found! Indexed fasta "
                             f"files require an index. See e.g. samtools faidx "
                             f"for help.")
        with open(self.index_filename) as f:
            records = [FastaIndexRecord.from_string(line) for line in f]
        self.records = {record.name: record for record in records}
        if len(self.records) != len(records):
            raise ValueError(f"Non-unique sequence names in {self.filename}")

    def __len__(self):
        return len(self.records)

    def __getitem__(self, item):
        """Get a sequence by name

        :param str item: name of the sequence
        :return: a featureio.Seq object
        """
        if isinstance(item, str):
            return self.get_sequence(item)
        else:
            raise ValueError("Can only getitem of type str")

    def __contains__(self, item: str) -> bool:
        """indicate whether a sequence is contained in the file

        :param item: A name of a sequence
        :return: True if present
        """
        return self.records.__contains__(item)

    def sequences(self):
        """Return all sequence names contained in the index"""
        return self.records.keys()

    def get_sequence(self, name: str) -> Seq:
        """Retrieve a sequence from an indexed fasta

        :param name: The name of the sequence in the file
        :return: A Seq object
        """
        record = self.records.get(name, None)
        if record is None:
            raise KeyError(f"No such sequence {name} in {self.filename}")
        with open(self.filename) as f:
            # FIXME: The fasta index from samtools points to the sequence
            #  itself, not the header. This is a problem if we want to parse
            #  the full header. Here we seek back in the file until we find
            #  a header carat, but this could go wrong for some files. Another
            #  possibility is to take until the next record, but that violates
            #  our code reuse

            offset = record.offset
            f.seek(offset)
            while f.read(1) != '>':
                offset -= 1
                f.seek(offset)
            f.seek(offset)
            return parse_fasta_record(f)


class IndexedFastaCollection(object):
    """A collection of indexed fasta sequences"""

    def __init__(self, files: List[str]):
        """Initialized an IndexedFastaCollection.

        :param files: a list of fasta files with associated ``.fai`` files.
        """
        indexed_fastas = [IndexedFasta(file) for file in files]
        self.index_map: Dict[str, IndexedFasta] = {}
        for indexed_fasta in indexed_fastas:
            for k in indexed_fasta.sequences():
                if k in self.index_map:
                    raise ValueError(f"Key collision: sequence name {k} appears more than once.")
                self.index_map[k] = indexed_fasta

    def __contains__(self, item) -> bool:
        """Indicate whether a sequence is in any file in the collection

        :param item: a sequence name
        :return: True if present in any file, False if not
        """
        return self.index_map.__contains__(item)

    def __getitem__(self, item):
        """Get a sequence by name

        :param str item: name of the sequence
        :return: a featureio.Seq object
        """
        if isinstance(item, str):
            return self.get_sequence(item)
        else:
            raise ValueError("Can only getitem of type str")

    def keys(self):
        """get the names of all sequences in the index"""
        return self.index_map.keys()

    def sequences_names(self):
        """get the names of all sequences in the index"""
        return self.keys()

    def get_sequence(self, name):
        """Retrieve a sequence object from the collection

        :param str name: The name of a sequence
        :return: a featureio.Seq object
        """
        index = self.index_map.get(name, None)
        if index is None:
            raise KeyError(f"No such sequence {name} found in any index!")
        return index.get_sequence(name)
