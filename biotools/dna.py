from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from biotools.oligo import Oligo

from biotools.bio_utils import rev_comp, validate_sequence, validate_annotations
from biotools.bio_enums import BioMolecule, BioProperty, BioOrientation
from biotools.bio_annotation import BioAnnotation, Block, BioPool, reverse_annotations
from biotools.bio_exceptions import (
    InvalidSequence, InvalidInstantiationException, InvalidAnnotationException
)
import bisect
from itertools import product
import re
from hashlib import sha256
import json

class DNA():
    """
    Class representing the BioMolecule DNA. DNA is created by passing in a string comprised 
    of one or more of the characters defined in NTs.

    Parameters
    ----------
    seq: str
        A string sequence of DNA. Must contain only valid nucleotide or degenerate base
        characters as specified in valid_chars
    type: BioMolecule
        An identifier Enum of the type DNA
    circular: BioProperty
        A boolean of whether or not the sequence is circular (a plasmid)
    strandedness: BioProperty
        A property stating whether or not the DNA is single- or double-stranded
    annotations: list[BioAnnotation | Block]
        A list of annotations present on the DNA
    offsets: tuple[tuple[int,int], tuple[int,int]]
        Offsets defining overhangs at either end of the DNA fragment for cut DNA
    parent: DNA
        The molecule of DNA from which self was generated, if applicable

    Examples
    --------
    >>> dna = DNA(seq="ATGCATGCATGC", name="foo")
    >>> dna.sequence()
    ATGCATGCATGC
    TACGTACGTACG

    >>> dna_slice = dna[:10]
    >>> dna_slice.sequence()
    ATGCATGCAT
    TACGTACGTA

    >>> foo
    """

    def __init__(
        self, 
        seq: str, 
        name: str = "", 
        type: BioMolecule = BioMolecule.DNA,
        circular: BioProperty = BioProperty.LINEAR,
        strandedness: BioProperty = BioProperty.DOUBLE_STRANDED,
        annotations: list[BioAnnotation | Block] | None = None,
        offsets: tuple[tuple[int,int], tuple[int,int]] = ((0,0),(0,0)),
        parent: DNA = None
    ):
        """Default constructor."""
        if validate_sequence(seq, type):
            self.seq = seq.upper() # NOTE: Sequences should be capitalized
        else:
            raise InvalidSequence("Must provide a valid DNA sequence!")
        self.name = name
        if type not in BioMolecule:
            raise InvalidInstantiationException("Must use a valid BioMolecule!")
        self.type = type
        if circular not in BioProperty:
            raise InvalidInstantiationException("Must use a valid BioProperty!")
        self.circular = circular
        if strandedness not in BioProperty:
            raise InvalidInstantiationException("Must use a valid BioProperty!")
        self.strandedness = strandedness
        if annotations is not None:
            if not isinstance(annotations, list):
                raise InvalidInstantiationException("Must provide a list of BioAnnotations!")
            for annot in annotations:
                if not isinstance(annot, BioAnnotation) and not isinstance(annot, Block):
                    raise InvalidInstantiationException("Must use a valid BioAnnotation!")
        if annotations is not None:
            if validate_annotations(self, annotations):
                self.annotations = sorted(annotations)
        else:
            self.annotations = []
        """
        Offsets is a tuple of ((top_5', top_3'), (bottom_3', bottom_5')), where any value other than
        ((0, 0), (0, 0)) implies that this was created from a reaction and/or contains overhangs.
        """
        self.offsets = offsets
        self.parent = parent

    @property
    def length(self) -> int:
        """The length of the DNA sequence."""
        return len(self.seq)
    
    @property
    def top_strand(self) -> str:
        """The top strand shown 5' -> 3'."""
        if self.offsets:
            return self.seq[self.offsets[0][0] : (self.length - self.offsets[0][1])]
        return self.seq
    
    @property
    def bottom_strand(self) -> str:
        """The bottom strand shown 5' -> 3'."""
        if self.strandedness == BioProperty.SINGLE_STRANDED:
            raise Exception("Cannot return bottom strand of a single-stranded sequence!")
        else:
            if self.offsets:
                return rev_comp(self.seq)[self.offsets[1][1] : (self.length - self.offsets[1][0])]
            return rev_comp(self.seq)

    def __repr__(self) -> str:
        """Representation of DNA, generically as name, BioMolecule, length, annotation number, 
        and strandedness."""
        if self.is_circular():
            return f"{{ {self.name} | {self.type.value} | [{self.length} bp] | [{len(self.annotations)} Annotation(s)] | O }}"
        elif self.is_double_stranded():
            return f"{{ {self.name} | {self.type.value} | [{self.length} bp] | [{len(self.annotations)} Annotation(s)] | ==> }}"
        elif self.is_cut():
            return f"{{ {self.name} | {self.type.value} | [{self.length} bp] | [{len(self.annotations)} Annotation(s)] | Cut }}"
        else:
            return f"{{ {self.name} | {self.type.value} | [{self.length} bp] | [{len(self.annotations)} Annotation(s)] | --> }}"
    
    def __getitem__(self, index) -> slice:
        """Grabs and returns a slice of DNA, preserving annotations and re-indexing them
        as appropriate.
        NOTE: Only reversing (-1) is supported as a step when slicing. All other values
        raise a ValueError."""
        annotations: list[BioAnnotation] = []
        if isinstance(index, int): # Single integer provided
            if index > self.length - 1:
                raise IndexError("Index provided is out of bounds of sequence length!")
            if index < 0: # Valid to have index like [-int] if abs(int) < sequence length
                if abs(index) > self.length:
                    raise IndexError("Index provided is out of bounds of sequence length!")
                index = len(self.seq) + index
            for annot in self.annotations:
                new_span = None
                if annot.span[0] < annot.span[1]:
                    if index > annot.span[0] and index < annot.span[1]:
                        new_span = (index, index + 1)
                else:
                    if index > annot.span[1] and index < annot.span[0]:
                        new_span = (index, index + 1)
                if new_span is not None:
                    if isinstance(annot, BioAnnotation):
                        annotations.append(BioAnnotation(new_span, annot.name, annot.orientation))
                    elif isinstance(annot, Block):
                        annotations.append(Block(new_span, annot.name, annot.orientation, annot.pool))
            return DNA(
                self.seq[index], 
                self.name, 
                self.type, 
                self.circular, 
                self.strandedness, 
                annotations, 
                self.offsets, 
                self
            )
        elif isinstance(index, slice):
            start = index.start if index.start is not None else 0
            if start < 0:
                if abs(start) > self.length:
                    raise IndexError("Index provided is out of bounds of sequence length!")
                start = len(self.seq) + start
            stop = index.stop if index.stop is not None else len(self.seq)
            new_length = stop - start

            if index.step is not None: # Check that step is -1
                if index.step != -1:
                    raise ValueError("Only reversing is allowed as a step for slicing (e.g. dna[::-1])!")
                
            # Check orientation, and continue as appropriate
            if not self.is_circular():
                if start > stop:
                    raise IndexError("For linear DNA, first index must be less than the second!")
                elif start < 0 or stop > len(self.seq):
                    raise IndexError("Cannot have indices beyond bounds of sequence length!")
                # Grab annotations that fall within slice
                for annot in self.annotations:
                    new_span = None
                    if annot.span[0] >= start and annot.span[0] < stop: # Start of span is within start and stop
                        if annot.span[1] <= stop:
                            new_span = (annot.span[0]-start, annot.span[1]-start)
                        else:
                            new_span = (annot.span[0]-start, new_length)
                    elif annot.span[0] < start:
                        if annot.span[1] > start and annot.span[1] <= stop:
                            if annot.span[1] <= stop:
                                new_span = (0, annot.span[1]-start)
                            else:
                                new_span = (0, new_length)
                        elif annot.span[1] > stop:
                            new_span = (0, new_length)
                    elif annot.span[1] > start and annot.span[1] < stop: # Start occurs earlier than new start
                        if annot.span[1] <= stop:
                            new_span = (0, annot.span[1]-start)
                        else:
                            new_span = (0, new_length)
                    if new_span is not None:
                        if new_span[0] != new_span[1]:
                            if isinstance(annot, BioAnnotation):
                                annotations.append(BioAnnotation(new_span, annot.name, annot.orientation))
                            elif isinstance(annot, Block):
                                annotations.append(Block(new_span, annot.name, annot.orientation, annot.pool))
                seq = self.seq[start:stop]
                if index.step is not None:
                    annotations = reverse_annotations(annotations, new_length)
                    seq = seq[::-1]
                return DNA(
                    seq, 
                    self.name, 
                    self.type, 
                    self.circular, 
                    self.strandedness, 
                    annotations, 
                    self.offsets, 
                    self
                )
            else: # Circular
                if start < 0 or stop > self.length or stop < 0 or start > self.length:
                    raise IndexError("Cannot have indices beyond bounds of sequence length!")
                if start > stop:
                    new_length = (self.length - start) + stop
                    for annot in self.annotations:
                        new_span: tuple[int, int] = None
                        if annot.span[0] >= start: # Span starts within slice bounds of circular sequence
                            if annot.span[1] <= stop:
                                new_span = (annot.span[0] - start, (self.length - start) + annot.span[1])
                            elif annot.span[1] > stop and annot.span[1] > start:
                                new_span = (annot.span[0] - start, annot.span[1] - start)
                            else: # Stop must be < length of seq
                                new_span = (annot.span[0] - start, new_length)
                        else: # Span's start is less than starting index) NOTE: Here we have to check if the spans are the same because of re-indexing weirdness
                            if annot.span[0] <= stop: # Case where we truncate the annotation
                                if annot.span[1] > start: # We don't cut off the annotation
                                    if annot.span[0] < start:
                                        if annot.span[0] < annot.span[1]:
                                            if new_length - (stop - annot.span[0]) == new_length:
                                                new_span = (0, annot.span[1] - start)
                                            else:
                                                new_span = (new_length - (stop - annot.span[0]), annot.span[1] - start)
                                        else:
                                            new_span = (annot.span[1] - start , new_length - (stop - annot.span[0]))
                                    else:
                                        new_span = (annot.span[1] - start, stop - annot.span[0])
                                else: # Cut off the annotation at the start
                                    if annot.span[1] >= stop:
                                        if annot.span[0] < annot.span[1]:
                                            new_span = (self.length - start + annot.span[0], self.length - start + annot.span[0] + (stop - annot.span[0]))
                                        else:
                                            new_span = (0, stop - annot.span[0])
                                    else:
                                        new_span = (self.length - start + annot.span[0], self.length - start + annot.span[0] + (annot.span[1] - annot.span[0]))
                            else: # We cut off the beginning of the annotation
                                if annot.span[1] > start: # We don't cut off the annotation
                                    if annot.span[1] > annot.span[0]:
                                        new_span = (0, annot.span[1] - start)
                                    else:    
                                        new_span = (annot.span[1] - start, new_length)
                                elif annot.span[0] > annot.span[1]: # Cut off the annotation at the start
                                    if annot.span[0] < start:
                                        if annot.span[1] > stop:
                                            new_span = (0, self.length - start + stop)
                                        else:
                                            new_span = (0, self.length - start + annot.span[1])
                                    else:
                                        new_span = (self.length - start, new_length)
                        if new_span is not None:
                            if new_span[0] != new_span[1]:
                                if isinstance(annot, BioAnnotation):
                                    annotations.append(BioAnnotation(new_span, annot.name, annot.orientation))
                                elif isinstance(annot, Block):
                                    annotations.append(Block(new_span, annot.name, annot.orientation, annot.pool))
                    seq = self.seq[start:]+self.seq[:stop]
                    if index.step is not None:
                        annotations = reverse_annotations(annotations, new_length)
                        seq = seq[::-1]
                    return DNA(
                        seq,
                        self.name,
                        self.type,
                        BioProperty.LINEAR,
                        self.strandedness,
                        annotations,
                        self
                    )
                else: # Regular workflow because the slice is like a linear fragment (start < stop)
                    for annot in self.annotations:
                        new_span: tuple[int, int] = None
                        if annot.span[0] < annot.span[1]:
                            if annot.span[0] >= start and annot.span[0] < stop: # Start of span is within start and stop
                                if annot.span[1] <= stop:
                                    new_span = (annot.span[0]-start, annot.span[1]-start)
                                else:
                                    new_span = (annot.span[0]-start, new_length)
                            elif annot.span[0] < start:
                                if annot.span[1] > start and annot.span[1] <= stop:
                                    if annot.span[1] <= stop:
                                        new_span = (0, annot.span[1]-start)
                                    else:
                                        new_span = (0, new_length)
                                elif annot.span[1] > stop:
                                    new_span = (0, new_length)
                            elif annot.span[1] > start and annot.span[1] < stop: # Start occurs earlier than new start
                                if annot.span[1] <= stop:
                                    new_span = (0, annot.span[1]-start)
                                else:
                                    new_span = (0, new_length)
                        else: # End span > beginning span
                            if annot.span[0] > start and annot.span[0] >= stop:
                                if annot.span[1] >= start and annot.span[1] <= stop:
                                    new_span = (0, annot.span[1] - start)
                                elif annot.span[1] >= start and annot.span[1] > stop: # Greater than stop
                                    new_span = (0, new_length)
                            elif annot.span[0] <= stop and annot.span[0] > start and annot.span[1] >= start:
                                if (annot.span[1] - start) == 0:
                                    new_span = (annot.span[0] - start, new_length)
                                else:
                                    new_span = (annot.span[0] - start, annot.span[1] - start)
                            elif annot.span[0] <= stop and annot.span[0] > start and annot.span[1] < start:
                                new_span = (annot.span[0] - start, stop - start)
                        if new_span is not None:
                            if new_span[0] != new_span[1]:
                                if isinstance(annot, BioAnnotation):
                                    annotations.append(BioAnnotation(new_span, annot.name, annot.orientation))
                                elif isinstance(annot, Block):
                                    annotations.append(Block(new_span, annot.name, annot.orientation, annot.pool))
                seq = self.seq[start:stop]
                if index.step is not None:
                    annotations = reverse_annotations(annotations, new_length)
                    seq = seq[::-1]
                return DNA(
                    seq,
                    self.name,
                    self.type,
                    BioProperty.LINEAR,
                    self.strandedness,
                    annotations,
                    self.offsets,
                    self
                )

    def __hash__(self) -> int:
        """Hashes the sequence, defined as the addition of all hashes of
        the name, sequence, BioProperties, and annotations."""
        hash_num: int = 0
        for key in self.__dict__.keys():
            if key == "annotations":
                for annot in self.__dict__[key]:
                    hash_num += hash(annot)
            elif key in ["type", "circular", "strandedness"]:
                hash_num += int(sha256(self.__dict__[key].name.encode()).hexdigest(), 16)
            elif key == "offsets":
                if self.__dict__[key]:
                    for tup in self.__dict__[key]:
                        hash_num += int(sha256(json.dumps(tup).encode("utf-8")).hexdigest(), 16)
            else:
                hash_num += hash(self.__dict__[key])
        return hash_num

    def __eq__(self, oth: DNA) -> bool:
        """Two DNA are equal if all their properties are equal."""
        if oth is self:
            return True
        if not isinstance(oth, DNA):
            return False
        if (self.seq, self.name, self.type, self.circular, self.strandedness, self.annotations, self.parent) == \
           (oth.seq, oth.name, oth.type, oth.circular, oth.strandedness, oth.annotations, oth.parent):
            return True
        if not ((self.seq, self.name, self.type, self.circular, self.strandedness, self.annotations, self.parent) == \
               (oth.seq, oth.name, oth.type, oth.circular, oth.strandedness, oth.annotations, oth.parent)) and \
               len(self) == len(oth):  # If they're the same length, could possibly be reindexes of each other
            if self.circular != oth.circular:
                return False
            elif not self.is_circular():  # Linear, so can't re-index
                return False
            return self._check_reindex_equality(oth)  # Same length with different annotations is same sequence
        return False

    def _check_reindex_equality(self, oth: DNA) -> bool:
        """Checks equality by reindexing oth length - 1 times to confirm the sequences
        are never identical.
        """
        for i in range(1, len(oth) - 1):
            if self.seq == oth.reindex(i, False).seq:
                return True
        return False

    def __add__(self, oth: DNA) -> DNA:
        """Generates a new DNA molecule just taking the sequence. 
        This works on linear and circular DNA.
        """
        return DNA(
            self.seq + oth.seq,
            self.name,
            self.type,
            self.circular,
            self.strandedness,
            self.annotations + oth.annotations,
            self.offsets,
            self
        )

    def __len__(self) -> int:
        return self.length

    def __lt__(self, oth: DNA) -> bool:
        """DNA are compared via their hash."""
        if not isinstance(oth, DNA):
            raise Exception(f"Cannot compare {type(self).__name__} with {type(oth).__name__}!")
        return hash(self) < hash(oth)

    def sequence(self) -> None:
        """Prints a string of the DNA sequence."""
        if self.strandedness == BioProperty.SINGLE_STRANDED:
            print(self.seq)
        elif self.strandedness == BioProperty.CUT:
            # NOTE: Products can probably only be double-stranded? Should I check this?
            print(self.offsets[0][0]*" "+self.seq[self.offsets[0][0]:self.length - self.offsets[0][1]])
            bottom = rev_comp(self.seq)[::-1]
            print(self.offsets[1][0]*" "+bottom[self.offsets[1][0]:len(bottom) - self.offsets[1][1]])
        else:
            print(self.seq)
            print(rev_comp(self.seq)[::-1])

    def info(self) -> None:
        """Prints out a summary of the DNA sequence."""
        # TODO: Update to improve, and return sequence instead of print?
        print("Name: ", self.name)
        print("Sequence: ", self.seq)
        print("Reverse Complement: ", rev_comp(self.seq))
        print("Length: ", self.length)
        print("Type: ", self.type.value)
        print(f"Circular: {self.is_circular()}")
        print(f"Strandedness: {self.strandedness.value}")

    def copy(self) -> DNA:
        """Returns a copy of the DNA sequence."""
        return DNA(self.seq, self.name, self.type, self.circular, \
                   self.strandedness, self.annotations, self.offsets, self.parent)
    
    def rev_comp(self) -> DNA:
        """Returns a DNA copy of the reverse complement of the sequence."""
        if not self.is_double_stranded():
            raise Exception("Cannot generate reverse complement of single-stranded molecule!")
        else:
            # Generate inverted annotations 
            annots: list[BioAnnotation] = reverse_annotations(self.annotations, self.length, True)
            return DNA(rev_comp(self.seq), self.name+"_revcomp", self.type, self.circular, self.strandedness, annots, self.offsets)
        
    def is_circular(self) -> bool:
        return self.circular == BioProperty.CIRCULAR
    
    def is_double_stranded(self) -> bool:
        return self.strandedness == BioProperty.DOUBLE_STRANDED
    
    def is_cut(self) -> bool:
        return self.strandedness == BioProperty.CUT
    
    def has_pool(self) -> bool:
        for annot in self.annotations:
            if isinstance(annot, Block):
                return True
        return False
    
    def get_pools(self) -> list[DNA] | None:
        """Returns all sequences generated from Blocks. If more than one Block exists,
        returns all possible combinations of sequences generated.
        """

        if not self.has_pool():
            return None

        next_idx: int = 0
        self_seqs: list[str] = []
        groups: list[list[DNA]] = []

        # Iterate over annotations, finding Blocks and adding sliced self.seq in between Blocks
        for i, annot in enumerate(self.annotations):
            if isinstance(annot, Block):
                if i == 0:
                    next_idx = annot.span[1]
                    seq = self.seq[: annot.span[0]]
                    self_seqs.append(seq)
                else:
                    seq = self.seq[next_idx : annot.span[0]]
                    self_seqs.append(seq)
                    next_idx = annot.span[1]
                groups.append(annot.pool.seqs)

        self_seqs.append(self.seq[next_idx :]) # Final sequence with no block seqs after

        # Generate all combinations of Blocks and resulting sequences
        seqs: list[DNA] = []
        block_combos = product(*groups)
        for combo in block_combos:
            new_seq: str = ""
            for i in range(len(combo)):
                new_seq += self_seqs[i]
                new_seq += combo[i].seq
            new_seq += self_seqs[len(combo)]
            seqs.append(
                DNA(
                    seq=new_seq,
                    circular=self.circular,
                    strandedness=self.strandedness,
                    offsets=self.offsets,
                    parent=self
                )
            )

        seqs.append(DNA(seq=self.seq)) # Add our own seq because it's not captured in the Block combinations
        return seqs

    def has_annotation(self, annotation: str) -> bool:
        """Checks if self contains anntotation."""
        annotation_names: list[str] = [a.name for a in self.annotations]
        if annotation in annotation_names:
            return True
        return False

    def add_annotations(self, *annotations: BioAnnotation | Block) -> None:
        """Adds an annotation or annotations to the DNA sequence.

        Parameters
        ----------
        annotation: BioAnnotation | Block
            One or more annotations to add to self

        Returns
        -------
        None
        """

        # Spans can't be same
        for annot in annotations:
            if annot.span[0] == annot.span[1]:
                raise InvalidAnnotationException("Must not have same span indices!")

        # Validate tuple based on circularity
        if not self.is_circular():
            for annot in annotations:
                if annot.span[0] > annot.span[1]:
                    raise InvalidAnnotationException("First index must be less than the second index on a linear sequence!")

        # Validate indices
        for annot in annotations:
            if not all(sp >= 0 and sp <= self.length for sp in annot.span):
                raise InvalidAnnotationException("Span integers must be valid!")

        # Annotations can't have the same name
        for annot in annotations:
            if self.has_annotation(annot.name):
                raise ValueError("Annotations must have unique names!")

        # Check that no Blocks will overlap, if applicable
        for annot in annotations:
            if isinstance(annot, Block):
                for a in self.annotations:
                    if isinstance(a, Block):
                        if (annot.span[0] > a.span[0] and annot.span[0] < a.span[1]) or \
                           (annot.span[1] > a.span[0] and annot.span[1] < a.span[1]):
                                raise InvalidAnnotationException(f"Block {annot.name} and Block {a.name} cannot overlap!")

        for annot in annotations:
            valid = True
            for a in self.annotations:
                if isinstance(a, type(annot)):
                    if annot == a:
                        valid = False
            if valid:
                bisect.insort(self.annotations, annot)

    def remove_annotations(self, *annotations: str) -> None:
        """Removes the specified annotation(s) from the DNA object.

        Parameters
        ----------
        annotations: str
            One or more names of annotations to remove

        Returns
        -------
        None

        Raises
        ------
        ValueError if one of the annotations can't be found
        """
        # TODO: Convert list to dictionary/set for annotations, should be faster

        for annot in annotations:
            if not self.has_annotation(annot):
                raise ValueError(f"Annotation {annot.name} not found!")
            
        for annot in annotations:
            for a in self.annotations:
                if annot == a.name:
                    self.annotations.remove(a)

    def from_annotation(self, annotation: str) -> DNA:
        """Returns a DNA sequence of a slice of the provided annotation.

        Parameters
        ----------
        annotation: str
            A string of the annotation provided to generate a slice of DNA
        
        Returns
        -------
        A DNA sequence created from the annotation provided

        Raises
        ------
        `ValueError` if the provided annotation is not a valid name of any annotation
        in self.
        """
        if annotation not in [annot.name for annot in self.annotations]:
            raise ValueError(f"Cannot find annotation: {annotation}")

        annot_span: tuple[int,int]
        for annot in self.annotations:
            if annot.name == annotation:
                return self[annot.span[0] : annot.span[1]]

    def concatenate(self, seq_list: list[DNA], name: str = "") -> DNA:
        """Concatenates two or more DNA sequences together, returning a new DNA sequence.

        Parameters
        ----------
        seq_list: list[DNA]
            A list of DNA to concatenate
        name: str
            The name of the new DNA molecule after concatenation

        Returns
        -------
        DNA resulting from concatenation

        Raises
        ------
        `InvalidSequence` error if sequences provided in seq_list are not all the same strandedness,
        if any sequence provided is circular, if not all are the same linearity, or not all the
        sequences provided are DNA.
        """
        # Validate all orientations/Bioproperties are the same
        if not (all(seq.is_double_stranded() for seq in seq_list) or \
                all(not seq.is_double_stranded() for seq in seq_list)):
            raise InvalidSequence("Not all provided sequences are the same strandedness!")

        # A circular sequence can't be concatenated. Where do we concatenate?
        if any(seq.circular == BioProperty.CIRCULAR for seq in seq_list):
            raise InvalidSequence("Can't concatenate a circular DNA sequence!")

        if not (all(seq.circular == BioProperty.LINEAR for seq in seq_list)):
            raise InvalidSequence("Not all provided sequences are the same circularity!")
        
        if not (all(seq.type == BioMolecule.DNA for seq in seq_list)):
            raise InvalidSequence("Not all provided sequences are DNA!")

        # Generate new DNA sequence
        new_seq = self.seq + "".join([sequence.seq for sequence in seq_list])

        # Generate new annotations, changing their indices as needed
        # TODO: Need to check if annotations are the same (by name probs) (not on slice, though)
        annotations = self.annotations.copy()
        offset = len(self.seq) # Used for calculating new offset in slices
        for seq in seq_list:
            if seq.annotations:
                for annot in seq.annotations:
                    new_span = (annot.span[0]+offset, annot.span[1]+offset)
                    if isinstance(annot, BioAnnotation):
                        annotations.append(BioAnnotation(new_span, annot.name, annot.orientation))
                    else:
                        annotations.append(Block(new_span, annot.name, annot.orientation, annot.pool))
            offset += len(seq)

        return DNA(
            new_seq, 
            name, 
            self.type, 
            self.circular, 
            self.strandedness, 
            annotations, 
            ((self.offsets[0][0], seq_list[-1].offsets[0][1]), (self.offsets[1][0], seq_list[-1].offsets[1][1]))
        )

    def to_oligo(self, slice: tuple[int, int] | None = None, reverse: bool = False) -> Oligo:
        """Converts a DNA piece to an oligo or a slice of the DNA if slice is provided.
        If reverse is True, takes the reverse complement of the sequence.

        Parameters
        ----------
        slice: tuple[int, int] | None (Optional)
            A tuple of integers indicating the slice of sequence to convert into
            a single-stranded Oligo
        reverse: bool (Default: False)
            A boolean of whether or not the sequence is the reverse complement of the
            sequence from slice

        Returns
        -------
        An Oligo of the sequence at the provided slice and of reverse orientation if 
        `reverse` is True, else forward
        """
        from biotools.oligo import Oligo
        if slice:
            # Validations
            if slice[0] == slice[1]:
                raise InvalidAnnotationException("Must not have same span indices!")

            # Validate tuple based on circularity
            if not self.circular:
                if slice[0] > slice[1]:
                    raise InvalidAnnotationException("First index must be less than the second index on a linear sequence!")

            # Validate indices
            if not all(sp >= 0 and sp <= self.length for sp in slice):
                raise InvalidAnnotationException("Spans integers must be valid!")

            if reverse:
                return Oligo(rev_comp(self.seq[slice[0]:slice[1]]), name=self.name+"_reverse_oligo")
            return Oligo(self.seq[slice[0]:slice[1]], name=self.name+"_oligo")

        if reverse:
            return Oligo(seq=rev_comp(self.seq), name=self.name+"_reverse_oligo")
        return Oligo(seq=self.seq, name=self.name+"_oligo")

    def print_annotations(self) -> None:
        """Prints annotations to visualize them sequentially"""
        # TODO: One thing could be something like the length of the annotation in bases is its width.
        # Maybe just print the annotation where it starts? And have its length in parentheses. Also, should
        # stack when needed.
        # Would be sick to have something like table of contents below the printed sequence so can identify 
        # annotation by number they are in table. In fact, this is what I'll do for all of them
        # Something like >>1>>, <<2222<<, etc. for fwd and rev, respectively and length 1 and 4
        # If there's ever overlapping annotations, they will appear stacked
        sorted_annots = self.annotations.copy()
        reprs = []

        for annot in sorted_annots:
            if annot.span[0] < annot.span[1]:
                length = annot.span[1] - annot.span[0]
            else:
                length = self.length - annot.span[0] + annot.span[1]

            if annot.orientation == BioOrientation.FORWARD:
                reprs.append(f">>>{annot.name} ({length} bp)>>>")
            else:
                reprs.append(f"<<<{annot.name} ({length} bp)<<<")

        if not self.is_circular():
            print(f"| {' --- '.join([r for r in reprs])} | ==>")
        else:
            print(f"| {' --- '.join([r for r in reprs])} | O")

    def reindex(self, index: int | str, inplace: bool = True) -> DNA | None:
        """Re-indexes the sequence so that index becomes index 0.
        NOTE: The default of this function is to perform an inplace reindex, but a new sequence may also
        be returned by setting inplace to False.

        Parameters
        ----------
        index: int | str
            If an integer, reindexes the sequence starting at that integer base.
            If a string, will reindex starting at the specified sequence.
        inplace: bool (Default: True)
            A boolean of whether or not to return the reindexed sequence as a new
            DNA object or perform the action in-place.

        Returns
        -------
        DNA of the reindexed sequence if inplace is True, None otherwise
        """
        useable_index: int = 0
        if isinstance(index, int):
            if index < 0 or index > self.length - 1:
                raise Exception("Invalid index provided. Must be within bounds of DNA sequence!")
            useable_index = index
            
        elif isinstance(index, str):
            f_matches = list(re.finditer(index, self.seq+self.seq[: len(index) - 1]))
            if len(f_matches) == 0:
                raise ValueError(f"Sequence {index} was not found in the DNA!")
            if len(f_matches) > 1:
                raise ValueError(f"Sequence {index} exists more than once! Must reindex to a unique sequence.")
            useable_index = f_matches[0].span()[0]
        
        if self.circular == BioProperty.LINEAR:
            raise Exception("Cannot re-index a linear fragment!")

        seq: str = ""
        if inplace: self.seq = self.seq[useable_index:] + self.seq[:useable_index]
        else: seq = self.seq[useable_index:] + self.seq[:useable_index]

        annotations: list[BioAnnotation] = []
        # Re-generate annotation indices
        for annot in self.annotations:
            start = annot.span[0] - useable_index
            stop = annot.span[1] - useable_index
            if start == 0:
                start = self.length
            elif start < 0:
                start += self.length

            if stop == 0:
                stop = self.length
            elif stop < 0:
                stop += self.length

            if inplace:
                annot.span = (start, stop)
            else:
                annotations.append(BioAnnotation((start, stop), annot.name, annot.orientation))

        if not inplace:
            return DNA(seq, self.name, self.type, self.circular, self.strandedness, annotations, self.offsets)

    def circularize(self) -> None:
        """Circularizes the sequence."""
        self.circular = BioProperty.CIRCULAR

    def get_overhangs(self) -> list[str]:
        """Returns the forward and reverse overhangs of the sequene as a list."""
        overhangs: list[str] = ["", ""]
        
        if self.strandedness != BioProperty.CUT:
            return overhangs
        
        if (diff := self.offsets[0][0] - self.offsets[1][0]) < 0:
            overhangs[0] = self.seq[0:-diff]
        elif diff != 0:
            overhangs[0] = self.bottom_strand[-diff:]

        if (diff := self.offsets[0][1] - self.offsets[1][1]) < 0:
            overhangs[1] = self.top_strand[(self.length + diff):]
        elif diff != 0:
            overhangs[1] = self.bottom_strand[:diff]
        return overhangs

    def insert(self, seq: DNA | str, index: int) -> None:
        """Inserts the provided sequence into self.seq at index.
        If seq is DNA, this will also add all annotations from seq to self.
        
        Parameters
        ----------
        seq: DNA | str
            The sequence to insert. Can be either a string or DNA of at least
            length one
        index: int
            The index location at which to insert seq
        
        Returns
        -------
        None
        """

        if len(seq) < 1:
            raise ValueError("Cannot insert an empty sequence!")
        
        if isinstance(seq, DNA):
            self.seq = self.seq[: index] + seq.seq + self.seq[index :]
        else:
            self.seq = self.seq[: index] + seq + self.seq[index :]

        # Add annotations from DNA, if it has any
        annotation_list: list[BioAnnotation] = []
        if isinstance(seq, DNA):
            for annot in seq.annotations:
                a = annot.copy()
                a.span = (a.span[0] + len(seq), a.span[1] + len(seq))
                annotation_list.append(a)
            if annotation_list: # Might not have any annotations
                self.add_annotations(*annotation_list)

        # Modify own annotations because we've inserted a sequence
        for annot in self.annotations:
            # index is within the annotation - expand annotation
            if index >= annot.span[0] and index < annot.span[1]:
                annot.span = (annot.span[0], annot.span[1] + len(seq))
            # index is before the annotation
            elif index < annot.span[0]:
                annot.span = (annot.span[0] + len(seq), annot.span[1] + len(seq))

    def find_sequence(self, seq: str) -> tuple[int, int] | list[tuple[int, int]] | None:
        """Finds the provided sequence within self.seq
        
        Parameters
        ----------
        seq: str
            The sequence to search for

        Returns
        -------
        A tuple of integers of the span or a list of such tuples if multiple of
        seq are present
        """

        matches = re.finditer(seq, self.seq)
        if matches:
            match_list: list[tuple[int, int]] = []
            for match in matches:
                match_list.append(match.span())
            return match_list[0] if len(match_list) == 1 else match_list
        return None
