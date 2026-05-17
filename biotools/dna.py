from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from oligo import Oligo

from bio_alphabet import NTs
from bio_utils import rev_comp, validate_sequence, validate_annotations
from bio_enums import *
from bio_annotation import *
from bio_exceptions import *
import bisect

class DNA():
    """
    Class representing DNA.
    
    DNA is a string comprised of one or more of the characters defined in
    valid_chars.

    DNA is a superclass from which more specific
    types of DNA such as Amplicon, Plasmid, Oligo, etc. inherit.

    Parameters
    ----------
    seq
        A string sequence of DNA. Must contain only valid nucleotide or degenerate base
        characters as specified in valid_chars

    type
        An identifier string of the type of DNA. Mostly used to specify
        classes that inherit 

    circular
        A boolean of whether or not the sequence is circular (a plasmid)
    """

    # TODO: ARE THERE ANY LIBRARIES THAT COMBINE SEQUENCE STUFF, MELTING TEMP, CLONING ABILITY (COMBINATORIAL),
    # IN ONE LIBRARY?
    # TODO: Reactions and their products are stored in something

    def __init__(
        self, 
        seq: str, 
        name: str = "", 
        type: str = BioMolecule.DNA,
        circular: str = BioProperty.LINEAR,
        strandedness: str = BioProperty.DOUBLE_STRANDED,
        annotations: list[BioAnnotation] | None = None,
        offsets: tuple[tuple[int,int], tuple[int,int]] | None = None,
        parent: DNA = None
    ):
        """Default constructor."""
        if validate_sequence(seq, type):
            self.seq = seq.upper() # NOTE: I want capitalized sequences
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
                if not isinstance(annot, BioAnnotation):
                    raise InvalidInstantiationException("Must use a valid BioAnnotation!")
        if annotations is not None:
            if validate_annotations(self, annotations):
                self.annotations = sorted(annotations)
        else:
            self.annotations = []
        """
        If offsets is None, then we assume that this is the original double- or single-stranded DNA
        sequence. If offsets is a tuple of ((top_5', top_3'), (bottom_3', bottom_5')), then
        we assume that this was created from a reaction and contains overhangs.
        Can grab both the exact way to print and also figure out the overhang of each using these tuples.
        """
        if offsets is not None:
            self.offsets = offsets
        else:
            self.offsets = None
        self.parent = parent

    @property
    def length(self) -> int:
        """The length of the DNA sequence."""
        return len(self.seq)
    
    @property
    def top_strand(self) -> str:
        """The top strand shown 5' -> 3'."""
        # TODO: Update to show correct top and bottom strands if cut w/ offsets
        return self.seq
    
    @property
    def bottom_strand(self) -> str:
        """The bottom strand shown 5' -> 3'."""
        if self.strandedness == BioProperty.SINGLE_STRANDED:
            raise Exception("Cannot return bottom strand of a single-stranded sequence!")
        else:
            return rev_comp(self.seq)

    def __repr__(self):
        if self.circular == BioProperty.CIRCULAR:
            return f"{self.name} | {self.type.value} | [{self.length} bp] | [{len(self.annotations)} Annotation(s)] | O"
        elif self.strandedness == BioProperty.DOUBLE_STRANDED:
            return f"{self.name} | {self.type.value} | [{self.length} bp] | [{len(self.annotations)} Annotation(s)] | ==>"
        elif self.strandedness == BioProperty.CUT:
            return f"{self.name} | {self.type.value} | [{self.length} bp] | [{len(self.annotations)} Annotation(s)] | Cut"
        else:
            return f"{self.name} | {self.type.value} | [{self.length} bp] | [{len(self.annotations)} Annotation(s)] | -->"
    
    def __getitem__(self, index): # TODO: Generate reverses as well?? [::-1] Might not be as hard as I think, can calculate everything and invert in some clever way
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
            return DNA(self.seq[index], self.name, self.type, self.circular, self.strandedness, annotations, self.offsets, self)
        elif isinstance(index, slice):
            start = index.start if index.start is not None else 0
            if start < 0:
                if abs(start) > self.length:
                    raise IndexError("Index provided is out of bounds of sequence length!")
                start = len(self.seq) + start
            stop = index.stop if index.stop is not None else len(self.seq)
            new_length = stop - start
            # Check orientation, and continue as appropriate
            if self.circular == BioProperty.LINEAR:
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
                return DNA(self.seq[start:stop], self.name, self.type, self.circular, self.strandedness, annotations, self.offsets, self)
            else: # Circular
                if start < 0 or stop > self.length or stop < 0 or start > self.length:
                    raise IndexError("Cannot have indices beyond bounds of sequence length!")
                if start > stop:
                    new_length = (self.length - start) + stop
                    for annot in self.annotations:
                        if annot.span[0] >= start: # Span starts within slice bounds of circular sequence
                            if annot.span[1] <= stop:
                                new_span = (annot.span[0] - start, (self.length - start) + annot.span[1])
                            elif annot.span[1] > stop and annot.span[1] > start:
                                new_span = (annot.span[0] - start, annot.span[1] - start)
                            else: # Stop must be < length of seq
                                new_span = (annot.span[0] - start, new_length)
                        else: # Span starts from start (0) (is less than starting index) NOTE: Here we have to check if the spans are the same because of re-indexing weirdness
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
                    return DNA(self.seq[start:]+self.seq[:stop], self.name, self.type, self.circular, self.strandedness, annotations, self)
                else: # Regular workflow because the slice is like a linear fragment (start < stop)
                    for annot in self.annotations:
                        new_span = None
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
                return DNA(self.seq[start:stop], self.name, self.type, self.circular, self.strandedness, annotations, self.offsets, self)

    def __hash__(self):
        # TODO: json.dump() hash, need a better one for gibson reaction
        return hash(self.seq)
    
    def __eq__(self, other):
        if other is self:
            return True
        return (self.seq, self.name, self.type, self.circular, self.strandedness, self.annotations) == \
               (other.seq, other.name, other.type, other.circular, other.strandedness, other.annotations)
    
    def __add__(self, other):
        """Generates a new DNA molecule just taking the sequence. This works on linear and circular DNA."""
        return DNA(self.seq+other.seq, self.name, self.type, self.circular, self.strandedness, self.annotations+other.annotations, self.offsets)

    def __len__(self):
        return self.length
    
    def __lt__(self, oth):
        if not isinstance(oth, DNA):
            raise Exception("Cannot compare between different objects!")
        return hash(self.seq) < hash(oth.seq)

    def sequence(self) -> None:
        """Prints a string of the DNA sequence."""
        if self.strandedness == BioProperty.SINGLE_STRANDED:
            print(self.seq) # TODO: Also make this cut by offsets?
        elif self.strandedness == BioProperty.CUT:
            # NOTE: Products can probably only be double-stranded? Should I check this?
            print(self.offsets[0][0]*" "+self.seq[self.offsets[0][0]:self.length - self.offsets[0][1]])
            bottom = rev_comp(self.seq)[::-1]
            print(self.offsets[1][0]*" "+bottom[self.offsets[1][0]:len(bottom) - self.offsets[1][1]])
        else:
            print(self.seq)
            print(rev_comp(self.seq)[::-1])

    def info(self) -> None:
        print("Name: ", self.name)
        print("Sequence: ", self.seq)
        print("Reverse Complement: ", rev_comp(self.seq))
        print("Length: ", self.length)
        print("Type: ", self.type.value)
        print(f"Circular: {self.circular == BioProperty.CIRCULAR}")
        print(f"Strandedness: {self.strandedness.value}")

    def copy(self) -> DNA:
        """Returns a copy of the DNA sequence."""
        return DNA(self.seq, self.name, self.type, self.circular, self.strandedness, self.offsets)
    
    def rev_comp(self) -> DNA:
        """Returns a DNA copy of the reverse complement of the sequence."""
        if not self.is_double_stranded():
            raise Exception("Cannot generate reverse complement of single-stranded molecule!")
        else:
            annots: list[BioAnnotation] = []
            # Generate inverted annotations 
            for a in self.annotations:
                new_span = (self.length - a.span[1], self.length - a.span[0])
                orient = BioOrientation.FORWARD if a.orientation == BioOrientation.REVERSE else BioOrientation.REVERSE
                if isinstance(a, BioAnnotation):
                    annots.append(BioAnnotation(new_span, a.name, orient))
                elif isinstance(a, Block):
                    annots.append(Block(new_span, a.name, orient, a.pool))
            return DNA(rev_comp(self.seq), self.name+"_revcomp", self.type, self.circular, self.strandedness, annots, self.offsets)
        
    def is_circular(self):
        return self.circular == BioProperty.CIRCULAR
    
    def is_double_stranded(self):
        return self.strandedness == BioProperty.DOUBLE_STRANDED
    
    def is_cut(self):
        return self.strandedness == BioProperty.CUT
        
    def add_annotation(self, annotation: BioAnnotation) -> None:
        """Adds an annotation to the DNA sequence."""
        # Spans can't be same
        if annotation.span[0] == annotation.span[1]:
            raise InvalidAnnotationException("Must not have same span indices!")

        # Validate tuple based on circularity
        if self.circular == BioProperty.LINEAR:
            if annotation.span[0] > annotation.span[1]:
                raise InvalidAnnotationException("First index must be less than the second index on a linear sequence!")

        # Validate indices
        if not all(sp >= 0 and sp <= self.length for sp in annotation.span):
            raise InvalidAnnotationException("Span integers must be valid!")

        if annotation not in self.annotations:
            bisect.insort(self.annotations, annotation)

    def concatenate(self, seq_list: list[DNA], name: str = "") -> DNA:
        """Concatenates two or more DNA sequences together, returning a new DNA sequence."""

        # Validate all orientations/Bioproperties are the same
        if not (all(seq.strandedness == BioProperty.DOUBLE_STRANDED for seq in seq_list) or \
                all(seq.strandedness == BioProperty.SINGLE_STRANDED for seq in seq_list)):
            raise InvalidSequence("Not all provided sequences are the same strandedness!")
        
        # A circular sequence can't be concatenated. Where do we concatenate?
        if any(seq.circular == BioProperty.CIRCULAR for seq in seq_list):
            raise InvalidSequence("Can't concatenate a circular DNA sequence!")

        if not (all(seq.circular == BioProperty.LINEAR for seq in seq_list)):
            raise InvalidSequence("Not all provided sequences are the same circularity!")
        
        if not (all(seq.type == BioMolecule.DNA for seq in seq_list)):
            raise InvalidSequence("Not all provided sequences are DNA!")

        # Generate new DNA sequence
        new_seq = self.seq + "".join([seq.top_strand for seq in seq_list])

        # Generate new annotations, changing their indices as needed
        annotations = self.annotations
        offset = len(self.seq) # Used for calculating new offset in slices
        for seq in seq_list:
            if seq.annotations:
                for annot in seq.annotations:
                    new_span = (annot.span[0]+offset, annot.span[1]+offset)
                    annotations.append(BioAnnotation(new_span, annot.name, annot.orientation))
            offset += len(seq)

        return DNA(new_seq, name, self.type, self.circular, self.strandedness, annotations, self.offsets)

    def to_oligo(self, slice: tuple[int, int] | None = None, reverse: bool = False) -> Oligo:
        """Converts a DNA piece to an oligo or a slice of the DNA if slice is provided."""
        from oligo import Oligo
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

    # TODO: Improve visualization for this
    def print_annotations(self) -> None:
        """Prints annotations to visualize them sequentially"""
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

        if self.circular == BioProperty.LINEAR:
            print(f"| {" --- ".join([r for r in reprs])} | ==>")
        else:
            print(f"| {" --- ".join([r for r in reprs])} | O")

    def reindex(self, index: int, inplace: bool = True) -> DNA | None:
        """Re-indexes the sequence so that index becomes index 0.
        NOTE: The default of this function is to perform an inplace reindex, but a new sequence may also
        be returned by setting inplace to False."""
        if index < 0 or index > self.length - 1:
            raise Exception("Invalid index provided. Must be within bounds of DNA sequence!")
        
        if self.circular == BioProperty.LINEAR:
            raise Exception("Cannot re-index a linear fragment!")

        seq: str = ""
        if inplace: self.seq = self.seq[index:] + self.seq[:index]
        else: seq = self.seq[index:] + self.seq[:index]

        annotations: list[BioAnnotation] = []
        # Re-generate annotation indices
        for annot in self.annotations:
            start = annot.span[0] - index
            stop = annot.span[1] - index
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

    # TODO: Do we infer sequence based on numpy array? Is there a fast,
    # memory-efficient way of doing this in python?
