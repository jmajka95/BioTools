from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from oligo import Oligo
    from dna import DNA

from protein import Protein
from bio_enums import BioMolecule, BioProperty
from bio_annotation import *
from bio_alphabet import *
from bio_utils import rev_comp, validate_sequence, validate_annotations

class RNA:
    """Class representing RNA."""

    # TODO: Make a Nucleotide superclass? I will probably be re-making some stuff that is consistent between? Maybe doesn't matter
    
    def __init__(
        self, 
        seq: str, 
        name: str = "", 
        type: str = BioMolecule.RNA,
        circular: str = BioProperty.LINEAR,
        strandedness: str = BioProperty.SINGLE_STRANDED,
        annotations: list[BioAnnotation] | None = None
    ):
        """Default constructor."""
        if validate_sequence(seq, type):
            self.seq = seq.upper() # NOTE: I want capitalized sequences
        else:
            raise InvalidSequence("Must provide a valid RNA sequence!")
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
                self.annotations = annotations
        else:
            self.annotations = []

    @property
    def length(self) -> int:
        """The length of the RNA sequence."""
        return len(self.seq)

    @property
    def top_strand(self) -> str:
        """The top strand shown 5' -> 3'."""
        return self.seq

    @property
    def bottom_strand(self) -> str:
        """The bottom strand shown 5' -> 3'."""
        if self.strandedness == BioProperty.SINGLE_STRANDED:
            raise Exception("Cannot return bottom strand of a single-stranded sequence!")
        else:
            return rev_comp(self.seq, False)

    def __repr__(self):
        if self.circular == BioProperty.CIRCULAR: 
            return f"{self.name} | {self.type.value} | [{self.length} bp] | [{len(self.annotations)} Annotation(s)] | O"
        elif self.strandedness == BioProperty.DOUBLE_STRANDED:
            return f"{self.name} | {self.type.value} | [{self.length} bp] | [{len(self.annotations)} Annotation(s)] | ==>"
        else:
            return f"{self.name} | {self.type.value} | [{self.length} bp] | [{len(self.annotations)} Annotation(s)] | -->"

    def __getitem__(self, index): # TODO: Generate reverses as well?? [::-1] Might not be as hard as I think, can calculate everything and invert in some clever way
        annotations: list[BioAnnotation] = []
        if isinstance(index, int): # Single integer provided
            if index < 0 or index > self.length - 1:
                raise Exception("Integer provided is out of bounds of sequence length!")
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
            return RNA(self.seq[index], self.name, self.type, self.circular, self.strandedness, annotations)
        elif isinstance(index, slice):
            start = index.start if index.start is not None else 0
            stop = index.stop if index.stop is not None else len(self.seq)
            new_length = stop - start
            # Check orientation, and continue as appropriate
            if self.circular == BioProperty.LINEAR:
                if start > stop:
                    raise Exception("For linear RNA, first index must be less than the second!")
                elif start < 0 or stop > len(self.seq):
                    raise Exception("Cannot have indices beyond bounds of sequence length!")
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
                return RNA(self.seq[start:stop], self.name, self.type, self.circular, self.strandedness, annotations)
            else: # Circular
                if start < 0 or stop > self.length or stop < 0 or start > self.length:
                    raise Exception("Cannot have indices beyond bounds of sequence length!")
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
                    return RNA(self.seq[start:]+self.seq[:stop], self.name, self.type, self.circular, self.strandedness, annotations)
                # Regular workflow because the slice is like a linear fragment (start < stop)
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
                            else:
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
                return RNA(self.seq[start:stop], self.name, self.type, self.circular, self.strandedness, annotations)

    def __hash__(self):
        # TODO: json.dump() hash of specific stuff related to this. Can use for equality
        return hash(self.seq)
    
    def __eq__(self, other):
        if other is self:
            return True
        return (self.seq, self.name, self.type, self.circular, self.strandedness, self.annotations) == \
               (other.seq, other.name, other.type, other.circular, other.strandedness, other.annotations)
    
    def __add__(self, other):
        """Generates a new RNA molecule just taking the sequence. This works on linear and circular DNA."""
        return RNA(self.seq+other.seq, self.name, self.type, self.circular, self.strandedness, self.annotations)
    
    def __len__(self):
        return self.length
    
    def sequence(self) -> None:
        """Prints a string of the RNA sequence."""
        if self.strandedness == BioProperty.SINGLE_STRANDED:
            print(self.seq)
        else:
            print(self.seq)
            print(rev_comp(self.seq)[::-1], False)

    def info(self) -> None:
        print("Name: ", self.name)
        print("Sequence: ", self.seq)
        print("Reverse Complement: ", rev_comp(self.seq, False))
        print("Length: ", self.length)
        print("Type: ", self.type.value)
        print(f"Circular: {self.circular == BioProperty.CIRCULAR:}")
        print(f"Strandedness: {self.strandedness.value}")

    def copy(self) -> RNA:
        """Returns a copy of the RNA sequence."""
        return RNA(self.seq, self.name, self.type, self.circular, self.strandedness)
    
    def rev_comp(self) -> RNA:
        """Returns an RNA copy of the reverse complement of the sequence."""
        if self.strandedness != BioProperty.DOUBLE_STRANDED:
            raise Exception("Cannot generate reverse complement of single-stranded molecule!")
        else:
            return RNA(rev_comp(self.seq, False), self.name+"_revcomp")

    def add_annotation(self, annotation: BioAnnotation) -> None:
        """Adds an annotation to the RNA sequence."""
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
            self.annotations.append(annotation)

    def concatenate(self, seq_list: list[RNA], name: str = "") -> RNA:
        """Concatenates two or more RNA sequences together, returning a new RNA sequence."""

        # Validate all orientations/Bioproperties are the same
        if not (all(seq.strandedness == BioProperty.DOUBLE_STRANDED for seq in seq_list) or \
                all(seq.strandedness == BioProperty.SINGLE_STRANDED for seq in seq_list)):
            raise InvalidSequence("Not all provided sequences are the same strandedness!")
        
        # A circular sequence can't be concatenated. Where do we concatenate?
        if any(seq.circular == BioProperty.CIRCULAR for seq in seq_list):
            raise InvalidSequence("Can't concatenate a circular RNA sequence!")

        if not (all(seq.circular == BioProperty.LINEAR for seq in seq_list)):
            raise InvalidSequence("Not all provided sequences are the same circularity!")
        
        if not (all(seq.type == BioMolecule.RNA for seq in seq_list)):
            raise InvalidSequence("Not all provided sequences are RNA!")

        # Generate new RNA sequence
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

        return RNA(new_seq, name, self.type, self.circular, self.strandedness, annotations)

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

    def reindex(self, index: int) -> None:
        """Re-indexes the sequence so that index becomes index 0."""
        if index < 0 or index > self.length - 1:
            raise Exception("Invalid index provided. Must be within bounds of RNA sequence!")
        
        if self.circular == BioProperty.LINEAR:
            raise Exception("Cannot re-index a linear fragment!")

        self.seq = self.seq[index:] + self.seq[:index]

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

            annot.span = (start, stop)

    def circularize(self) -> None:
        """Circularizes the sequence."""
        self.circular = BioProperty.CIRCULAR
    
