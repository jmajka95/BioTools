from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from biotools.oligo import Oligo

from biotools.bio_utils import rev_comp, validate_sequence, validate_annotations
from biotools.bio_enums import BioMolecule, BioProperty, BioOrientation
from biotools.bio_annotation import BioAnnotation, Block, reverse_annotations
from biotools.bio_exceptions import (
    InvalidSequence, InvalidInstantiationException, InvalidAnnotationException
)
import bisect

class Protein:
    """
    Class representing the BioMolecule Protein. Protein is created by passing in a string comprised 
    of one or more of the characters defined in AAs

    Parameters
    ----------
    seq: str
        A string sequence of Protein. Must contain only valid nucleotide or degenerate
        base characters as specified in valid_chars
    type: BioMolecule
        An identifier Enum of the type Protein
    circular: BioProperty
        A boolean of whether or not the sequence is circular (a plasmid)
    strandedness: BioProperty
        A property stating whether or not the Protein is single- or double-stranded
    annotations: list[BioAnnotation]
        A list of annotations present on the Protein
    parent: Protein
        The Protein from which self was generated, if applicable
    """
    
    def __init__(
        self, 
        seq: str, 
        name: str = "", 
        type: str = BioMolecule.PROTEIN,
        circular: str = BioProperty.LINEAR,
        strandedness: str = BioProperty.SINGLE_STRANDED,
        annotations: list[BioAnnotation] | None = None,
        parent: Protein | None = None
    ):
        """Default constructor."""
        if validate_sequence(seq, type):
            self.seq = seq.upper() # NOTE: I want capitalized sequences
        else:
            raise InvalidSequence("Must provide a valid Protein sequence!")
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
                self.annotations = annotations
        else:
            self.annotations = []
        self.parent = parent if parent is not None else None

    @property
    def length(self) -> int:
        """The length of the Protein sequence."""
        return len(self.seq)

    def __repr__(self):
        return f"{{ {self.name} | {self.type.value} | [{self.length} amino acid(s)] | [{len(self.annotations)} Annotation(s)] }}"

    def __getitem__(self, index):
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
            return Protein(
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
                return Protein(
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
                    seq = self.seq[start:]+self.seq[:stop]
                    if index.step is not None:
                        annotations = reverse_annotations(annotations, new_length)
                        seq = seq[::-1]
                    return Protein(
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
                return Protein(
                    seq, 
                    self.name, self.type, 
                    BioProperty.LINEAR, 
                    self.strandedness, 
                    annotations, 
                    self.offsets, 
                    self
                )

    def __hash__(self):
        # TODO: json.dump() hash
        return hash(self.seq)
    
    def __eq__(self, other):
        if other is self:
            return True
        return (self.seq, self.name, self.type, self.circular, self.strandedness, self.annotations, self.parent) == \
               (other.seq, other.name, other.type, other.circular, other.strandedness, other.annotations, self.parent)

    def __len__(self):
        return self.length

    def sequence(self) -> str:
        """Prints a string of the Protein sequence."""
        return self.seq

    