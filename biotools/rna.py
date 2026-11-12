from __future__ import annotations
from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from biotools.oligo import Oligo
    from biotools.protein import Protein

from biotools.bio_alphabet import NTs
from biotools.bio_utils import rev_comp, validate_sequence, validate_annotations
from biotools.bio_enums import BioMolecule, BioProperty, BioOrientation
from biotools.bio_annotation import BioAnnotation, Block, reverse_annotations
from biotools.bio_exceptions import (
    InvalidSequence, InvalidInstantiationException, InvalidAnnotationException
)
from biotools.bio_alphabet import ANNOTATION_CHARS, RNA_CODON_AAs
import bisect
from IPython.display import display
import ipywidgets as widgets

class RNA:
    """
    Class representing the BioMolecule RNA. RNA is created by passing in a string comprised 
    of one or more of the characters defined in AAs

    Parameters
    ----------
    seq: str
        A string sequence of RNA. Must contain only valid nucleotide or degenerate
        base characters as specified in valid_chars
    type: BioMolecule
        An identifier Enum of the type RNA
    circular: BioProperty
        A boolean of whether or not the sequence is circular (a plasmid)
    strandedness: BioProperty
        A property stating whether or not the RNA is single- or double-stranded
    annotations: list[BioAnnotation]
        A list of annotations present on the RNA
    parent: RNA
        The RNA from which self was generated, if applicable
    """
    
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
                if not isinstance(annot, BioAnnotation) and not isinstance(annot, Block):
                    raise InvalidInstantiationException("Must use a valid BioAnnotation!")
        if annotations is not None:
            if validate_annotations(self, annotations):
                self.annotations = annotations
        else:
            self.annotations = []

        self._has_pool = False
        self._block_count = 0
        self._hash = None

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
            return f"{{ {self.name} | {self.type.value} | [{self.length} bp] | [{len(self.annotations)} Annotation(s)] | O }}"
        elif self.strandedness == BioProperty.DOUBLE_STRANDED:
            return f"{{ {self.name} | {self.type.value} | [{self.length} bp] | [{len(self.annotations)} Annotation(s)] | ==> }}"
        else:
            return f"{{ {self.name} | {self.type.value} | [{self.length} bp] | [{len(self.annotations)} Annotation(s)] | --> }}"

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
                        annotations.append(BioAnnotation(new_span, annot.name, annot.orientation, self))
                    elif isinstance(annot, Block):
                        annotations.append(Block(new_span, annot.name, annot.orientation, annot.pool, self))
            return RNA(
                self.seq[index],
                self.name,
                self.type,
                self.circular,
                self.strandedness,
                annotations
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
                                annotations.append(BioAnnotation(new_span, annot.name, annot.orientation, self))
                            elif isinstance(annot, Block):
                                annotations.append(Block(new_span, annot.name, annot.orientation, annot.pool, self))
                seq = self.seq[start:stop]
                if index.step is not None:
                    annotations = reverse_annotations(annotations, new_length)
                    seq = seq[::-1]
                return RNA(
                    seq,
                    self.name,
                    self.type,
                    self.circular,
                    self.strandedness,
                    annotations
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
                                    annotations.append(BioAnnotation(new_span, annot.name, annot.orientation, self))
                                elif isinstance(annot, Block):
                                    annotations.append(Block(new_span, annot.name, annot.orientation, annot.pool, self))
                    seq = self.seq[start:]+self.seq[:stop]
                    if index.step is not None:
                        annotations = reverse_annotations(annotations, new_length)
                        seq = seq[::-1]
                    return RNA(
                        seq,
                        self.name,
                        self.type,
                        BioProperty.LINEAR,
                        self.strandedness,
                        annotations
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
                                    annotations.append(BioAnnotation(new_span, annot.name, annot.orientation, self))
                                elif isinstance(annot, Block):
                                    annotations.append(Block(new_span, annot.name, annot.orientation, annot.pool, self))
                seq = self.seq[start:stop]
                if index.step is not None:
                    annotations = reverse_annotations(annotations, new_length)
                    seq = seq[::-1]
                return RNA(
                    seq,
                    self.name, self.type,
                    BioProperty.LINEAR,
                    self.strandedness,
                    annotations
                )

    def __hash__(self):
        """Hashes the sequence, defined as the addition of all hashes of
        the name, sequence, BioProperties, and annotations. Caches the hash as well
        for performance."""
        if self._hash is None:
            self._hash_tup = (
                self.seq, self.name, self.type.name, self.circular,
                self.strandedness, tuple(self.annotations)
            )
            self._hash = hash(self._hash_tup)
        else:
            if self._hash_tup != (
                self.seq, self.name, self.type.name, self.circular,
                self.strandedness, tuple(self.annotations)
            ):
                self._hash_tup = (
                    self.seq, self.name, self.type.name, self.circular,
                    self.strandedness, tuple(self.annotations)
                )
                self._hash = hash(self._hash_tup)
        return self._hash
    
    def __eq__(self, oth: RNA) -> bool:
        """Two RNA are equal if all their properties are equal."""
        if oth is self:
            return True
        if not isinstance(oth, RNA):
            return False
        if (self.seq, self.name, self.type, self.circular, self.strandedness, self.offsets, self.annotations) == \
            (oth.seq, oth.name, oth.type, oth.circular, oth.strandedness, oth.offsets, oth.annotations):
            return True
        if not ((self.seq, self.name, self.type, self.circular, self.strandedness, self.offsets, self.annotations) == \
                (oth.seq, oth.name, oth.type, oth.circular, oth.strandedness, oth.offsets, oth.annotations)) and \
                len(self) == len(oth):  # If they're the same length, could possibly be reindexes of each other
            if self.circular != oth.circular:
                return False
            elif not self.is_circular():  # Linear, so can't re-index
                return False
            return self._check_reindex_equality(oth)  # Same length with different annotations is same sequence
        return False
    
    def __add__(self, other):
        """Generates a new RNA molecule just taking the sequence. This works on linear and circular DNA."""
        return RNA(self.seq+other.seq, self.name, self.type, self.circular, self.strandedness, self.annotations)
    
    def __len__(self):
        return self.length

    def translate(self) -> Protein:
        """Translates self into a Protein."""
        if self.length % 3 != 0:
                    raise InvalidSequence(f"Cannot transcribe into Protein with {self.length} bases!")
        return Protein(seq="".join([RNA_CODON_AAs[self.seq[i:i+3]] for i in range(0, self.length, 3)]))
    
    def sequence(self) -> None:
        """Prints a string of the RNA sequence."""
        if self.strandedness == BioProperty.SINGLE_STRANDED:
            print(self.seq)
        else:
            print(self.seq)
            print(rev_comp(self.seq)[::-1], False)

    def copy(self) -> RNA:
        """Returns a copy of the RNA sequence."""
        return RNA(self.seq, self.name, self.type, self.circular, \
                   self.strandedness, self.annotations)
    
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
                    annotations.append(BioAnnotation(new_span, annot.name, annot.orientation, self))
            offset += len(seq)

        return RNA(new_seq, name, self.type, self.circular, self.strandedness, annotations)

    def print_annotations(self) -> None:
        """Prints annotations to visualize them by position on self.sequence().
        Renders a legend that maps chars to annotations."""
        reprs: str = ""
        annot_array: list[str] = []

        if len(self.annotations) > 77:
            print("WARNING!\nANNOTATIONS 78 AND ONWARD WILL NOT BE DISPLAYED!\n")  # TODO: View in web viewer instead!

        for i, annot in enumerate(self.annotations):
            if i > 76:  # No more than 77 (0-indexed) annotations can be printed
                break
            annot_array.append([*["X"] * self.length])
            if annot.span[0] < annot.span[1]:
                length = annot.span[1] - annot.span[0]
                for n in range(annot.span[0], annot.span[1]):
                    annot_array[i][n] = ANNOTATION_CHARS[i]
            else:  # Wraps around
                length = self.length - annot.span[0] + annot.span[1]
                for n in range(0, annot.span[1]):  # Beginning to annot.span[1]
                    annot_array[i][n] = ANNOTATION_CHARS[i]
                for n in range(annot.span[0], self.length):  # annot.span[0] to end
                    annot_array[i][n] = ANNOTATION_CHARS[i]

            if annot.is_forward():
                reprs += f"<text>{ANNOTATION_CHARS[i]}:  >>> {annot.name} ({length} bp) >>><br></text>"
            else:
                reprs += f"<text>{ANNOTATION_CHARS[i]}:  <<< {annot.name} ({length} bp) <<<<br></text>"

        self.sequence()

        max_depth: int = 0
        num_annots: int = len(annot_array)
        # Calculate max depth
        for i in range(self.length):
            pos_max_depth: int = 0
            for n in range(num_annots):
                if annot_array[n][i] != "X":
                    pos_max_depth += 1
            if pos_max_depth > max_depth:
                max_depth = pos_max_depth

        # Collapse columns together based on max_depth
        annot_strings: list[str] = [" " * self.length for n in range(max_depth)]
        span_list: list[list[tuple[int, int]]] = [[] for n in range(max_depth)]
        for i, annot in enumerate(self.annotations):
            annot_span: tuple[int, int] = annot.span
            for j, spans in enumerate(span_list):
                overlap: bool = False
                if spans:  # We have spans to check
                    for span in spans:
                        if (annot_span[0] > span[0] and annot_span[0] < span[1]) or \
                            (annot_span[1] > span[0] and annot_span[1] < span[1]):  # Overlap
                            overlap = True
                            break
                    if not overlap:
                        if annot_span[0] < annot_span[1]:
                            annot_strings[j] = annot_strings[j][: span_list[j][-1][1]] + annot_strings[j][span_list[j][-1][1] : annot_span[0]] + \
                                            ANNOTATION_CHARS[i] * annot.length + annot_strings[j][annot_span[1] :]
                        else:  # wrap-around annotation
                            annot_strings[j] = ANNOTATION_CHARS[i] * annot_span[1] + annot_strings[j][annot_span[1] : span_list[j][-1][1]] + \
                                                annot_strings[j][span_list[j][-1][1] : annot_span[0]] + ANNOTATION_CHARS[i] * (annot.parent.length - annot_span[0])
                        spans.append(annot_span)
                        break
                else:  # First span added in that position
                    if annot_span[0] < annot_span[1]:
                        annot_strings[j] = annot_strings[j][: annot_span[0]] + \
                                            ANNOTATION_CHARS[i] * annot.length + annot_strings[j][annot_span[1] :]
                    else:  # wrap-around annotation
                        annot_strings[j] = ANNOTATION_CHARS[i] * annot_span[1] + annot_strings[j][annot_span[1] : annot_span[0]] + \
                                            ANNOTATION_CHARS[i] * (annot.parent.length - annot_span[0])
                    spans.append(annot_span)
                    break

        # Create and print annotation strings under sequence
        for annot_string in annot_strings:
            print(annot_string)
        print()

        display(
            widgets.HTML(f"""
            <div style="
                background: #1e1e1e;
                color: #f0f0f0;
                border: 1px solid #444;
                border-radius: 8px;
                padding: 12px;
                max-height: 300px;
                overflow-y: auto;
                font-family: Arial, sans-serif;">
            {reprs}
            </div>"""
            )
        )

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
    
