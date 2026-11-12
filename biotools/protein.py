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
from biotools.bio_alphabet import ANNOTATION_CHARS
import bisect
from IPython.display import display
import ipywidgets as widgets

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
        annotations: list[BioAnnotation] | None = None
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

        self._has_pool = False
        self._block_count = 0
        self._hash = None

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
                        annotations.append(BioAnnotation(new_span, annot.name, annot.orientation, self))
                    elif isinstance(annot, Block):
                        annotations.append(Block(new_span, annot.name, annot.orientation, annot.pool, self))
            return Protein(
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
            return Protein(
                seq,
                self.name,
                self.type,
                self.circular,
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
    
    def __eq__(self, oth: Protein) -> bool:
        """Two Proteins are equal if all their properties are equal."""
        if oth is self:
            return True
        if not isinstance(oth, Protein):
            return False
        if (self.seq, self.name, self.type, self.circular, self.strandedness, self.annotations) == \
            (oth.seq, oth.name, oth.type, oth.circular, oth.strandedness, oth.annotations):
            return True
        if not ((self.seq, self.name, self.type, self.circular, self.strandedness, self.annotations) == \
                (oth.seq, oth.name, oth.type, oth.circular, oth.strandedness, oth.annotations)) and \
                len(self) == len(oth):  # If they're the same length, could possibly be reindexes of each other
            if self.circular != oth.circular:
                return False
        return False

    def extinction_coefficient(self) -> float:
        """Returns the extinction coefficient calculated using the 
        Gill and von Hippel formula:
        
        EC = (n_W * 5500) + (n_Y * 1490) + (n_C * 125)
        n_W = Number of trytophans in the protein
        n_Y = Number of tyrosines in the protein
        n_C = Number of cystines in the protein
        """
        return self.seq.count("W") * 5500 + \
               self.seq.count("Y") * 1490 + self.seq.count("C") * 125

    def copy(self) -> Protein:
        """Returns a copy of the RNA sequence."""
        return Protein(self.seq, self.name, self.type, self.circular, \
                       self.strandedness, self.annotations)

    def __len__(self):
        return self.length

    def sequence(self) -> str:
        """Prints a string of the Protein sequence."""
        return self.seq

    def add_annotations(self, *annotations: BioAnnotation | Block) -> None:
        """Adds an annotation or annotations to the DNA sequence.

        Parameters
        ----------
        annotation: BioAnnotation | Block
            One or more annotations to add to self
        """
        # Spans can't be same
        for annot in annotations:
            if annot.span[0] == annot.span[1]:
                raise InvalidAnnotationException("Must not have same span indices!")

        for annot in annotations:
            if annot.span[0] > annot.span[1]:
                raise InvalidAnnotationException("First index must be less than the second index on a linear sequence!")

        # Validate indices
        for annot in annotations:
            if not all(sp >= 0 and sp <= self.length for sp in annot.span):
                raise InvalidAnnotationException("Span integers must be valid!")

        # Check valid name because of from_annotation()
        for annot in annotations:
            if annot.name in [a.name for a in self.annotations]:
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

    def concatenate(self, seq_list: list[Protein], name: str = "") -> Protein:
        """Concatenates two or more Protein sequences together, returning a new Protein sequence."""
        
        # A circular sequence can't be concatenated. Where do we concatenate?
        if any(seq.circular == BioProperty.CIRCULAR for seq in seq_list):
            raise InvalidSequence("Can't concatenate a circular Protein sequence!")

        if not (all(seq.circular == BioProperty.LINEAR for seq in seq_list)):
            raise InvalidSequence("Not all provided sequences are the same circularity!")
        
        if not (all(seq.type == BioMolecule.Protein for seq in seq_list)):
            raise InvalidSequence("Not all provided sequences are RNA!")

        # Generate new RNA sequence
        new_seq = self.seq + "".join([seq.top_strand for seq in seq_list])

        # Generate new annotations, changing their indices as needed
        annotations = self.annotations
        offset = len(self.seq) # Used for calculating new offset in slices
        for seq in seq_list:
            if seq.annotations:
                for annot in seq.annotations:
                    new_span = (annot.span[0] + offset, annot.span[1] + offset)
                    annotations.append(BioAnnotation(new_span, annot.name, annot.orientation, self))
            offset += len(seq)

        return Protein(new_seq, name, self.type, self.circular, self.strandedness, annotations)

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
