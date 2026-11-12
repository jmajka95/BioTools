from __future__ import annotations

from biotools.bio_enums import BioMolecule, BioProperty
from biotools.bio_annotation import BioAnnotation, Block, reverse_annotations
from biotools.bio_exceptions import (
    InvalidAnnotationException, InvalidSequence, InvalidInstantiationException
)
from biotools.bio_utils import validate_sequence, validate_annotations, rev_comp
from biotools.bio_alphabet import ANNOTATION_CHARS

import json
import requests  # type: ignore
import bisect
from itertools import product
from IPython.display import display
import ipywidgets as widgets

from polymerase_tm.core import tm

class Oligo():
    """Class representing single-stranded oligos, typically used as primers.
    Oligos may also be annealed to form double-stranded DNA with or without
    overhangs.

    Parameters
    ----------
    seq: str
        A string sequence of DNA. Must contain only valid nucleotide or degenerate base
        characters as specified in valid_chars
    type: BioMolecule
        An identifier Enum of the type Oligo
    circular: BioProperty
        A boolean of whether or not the sequence is circular (a plasmid)
    strandedness: BioProperty
        A property stating whether or not the Oligo is single- or double-stranded
    annotations: list[BioAnnotation | Block]
        A list of annotations present on the Oligo
    """

    def __init__(
        self, 
        seq: str, 
        name: str = "",
        type: BioMolecule = BioMolecule.OLIGO,
        circular: BioProperty = BioProperty.LINEAR,
        strandedness: BioProperty = BioProperty.SINGLE_STRANDED,
        annotations: list[BioAnnotation | Block] | None = None
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

        self._has_pool = False
        self._block_count = 0
        self._hash = None

    @property
    def length(self) -> int:
        """The length of the DNA sequence."""
        return len(self.seq)

    @property
    def tm(self) -> int:
        return tm(self.seq)

    def __hash__(self) -> int:
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

    def __eq__(self, oth: Oligo) -> bool:
        """Two Oligo are equal if all their properties are equal."""
        if oth is self:
            return True
        if not isinstance(oth, Oligo):
            return False
        if (self.seq, self.name, self.type, self.circular, self.strandedness, self.annotations) == \
            (oth.seq, oth.name, oth.type, oth.circular, oth.strandedness, oth.annotations):
            return True
        if ((self.seq, self.name, self.type, self.circular, self.strandedness, self.annotations) == \
                (oth.seq, oth.name, oth.type, oth.circular, oth.strandedness, oth.annotations)):
            return True
        return False

    def __repr__(self) -> str:
        """Representation of an Oligo, generically as name, BioMolecule, length, annotation number, 
        and strandedness."""
        return f"{{ {self.name} | {self.type.value} | [{self.length} bp] | [{len(self.annotations)} Annotation(s)] | --> }}"

    def __add__(self, oth: Oligo) -> Oligo:
        """Generates a new DNA molecule just taking the sequence. 
        This works on linear and circular DNA.
        """
        return Oligo(
            self.seq + oth.seq,
            self.name,
            self.type,
            self.circular,
            self.strandedness,
            self.annotations + oth.annotations
        )

    def __len__(self) -> int:
        return self.length

    def __lt__(self, oth: Oligo) -> bool:
        """Oligos are compared via their hash."""
        if not isinstance(oth, Oligo):
            raise Exception(f"Cannot compare {type(self).__name__} with {type(oth).__name__}!")
        return hash(self) < hash(oth)

    def __getitem__(self, index) -> slice:
        """Grabs and returns a slice of an Oligo, preserving annotations and re-indexing them
        as appropriate.
        NOTE: Only reversing (-1) is supported as a step when slicing. All other values
        raise a ValueError.
        """
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
            return Oligo(
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
                raise IndexError("For Oligos, first index must be less than the second!")
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
            return Oligo(
                seq,
                self.name,
                self.type,
                self.circular,
                self.strandedness,
                annotations
            )

    def has_pool(self) -> bool:
        """Whether or not the Oligo has a pool"""
        for annot in self.annotations:
            if isinstance(annot, Block):
                return True
        return False

    def copy(self) -> Oligo:
        """Returns a copy of self."""
        return Oligo(
            self.seq, self.name, self.type, self.circular, 
            self.strandedness, self.annotations
        )

    def concatenate(self, seq_list: list[Oligo], name: str = "") -> Oligo:
        """Concatenates two or more Oligo sequences together, returning a new Oligo sequence."""

        # A circular sequence can't be concatenated. Where do we concatenate?
        if any(seq.circular == BioProperty.CIRCULAR for seq in seq_list):
            raise InvalidSequence("Can't concatenate a circular Oligo sequence!")

        if not (all(seq.circular == BioProperty.LINEAR for seq in seq_list)):
            raise InvalidSequence("Not all provided sequences are the same circularity!")
        
        if not (all(seq.type == BioMolecule.Oligo for seq in seq_list)):
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

        return Oligo(new_seq, name, self.type, self.circular, self.strandedness, annotations)

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

        print(self.seq)

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

    def get_pools(self) -> list[Oligo] | None:
        """Returns all sequences generated from Blocks. If more than one Block exists,
        returns all possible combinations of sequences generated.

        Returns
        -------
        list[DNA] if self contains pools, otherwise None
        """
        if not self.has_pool():
            return None

        next_idx: int = 0
        self_seqs: list[str] = []
        groups: list[list[Oligo]] = []

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
        seqs: list[Oligo] = []
        block_combos = product(*groups)
        for combo in block_combos:
            new_seq: str = ""
            for i in range(len(combo)):
                new_seq += self_seqs[i]
                new_seq += combo[i].seq
            new_seq += self_seqs[len(combo)]
            seqs.append(
                Oligo(seq=new_seq)
            )

        seqs.append(self.copy())  # Add our own seq because it's not captured in the Block combinations
        return seqs

    def gc(self) -> float:
        """Calculates the fraction of the sequence that is GC.

        Returns
        -------
        A float of the fraction of the sequence that is GC.
        """

        return (self.seq.count("G") + self.seq.count("C")) / self.length
