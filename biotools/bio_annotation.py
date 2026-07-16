from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from biotools.dna import DNA

from biotools.bio_enums import BioOrientation
from biotools.bio_exceptions import InvalidInstantiationException
from biotools.bio_pool import BioPool

class BioAnnotation():
    """Class for representing annotations on a Biomolecule sequence. BioAnnotations
    are created by passing in a string comprised of one or more of the characters 
    defined in NTs.

    Parameters
    ----------
    span: tuple[int, int]
        A tuple of integers corresponding to the span of the 0-indexed
        sequence of the BioMolecule on which the annotation resides
    name: str
        An identifying name of the BioAnnotation
    orientation: BioOrientation
        The orientation, either forward (top strand) or reverse (bottom strand)
    """

    def __init__(
        self,
        span: tuple[int, int],
        name: str,
        orientation: BioOrientation
    ):
        """Default constructor."""
        if not isinstance(span, tuple):
            raise InvalidInstantiationException("Must provide a tuple of ints!")
        elif not isinstance(span[0], int) or not isinstance(span[1], int):
            raise InvalidInstantiationException("Must provide ints in the tuple!")
        elif span[0] < 0 or span[1] < 0:
            raise InvalidInstantiationException("Can't provide negative span indices!")
        elif span[0] == span[1]:
            raise InvalidInstantiationException("Spans can't be the same value!")
        elif len(span) != 2:
            raise InvalidInstantiationException("Spans must be exactly 2 values!")
        self.span = span # NOTE: Span is 0-indexed
        self.name = name
        if orientation not in BioOrientation:
            raise InvalidInstantiationException("Must provide a valid BioOrientation!")
        self.orientation = orientation

    def __repr__(self):
        """Self is represented as name, span, and length."""
        if self.orientation == BioOrientation.FORWARD:
            return f">>> Annotation | {self.name} | [Span: {self.span}] | [Length: {self.span[1] - self.span[0]}] >>>"
        else:
            return f"<<< Annotation | {self.name} | [Span: {self.span}] | [Length: {self.span[1] - self.span[0]}] <<<"

    def __eq__(self, other: BioAnnotation) -> bool:
        """Two BioAnnotations are equal if every field is equal."""
        return (self.span, self.name, self.orientation) == \
               (other.span, other.name, other.orientation)
    
    def __leq__(self, other: BioAnnotation) -> bool:
        return self.span[0] <= other.span[0]
    
    def __lt__(self, other: BioAnnotation) -> bool:
        return self.span[0] < other.span[0]
    
    def __hash__(self) -> int:
        """Hashes are calculated by hashing the annotation's name, span, and orientation."""
        return hash(self.name) + hash(self.span) + hash(self.orientation)
    
    def copy(self) -> BioAnnotation:
        """Returns a copy of self"""
        return BioAnnotation(self.span, self.name, self.orientation)

class Block():
    """Annotation class representing a Pool of sequences..
    
    Parameters
    ----------
    span: tuple[int, int]
        THE
    name: str,
        THE
    orientation: BioOrientation:
        THE
    pool: BioPool
        THE
    """
    
    def __init__(
        self,
        span: tuple[int, int],
        name: str,
        orientation: BioOrientation,
        pool: BioPool
    ):
        """Default constructor."""
        self.span = span
        self.name = name 
        self.orientation = orientation
        self.pool = pool

    def __repr__(self):
        if self.orientation == BioOrientation.FORWARD:
            return f">>> Block | {self.name} | [Span: {self.span}] | Length: {self.span[1] - self.span[0]} | [{self.pool.length} Sequence(s)] >>>"
        else:
            return f"<<< Block | {self.name} | [Span: {self.span}] | Length: {self.span[1] - self.span[0]} | [{self.pool.length} Sequence(s)] <<<"
        
    def __eq__(self, other: BioAnnotation) -> bool:
        """Two BioAnnotations are equal if every field is equal."""
        return (self.span, self.name, self.orientation, self.pool) == \
               (other.span, other.name, other.orientation, other.pool)
    
    def __leq__(self, other: BioAnnotation) -> bool:
        return self.span[0] <= other.span[0]
    
    def __lt__(self, other: BioAnnotation) -> bool:
        return self.span[0] < other.span[0]
    
    def __hash__(self) -> int:
        """Hashes are calculated by hashing the annotation's name, span, and orientation."""
        seq_hash: int = 0
        for seq in self.pool.seqs:
            seq_hash += hash(seq)
        return seq_hash + hash(self.name) + hash(self.span) + hash(self.orientation)
    
    def __len__(self) -> int:
        return self.span[1] - self.span[0]
    
    def print_pool(self, print_seq: bool = False) -> None:
        """Prints the sequences in self.pool
        
        Parameters
        ----------
        print_seq: bool
            Whether or not to print the sequence instead of the representation of the pool
            
        Returns
        -------
        None
        """
        for i, seq in enumerate(self.pool.seqs):
            if print_seq:
                print(f"{i + 1}:" + " " * (8 - len(str(i + 1))) + f"{seq.seq}")
            else:
                print(seq)

    def copy(self) -> BioAnnotation:
        """Returns a copy of self"""
        return Block(self.span, self.name, self.orientation, self.pool)

######################### FUNCTIONS #########################

def reverse_annotations(
    annotations: list[BioAnnotation], length: int, rev_comp: bool = False
) -> list[BioAnnotation]:
    """Reverses the provided annotations based on the length of the input sequence.
    
    Parameters
    ----------
    annotations: list[BioAnnotation]
    length: int
    rev_comp: bool"""
    reversed_annotations: list[BioAnnotation] = []
    for a in annotations:
        new_span = (length - a.span[1], length - a.span[0])
        if rev_comp:
            orient = BioOrientation.FORWARD if a.orientation == BioOrientation.REVERSE else BioOrientation.REVERSE
        else:
            orient = a.orientation

        if isinstance(a, BioAnnotation):
            reversed_annotations.append(BioAnnotation(new_span, a.name, orient))
        elif isinstance(a, Block):
            reversed_annotations.append(Block(new_span, a.name, orient, a.pool)) # TODO: Reverse the pool as well? Yes!

    return reversed_annotations
