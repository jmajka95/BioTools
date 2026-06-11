from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from dna import DNA

from bio_enums import BioOrientation
from bio_exceptions import *
from bio_pool import BioPool

class BioAnnotation():
    """Class for representing annotations on a Biomolecule sequence"""

    # Block annotation that contains sequences?
    # Regular annotation for simple labeling?

    def __init__(
        self,
        span: tuple[int, int],
        name: str,
        orientation: BioOrientation
    ):
        """Default constructor for BioAnnotation."""
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
        if self.orientation == BioOrientation.FORWARD:
            return f">>> Annotation | {self.name} | [Span: {self.span}] >>>" # TODO: More elaborate span length calculation??
        else:
            return f"<<< Annotation | {self.name} | [Span: {self.span}] <<<" 

    def __eq__(self, other: BioAnnotation) -> bool:
        return (self.span, self.name, self.orientation) == \
               (other.span, other.name, other.orientation)
    
    def __leq__(self, other: BioAnnotation) -> bool:
        return self.span[0] <= other.span[0]
    
    def __lt__(self, other: BioAnnotation) -> bool:
        return self.span[0] < other.span[0]

class Block(BioAnnotation):
    """Annotation class representing a Pool of sequences."""
    
    def __init__(
        self, 
        span: tuple[int, int],
        name: str,
        orientation: BioOrientation,
        pool: BioPool
    ):
        super().__init__(span, name, orientation)
        self.pool = pool

    def __repr__(self):
        if self.orientation == BioOrientation.FORWARD:
            return f">>> BLOCK | {self.name} | [Span: {self.span}] | Length: {self.span[1] - self.span[0]} | [{self.pool.length} Sequence(s)] >>>"
        else:
            return f"<<< BLOCK | {self.name} | [Span: {self.span}] | Length: {self.span[1] - self.span[0]} | [{self.pool.length} Sequence(s)] <<<"

######################### FUNCTIONS #########################

def reverse_annotations(annotations: list[BioAnnotation], length: int, rev_comp: bool = False) -> list[BioAnnotation]:
    """Reverses the provided annotations"""
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
            reversed_annotations.append(Block(new_span, a.name, orient, a.pool)) # TODO: Reverse the pool as well?

    return reversed_annotations
