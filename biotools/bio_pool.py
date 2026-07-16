from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from biotools.dna import DNA
    from biotools.rna import RNA
    from biotools.protein import Protein

class BioPool():
    """Class representing a pool of sequences. Pools are able to represent a
    library/group of sequences. These groups can have """
    
    def __init__(
        self,
        seqs: list[DNA] | list[RNA] | list[Protein] | None = None
    ):
        """Default constructor."""
        self.seqs = [] if seqs is None else seqs

    def __repr__(self):
        return f"BioPool | [{len(self.seqs)} Sequence(s)]"

    @property
    def length(self):
        return len(self.seqs)

    # TODO: What interesting things can we do with this?