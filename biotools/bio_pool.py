from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from biotools.dna import DNA
    from biotools.rna import RNA
    from biotools.protein import Protein
    from biotools.oligo import Oligo

class BioPool():
    """Class representing a pool of sequences. Pools are able to represent a
    library/group of sequences. These groups can have varying lengths
    but should be the same type of BioMolecule.
    
    Parameters
    ----------
    seqs: list[DNA] | list[RNA] | list[Protein] | list[Oligo] | None (Default: None)
        The sequences of the pool
    """
    
    def __init__(
        self,
        seqs: list[DNA] | list[RNA] | list[Protein] | list[Oligo] | None = None
    ):
        """Default constructor."""
        self.seqs = [] if seqs is None else seqs

    def __repr__(self):
        return f"{{ BioPool | [{len(self.seqs)} Sequence(s)] }}"

    @property
    def length(self):
        return len(self.seqs)

    # TODO: What interesting things can we do with this?