from __future__ import annotations

from typing import TYPE_CHECKING
if TYPE_CHECKING:
    from dna import DNA

class BioPool():
    """Class representing a pool of sequences."""
    
    def __init__(
        self,
        seqs: list[DNA] | None = None
    ):
        self.seqs = [] if seqs is None else seqs

    @property
    def length(self):
        return len(self.seqs)

    # TODO: What interesting things can we do with this?