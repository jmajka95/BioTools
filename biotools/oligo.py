from dna import DNA
from bio_enums import *
from bio_annotation import BioAnnotation
import primer3


class Oligo(DNA):
    """Class representing single-stranded oligos, typically used as primers."""
    
    def __init__(
        self, 
        seq: str, 
        name: str = "",
        type: str = BioMolecule.OLIGO,
        circular: str = BioProperty.LINEAR,
        strandedness: str = BioProperty.SINGLE_STRANDED,
    ):
        """Default constructor from superclass. 
        NOTE: Strandedness for oligos is single-stranded.
        NOTE: Oligos cannot have annotations.
        """

        super().__init__(
            seq,
            name,
            type,
            circular,
            strandedness
        )

    def info(self) -> None:
        """Prints oligo information."""
        print("Name: ", self.name)
        print("Sequence: ", self.seq)
        print("Length: ", self.length)
        print("Type: ", self.type)
        print(f"Circular: {self.circular}")
        print(f"Strandedness: {self.strandedness}")
        print(f"Tm: {primer3.calc_tm(self.seq):.2f}")

    #TODO: Generate an NEB-type melting temp calc
    def _calc_tm(self):
        raise NotImplementedError
    