from dna import DNA
from bio_enums import *
from bio_annotation import BioAnnotation
import json
import requests # type: ignore

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
        # TODO: Make this its own class. Starting to make less sense having it be subclass of DNA

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
        print(f"Tm: {self.calc_tm()}")

    def calc_tm(self, prod_code: str = "q5hs-1"):
        """Calculates a melting temperature for the primer"""
        res = requests.get(f"https://tmapi.neb.com/tm/{prod_code}/0.5/{self.seq}")

        r = json.loads(res.content)
        if not r['success']:
            raise Exception(f"Failed to retrieve Tm. Error code: {r['error'][0]}")
        return r['data']['tm1']
