from biotools.dna import DNA
from biotools.bio_enums import BioMolecule, BioProperty
from biotools.bio_annotation import BioAnnotation, Block 
from biotools.bio_exceptions import InvalidAnnotationException

import json
import requests # type: ignore
import bisect
from itertools import product

class Oligo(DNA):
    """Class representing single-stranded oligos, typically used as primers.
    Oligos may also be annealed to form double-stranded DNA.
    
    Parameters
    ----------
    seq: str
        A string sequence of DNA. Must contain only valid nucleotide or degenerate base
        characters as specified in valid_chars
    type: BioMolecule
        An identifier Enum of the type DNA
    circular: BioProperty
        A boolean of whether or not the sequence is circular (a plasmid)
    strandedness: BioProperty
        A property stating whether or not the DNA is single- or double-stranded
    annotations: list[BioAnnotation | Block]
        A list of annotations present on the DNA
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
        """Default constructor. 
        NOTE: Strandedness for oligos is single-stranded.
        NOTE: Oligos cannot have annotations.
        """
        # TODO: Make this its own class. Starting to make less sense having it be subclass of DNA

        super().__init__(
            seq,
            name,
            type,
            circular,
            strandedness,
            annotations
        )

    @property
    def length(self) -> int:
        """The length of the DNA sequence."""
        return len(self.seq)

    def has_pool(self) -> bool:
        for annot in self.annotations:
            if isinstance(annot, Block):
                return True
        return False

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

    def add_annotations(self, *annotations: BioAnnotation | Block) -> None:
        """Adds an annotation or annotations to the DNA sequence.
        
        Parameters
        ----------
        annotation: BioAnnotation | Block
            One or more annotations to add to self
            
        Returns
        -------
        None
        """

        # Spans can't be same
        for annot in annotations:
            if annot.span[0] == annot.span[1]:
                raise InvalidAnnotationException("Must not have same span indices!")

        # Validate tuple based on circularity
        if not self.is_circular():
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

    def get_pools(self) -> list[DNA] | None:
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
        groups: list[list[DNA]] = []

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
        seqs: list[DNA] = []
        block_combos = product(*groups)
        for combo in block_combos:
            new_seq: str = ""
            for i in range(len(combo)):
                new_seq += self_seqs[i]
                new_seq += combo[i].seq
            new_seq += self_seqs[len(combo)]
            seqs.append(
                DNA(
                    seq=new_seq,
                    circular=self.circular,
                    strandedness=self.strandedness,
                    offsets=self.offsets,
                    parent=self
                )
            )

        seqs.append(self.copy()) # Add our own seq because it's not captured in the Block combinations
        return seqs
    
