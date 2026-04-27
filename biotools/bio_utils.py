from __future__ import annotations

from typing import TYPE_CHECKING, Iterator
if TYPE_CHECKING:
    from dna import DNA
    from rna import RNA
    from protein import Protein
    from oligo import Oligo
    from re import Match
    
from bio_alphabet import *
import numpy as np
from random import choice, sample
from bio_enums import *
from bio_exceptions import *
from bio_pool import BioPool
from bio_annotation import BioAnnotation
import random
import re

# TODO: Implement all of this for DNA instead of str?

# File for defining utility functions #

def rev_comp(seq: str, is_dna: bool = True) -> str:
    """Generates the reverse complement of an input."""
    seq = seq.upper()
    if is_dna:
        return "".join([nt_complements[char] for char in seq[::-1]])
    return "".join([rna_complements[char] for char in seq[::-1]])

def calc_mut_bounds(seq: str, ref: str) -> tuple[int, int, int]:
    """Function that generates mutation bounds comparing a sequence to a reference."""
    seq = seq.upper()
    ref = ref.upper()

    if (seq == ref):
        return (np.nan, np.nan, np.nan)

    # X ensures insertions check out because it adds padding that will always catch the sequence difference
    for pos, (seq1, seq2) in enumerate(zip(seq + "X", ref + "X")):
        if (seq1 != seq2):
            start_pos = pos
            break
    for pos, (seq1, seq2) in enumerate(zip(seq[start_pos:][::-1] + "X", ref[start_pos:][::-1] + "X")):
        if (seq1 != seq2):
            seq_end_pos = len(seq) - pos
            ref_end_pos = len(ref) - pos
            break

    # Sanity checks
    assert start_pos <= seq_end_pos
    assert start_pos <= ref_end_pos
    return start_pos, seq_end_pos, ref_end_pos
    
def validate_sequence(seq: str, type: BioMolecule) -> bool:
    """Helper function that validates an input sequence based on valid characters."""
    if type == BioMolecule.DNA:
        return all(char in NTs for char in seq)
    elif type == BioMolecule.RNA:
        return all(char in RNA_NTs for char in seq)
    return all(char in AAs for char in seq)

def validate_annotations(seq: DNA, annotations: list[BioAnnotation]) -> bool:
    """Validates annotation spans for annotations being added to a sequence."""
    for annot in annotations:
        if annot.span[0] < 0 or annot.span[1] < 0:
            raise InvalidInstantiationException(f"Invalid annotation provided! One span is less than zero!")
        elif annot.span[0] > len(seq) or annot.span[1] > seq.length:
            raise InvalidInstantiationException(f"Invalid annotation provided! One span is higher than the sequence length!")
    return True

def transcribe(dna: DNA) -> RNA:
    """Transcribes DNA into RNA."""
    from rna import RNA
    return RNA(seq="".join([dna_to_rna[char] for char in dna.seq]))

def translate(rna: RNA) -> Protein:
    """Translates RNA into protein."""
    from protein import Protein

    if rna.length % 3 != 0:
        raise Exception("Cannot translate RNA sequence that contains a partial codon.")
    return Protein(seq="".join([rna_codon_AAs[rna.seq[i:i+3]] for i in range(0, rna.length, 3)]))

def reverse_translate(aa_seq: Protein) -> str:
    """Reverse translates an amino acid sequence by random codon choice."""
    from dna import DNA
    random.seed(1) # For consistent behavior
    return DNA(seq="".join([choice(aa_codons[aa_seq.seq[i]]) for i in range(0, len(aa_seq))]))

def contains_re_site(seq: DNA, site: str) -> bool:
    """Returns True if the DNA sequence contains an re site, else False."""
    return (site in seq.top_strand) or (site in seq.bottom_strand)

def generate_random_sequence(length: int, nt: bool = False) -> str:
    """Generates a random sequence of either amino acids or nucleotides of a specified length."""
    seq: str = ""
    keys = list(aa_codons.keys())

    for i in range(length):
        aa = choice(keys)
        while (aa == "*"):
            aa = choice(keys)
        seq += aa

    if nt:
        seq = reverse_translate(seq)

    return seq

def generate_random_mutations(seq: str, mut_range: tuple[int, int], num_muts: int) -> str:
    """Generates random mutants within the specified indices of the provided sequence."""

    if (mut_range[1] <= mut_range[0]):
        raise RuntimeError("Provided range must be in form (<initial index>, <final index>) with a difference of at least 1.")
    if (num_muts > (mut_range[1] - mut_range[0]) + 1) or (num_muts <= 0):
        raise RuntimeError("Number of mutants must be between 1 and the number of possible mutation locations.")
    if len(seq) <= 0:
        raise RuntimeError("Must submit a valid sequence.")

    indices = sample(range(mut_range[0], mut_range[1]+1), num_muts)
    keys = list(aa_codons.keys())

    for idx in indices:
        aa = choice(keys)
        while (aa == "*" or aa == seq[idx]):
            aa = choice(keys)
        seq = seq[:idx] + aa + seq[idx+1:]

    return seq

def find_binding_site(
    primer: Oligo,
    template: DNA,
    min_binding_length: int = 15,
    single: bool = False
) -> tuple[int, int] | Iterator[Match[str]] | None:
    """Finds a binding site(s) for a primer."""
    match = re.finditer(primer.seq[-min_binding_length:], template.seq)
    
    if single:
        while len(list(match)) != 1:
            if (min_binding_length := min_binding_length + 1) > primer.length:
                return None
            match = re.finditer(primer.seq[-min_binding_length:], template.seq)

        match = re.finditer(primer.seq[-min_binding_length:], template.seq)
        return next(match).span()
    else:
        while len(match) == 0:
            if (min_binding_length := min_binding_length + 1) > primer.length:
                return None
            match = re.finditer(primer.seq[-min_binding_length:], template.seq)

    match = re.finditer(primer.seq[-min_binding_length:], template.seq)
    return match

# TODO: Do we need Enzyme class for anything specific?
def find_re_sites(seq: DNA, *enzymes: str) -> dict:
    """Finds all instances of the provided enzymes on the forward and reverse strand
    of a DNA sequence. Returns a dictionary mapping the site to a list of tuples of their
    spans for both forward and reverse strand.
    """
    enzyme_seqs = [re_enzymes[enzyme] for enzyme in enzymes]

    site_span_dict = {}
    for i, enzyme in enumerate(enzymes):
        site_span_dict[enzyme] = {}
        site_span_dict[enzyme]["Forward"] = []
        site_span_dict[enzyme]["Reverse"] = []

        if seq.is_circular(): # Pad sequence up to (enzyme length - 1) because wraparound is valid
            sequence = seq.top_strand + seq.top_strand[:len(enzyme) - 1]
            rev_sequence = seq.bottom_strand + seq.bottom_strand[:len(enzyme) - 1]
        else:
            sequence = seq.top_strand
            rev_sequence = seq.bottom_strand

        f_matches = re.finditer(enzyme_seqs[i], sequence)
        for match in f_matches:
            f_span = match.span()
            if f_span[1] >= seq.length:
                f_span = (f_span[0], f_span[1] - seq.length)
            site_span_dict[enzyme]["Forward"].append(f_span)

        r_matches = re.finditer(enzyme_seqs[i], rev_sequence)
        for match in r_matches:
            r_span = match.span()
            if r_span[1] > seq.length:
                r_span = (r_span[0], r_span[1] - seq.length)
            site_span_dict[enzyme]["Reverse"].append(r_span)

    return site_span_dict

def validate_sites(site_dict: dict, seq_len: int, circular: bool, padding: int = 6) -> bool:
    """Helper function for digest().
    Validates whether or not a cut site exists and can be used."""
    # TODO: Also have this validate that multiple cut sites aren't too close
    valid: bool = False
    spans: list[tuple[int, int]] = []
    for val in site_dict.values():
        for orientation in val.keys():
            for lst in val[orientation]:
                if not circular:
                    if lst:
                        valid = True
                        for span in lst:
                            if span[0] < padding: # Must have at least padding bases for digestion to work if linear
                                return False
                else: # Padding not an issue if circular
                    if lst:
                        valid = True
                if orientation == "Reverse":
                    lst = (seq_len - lst[0], seq_len - lst[1]) # Convert to forward orientation
                    spans.append(lst)
                else:
                    spans.append(lst)
    if not spans: return False

    # Check span overlaps?
    for s1 in spans:
        for s2 in spans:
            if s1 != s2:
                pass # re_enzymes_offsets

    if valid: return True
    return False

# TODO: Generate random pool of DNA seqs?

# TODO: Alignment algorithm? Add to bio_io.py instead?
def align_sequences(seq: str, ref: str) -> str:
    """Aligns two sequences together using the X format."""
    raise NotImplementedError

# TODO: Melting temp (NEB calculator, but also allow for Taq annealing etc.)
