"""File defining utility functions to be used for simulation.
"""

from __future__ import annotations

from typing import TYPE_CHECKING, Iterator, Any
if TYPE_CHECKING:
    from biotools.dna import DNA
    from biotools.rna import RNA
    from biotools.protein import Protein
    from biotools.oligo import Oligo
   
from biotools.bio_alphabet import (
    NT_COMPLEMENTS, RNA_COMPLEMENTS, NTs, AAs,
    RNA_NTs, DNA_TO_RNA, RNA_CODON_AAs, AA_CODONS,
    RE_ENZYMES_REGEX, RE_ENZYMES_OFFSETS
)
from biotools.bio_enums import BioMolecule
from biotools.bio_exceptions import InvalidInstantiationException, ReactionError
from biotools.bio_annotation import BioAnnotation, BioOrientation

import numpy as np
import random
from random import choice, sample, randint
import re
from re import Match
import requests # type: ignore
import json

def rev_comp(seq: str, is_dna: bool = True) -> str:
    """Generates the reverse complement of an input."""
    seq = seq.upper()
    if is_dna:
        return "".join([NT_COMPLEMENTS[char] for char in seq[::-1]])
    return "".join([RNA_COMPLEMENTS[char] for char in seq[::-1]])

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
    return all(char in AAs for char in seq) # Default to amino acids

def validate_annotations(seq: DNA, annotations: list[BioAnnotation]) -> bool:
    """Validates annotation spans for annotations being added to a sequence."""
    for annot in annotations:
        if annot.span[0] < 0 or annot.span[1] < 0:
            raise InvalidInstantiationException(f"Invalid annotation provided! One span is less than zero!")
        elif annot.span[0] > len(seq) or annot.span[1] > seq.length:
            print(annot.span) # DEBUG
            raise InvalidInstantiationException(f"Invalid annotation provided! One span is higher than the sequence length!")
    return True

def transcribe(dna: DNA) -> RNA:
    """Transcribes DNA into RNA."""
    from biotools.rna import RNA
    return RNA(seq="".join([DNA_TO_RNA[char] for char in dna.seq]))

def translate(rna: RNA) -> Protein:
    """Translates RNA into protein."""
    from biotools.protein import Protein

    if rna.length % 3 != 0:
        raise Exception("Cannot translate RNA sequence that contains a partial codon.")
    return Protein(seq="".join([RNA_CODON_AAs[rna.seq[i:i+3]] for i in range(0, rna.length, 3)]))

def reverse_translate(aa_seq: Protein, rand: bool = False) -> DNA:
    """Reverse translates an amino acid sequence by random codon choice."""
    from biotools.dna import DNA
    if not rand:
        random.seed(1) # For consistent behavior
    return DNA(seq="".join([choice(AA_CODONS[aa_seq.seq[i]]) for i in range(0, len(aa_seq))]))

def contains_re_site(seq: DNA, site: str) -> bool:
    """Returns True if the DNA sequence contains an re site, else False."""
    if seq.is_circular():
        return (site in seq.top_strand+seq.top_strand[: len(site) - 1]) \
            or (site in seq.bottom_strand+seq.bottom_strand[: len(site) - 1])
    return (site in seq.top_strand) or (site in seq.bottom_strand)

def generate_random_sequence(
    length: int | tuple[int, int],
    seq_type: BioMolecule,
    n_seqs: int = 1,
    enzyme_blacklist: list[str] | None = None,
    kwargs: dict[str, Any] = {}
) -> DNA | list[DNA] | Protein | RNA | list[RNA] | list[Protein]:
    """Generates a random sequence of either amino acids or nucleotides of a specified length.
    
    Parameters
    ----------
    length: int | tuple[int, int]
        If an int, it is the length of the sequence to generate.
        If an int tuple, the range of values from which to generate sequences randomly
    seq_type: BioMolecule
        A BioMolecule to generate. Must be one of DNA, RNA, or Protein
    n_seqs: int (Optional)
        An integer of the number of sequences to generate 
    enzyme_blacklist: list[str] (Optional)
        A list of enzymes to avoid
    kwargs: dict[Any] (Optional)
        A dictionary of additional arguments for constructing DNA, RNA, or Proteins

    Returns
    -------
    DNA, a list of DNA, RNA, a list of RNA, Protein, or a list of Proteins
    """
    from biotools.dna import DNA
    from biotools.rna import RNA
    from biotools.protein import Protein

    # TODO: Make recursive to add multiprocessing?
    nt_keys = list(DNA_TO_RNA.keys())
    rna_keys = list(DNA_TO_RNA.values())
    aa_keys = list(AA_CODONS.keys())

    seq_list = []
    for n in range(n_seqs):
        seq: str = ""
        seq_length: int = 0
        if isinstance(length, tuple):
            seq_length = randint(length[0], length[1])
        else:
            seq_length = length
        if enzyme_blacklist: # Check for restriction enzymes
            valid = False
            while not valid:
                seq: str = ""
                valid = True
                for i in range(seq_length):
                    if seq_type == BioMolecule.DNA:
                        p = choice(nt_keys)
                    elif seq_type == BioMolecule.RNA:
                        p = choice(rna_keys)
                    elif seq_type == BioMolecule.PROTEIN:
                        p = choice(aa_keys)
                        while (p == "*"):
                            p = choice(aa_keys)
                    seq += p
                for e in enzyme_blacklist:
                    m = re.search(RE_ENZYMES_REGEX[e], seq)
                    if m:
                        valid = False
                        break
                    m = re.search(RE_ENZYMES_REGEX[e], rev_comp(seq))
                    if m:
                        valid = False
                        break
        else: # No restriction enzyme validation
            for i in range(seq_length):
                    if seq_type == BioMolecule.DNA:
                        p = choice(nt_keys)
                    elif seq_type == BioMolecule.RNA:
                        p = choice(rna_keys)
                    elif seq_type == BioMolecule.PROTEIN:
                        p = choice(aa_keys)
                        while (p == "*"):
                            p = choice(aa_keys)
                    seq += p
        if seq_type == BioMolecule.DNA:
            seq_list.append(DNA(seq=seq, **kwargs))
        elif seq_type == BioMolecule.RNA:
            seq_list.append(RNA(seq=seq, **kwargs))
        elif seq_type == BioMolecule.PROTEIN:
            seq_list.append(Protein(seq=seq, **kwargs))
    return seq_list[0] if len(seq_list) == 1 else seq_list

def generate_random_mutations(seq: str, mut_range: tuple[int, int], num_muts: int) -> str:
    """Generates random mutants within the specified indices of the provided sequence."""

    if (mut_range[1] <= mut_range[0]):
        raise RuntimeError("Provided range must be in form (<initial index>, <final index>) with a difference of at least 1.")
    if (num_muts > (mut_range[1] - mut_range[0]) + 1) or (num_muts <= 0):
        raise RuntimeError("Number of mutants must be between 1 and the number of possible mutation locations.")
    if len(seq) <= 0:
        raise RuntimeError("Must submit a valid sequence.")

    indices = sample(range(mut_range[0], mut_range[1]+1), num_muts)
    keys = list(AA_CODONS.keys())

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

def find_re_sites(seq: DNA, *enzymes: str) -> dict:
    """Finds all instances of the provided enzymes on the forward and reverse strand
    of a DNA sequence. Returns a dictionary mapping the site to a list of tuples of their
    spans for both the forward and the reverse strand.
    """
    try:
        enzyme_seqs = [RE_ENZYMES_REGEX[enzyme] for enzyme in enzymes]
    except KeyError:
        raise ReactionError(f"One of provided enzymes not found. Must use one of: {sorted(list(RE_ENZYMES_REGEX.keys()))}")

    site_span_dict = {}
    for i, enzyme in enumerate(enzymes):
        site_span_dict[enzyme] = {}
        site_span_dict[enzyme][BioOrientation.FORWARD] = []
        site_span_dict[enzyme][BioOrientation.REVERSE] = []

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
            site_span_dict[enzyme][BioOrientation.FORWARD].append(f_span)

        r_matches = re.finditer(enzyme_seqs[i], rev_sequence)
        for match in r_matches:
            r_span = match.span()
            if r_span[1] > seq.length:
                r_span = (r_span[0], r_span[1] - seq.length)
            f_span = (seq.length - r_span[1], seq.length - r_span[0])
            if not f_span in site_span_dict[enzyme][BioOrientation.FORWARD]: # Don't want duplicate spans if there's a palindromic recognition site
                site_span_dict[enzyme][BioOrientation.REVERSE].append(r_span)

    return site_span_dict

def validate_sites(
    site_dict: dict, seq_len: int, circular: bool, padding: int = 6
) -> tuple[int, int, int, tuple[int, int], int]:
    """Helper function for digest().
    Validates whether or not a cut site(s) exists and can be used.
    Returns a dictionary mapping the cut site tuples with their respective enzyme."""
    spans: list[tuple[int, int]] = []
    tuple_list: list[tuple[tuple[int,int,BioOrientation], str]] = []
    enzymes = list(site_dict.keys())

    for i, val in enumerate(enzymes):
        # A site must exist for each enzyme for either Forward or Reverse
        if not site_dict[val][BioOrientation.FORWARD] and not site_dict[val][BioOrientation.REVERSE]:
            raise ReactionError(f"No sites detected for {val}!")

        for orientation in site_dict[val].keys():
            for span in site_dict[val][orientation]:
                if not circular:
                    if span:
                        if orientation == BioOrientation.FORWARD:
                            if span[0] < padding:
                                raise ReactionError(f"Must have at least {padding} bases for digestion to work if input is linear!")
                        else: # Checking reverse orientation
                            if span[1] >= seq_len - padding:
                                raise ReactionError(f"Must have at least {padding} bases for digestion to work if input is linear!")
                if orientation == BioOrientation.REVERSE:
                    span = (seq_len - span[1], seq_len - span[0]) # Convert to forward orientation
                spans.append(span)

                # Generate list for ( (int, int, int, BioOrientation) , (int, int, str, BioOrientation) )
                if isinstance(RE_ENZYMES_OFFSETS[val], list):
                    offsets = RE_ENZYMES_OFFSETS[val] if orientation == BioOrientation.FORWARD else RE_ENZYMES_OFFSETS[val][::-1]
                    for tup in offsets:
                        tuple_list.append((tup, span, val, orientation))
                else:
                    tuple_list.append((RE_ENZYMES_OFFSETS[val], span, val, orientation))
                                       
    if not spans: 
        raise ReactionError("No valid spans found!")
    
    seen: set[tuple[tuple[int,int,BioOrientation], str]] = set()
    tuple_list = [tup for tup in tuple_list if not (tup in seen or seen.add(tup))] # Clean up repeats

    cut_spans: list[tuple[int, int]] = []
    cuts_tuple_list: tuple[int, int, int, tuple[int, int], int] = []
    # Validate span overlaps
    for tup_1 in tuple_list:
        for tup_2 in tuple_list:
            if tup_1 != tup_2:
                # Get recognition site tuples
                site_1 = tup_1[1]
                site_2 = tup_2[1]
                if tup_1[-1] == BioOrientation.FORWARD:
                    e1_cut_span = (site_1[0] + tup_1[0][0], tup_1[0][0] + tup_1[0][0] + tup_1[0][1])
                else: # Must reverse from span[1] because it's 5' bottom strand
                    e1_cut_span = (site_1[1] - tup_1[0][0], site_1[1] - tup_1[0][0] - tup_1[0][1])
                cuts_tuple_list.append((tup_1[0], tup_1[1], tup_1[2], e1_cut_span, tup_1[3]))
                cut_spans.append(e1_cut_span)
    
                if tup_2[-1] == BioOrientation.FORWARD:
                    e2_cut_span = (site_2[0] + tup_2[0][0], site_2[0] + tup_2[0][0] + tup_2[0][1])
                else:
                    e2_cut_span = (site_2[1] - tup_2[0][0], site_2[1] - tup_2[0][0] - tup_2[0][1])

                # Check for cutting the re site
                if (e1_cut_span[0] > site_2[0] and e1_cut_span[0] < site_2[1]) or (e1_cut_span[1] > site_2[0] and e1_cut_span[1] < site_2[1]) \
                or (e1_cut_span[0] < site_2[0] and e1_cut_span[1] > site_2[1]):
                    raise ReactionError("Enzyme cuts re site!")
                # Check for cutting the cut site
                if (e1_cut_span[0] > e2_cut_span[0] and e1_cut_span[0] < e2_cut_span[1]) or (e1_cut_span[1] > e2_cut_span[0] and e1_cut_span[1] < e2_cut_span[1]):
                    raise ReactionError("Enzyme cuts cut site!")
                
    # Sort based on cut tuple
    sorted_tuple_list: list[tuple[int, int, int, tuple[int, int], int]] = []
    for s in sorted(cut_spans):
        for tup in cuts_tuple_list:
            if tup[3] == s and tup not in sorted_tuple_list:
                sorted_tuple_list.append(tup)

    return sorted_tuple_list

def check_homology(p: DNA, q: DNA, min_len: int, max_len: int) -> int | None:
    """Helper function for gibson().
    Checks if there is any homology from min_len to max_len. Compares p's 5' end with q's 3' end.
    Returns the length of homology that exists between p and q.
    """
    max_len = min(max_len, len(p), len(q))
    for i in range(min_len, max_len + 1):
        if p[:i].top_strand == q[-i:].top_strand:
            if len(p) == i or len(q) == i:  # Edge case where the entire sequence matches if short enough
                return None
            return i
    return None

def get_annealing_temp(
    primer_1: list[Oligo], 
    primer_2: list[Oligo], 
    conc: float = 0.5, 
    prod_code: str = "q5hs-1"
) -> int | list[int]:
    """Generates an annealing temperature of the provided primers. This uses NEB's Tm API.
    More information can be found at https://tmapi.neb.com/
    NOTE: prod_code has been defaulted to that corresponding to Q5 2X Hot Start."""
    
    url = 'https://tmapi.neb.com/tm/batch'
    primer_pairs = [(p1.seq, p2.seq) for p1, p2 in zip(primer_1, primer_2)]
    input = {
        "seqpairs": primer_pairs,
        "conc": conc,
        "prodcode": prod_code
    }
    headers = {"content-type" : "application/json"}
    res = requests.post(url, data=json.dumps(input), headers=headers)
    r = json.loads(res.content)

    annealing_temps: list[int] = []
    if r['success']:
        for row in r['data']:
            annealing_temps.append(row['ta'])
    else:
        print(f"Request failure. Error code: {r['error'][0]}")

    return annealing_temps[0] if len(annealing_temps) == 1 else annealing_temps

# TODO: Alignment algorithm? Add to bio_io.py instead?
def align_sequences(seq: str, ref: str) -> str:
    """Aligns two sequences together using the X format."""
    raise NotImplementedError

# TODO: Melting temp (NEB calculator, but also allow for Taq annealing etc.)
