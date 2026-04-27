# Functions supporting reactions using DNA sequences
from dna import DNA
from oligo import Oligo
from bio_exceptions import ReactionError
from bio_enums import *
from bio_utils import find_binding_site, rev_comp, find_re_sites, validate_sites
from bio_alphabet import re_enzymes, re_enzymes_offsets
import re

# TODO: Would be cool to simulate reactions in graphs similar to what Nirmit did.
# Probably in a different class that can set up reactions in
# Probably can just make a compiler class that utilizes some functionality from the CombinatorialCompiler
# Should put assign_codons() in utils

def amplify(
    f_primer: list[Oligo], 
    r_primer: list[Oligo], 
    template: DNA, 
    min_binding_length: int = 15
) -> DNA | list[DNA]:
    """Amplifies DNA using provided primers that must bind at least min_binding_length bases."""
    if any(len(f) < min_binding_length for f in f_primer) or any(len(r) < min_binding_length for r in r_primer):
        raise ReactionError(f"Minimum length of primer must be at least {min_binding_length} bp!")
    
    f_binding_spans = []
    for f in f_primer:
        site = find_binding_site(primer=f, template=template, single=True)
        if site:
            f_binding_spans.append(site)
        else:
            raise ReactionError(f"No single binding site found for f primer {f.name}!")

    r_binding_spans = []
    for r in r_primer:
        site = find_binding_site(primer=r, template=template.rev_comp(), single=True)
        if site:
            # Convert to forward span indices
            new_r_span = (template.length - site[1], template.length - site[0])
            r_binding_spans.append(new_r_span)
        else:
            raise ReactionError(f"No single binding site found for f primer {r.name}!")
        
    # Check if all f priming sites are before r priming sites and not less than min_binding_length away
    if not template.is_circular():
        if not all(f_span[1] - r_span[0] >= 15 for f_span in f_binding_spans for r_span in r_binding_spans):
            raise ReactionError("Forward primers must be upstream of reverse primers and must not overlap!")
    else:
        if not all(f_span[0] < r_span[1] for f_span in f_binding_spans for r_span in r_binding_spans): # TODO: Can span[0] and span[1] be the same?
            raise ReactionError("Forward primers must be not overlap with reverse primers!")
        
    final_products: list[DNA] = []
    # Generate PCR product(s)
    for i, f_span in enumerate(f_binding_spans):
        for j, r_span in enumerate(r_binding_spans):
            product = template[f_span[0]:r_span[1]]
            if f_span[1] - f_span[0] < f_primer[i].length: # Need to take difference in primer
                product = DNA(seq=f_primer[i][:f_primer[i].length - (f_span[1] - f_span[0])].seq) + product
            if r_span[1] - r_span[0] < r_primer[j].length: # Need to take difference in primer
                product = product + DNA(seq=rev_comp(r_primer[j][:r_primer[j].length - (r_span[1] - r_span[0])].seq))
            product.name = template.name+f"_PCR_Product_{i+j+1}"
            final_products.append(product)

    if len(final_products) == 1:
        return final_products[0]
    return final_products

def kld(input: DNA) -> DNA:
    """Performs KLD on a linear DNA fragment.
    For more information, see https://www.neb.com/en-ca/protocols/kld-enzyme-mix-reaction-protocol-m0554
    """
    # TODO: Length check necessary?
    if input.is_circular():
        raise ReactionError("Cannot perform KLD on a circular input!")
    return DNA(input.seq, input.type, BioProperty.CIRCULAR, input.strandedness, input.annotations)

def digest(input: DNA, enzymes: list[str]) -> DNA | list[DNA]:
    """Digests the input with the provided enzyme and returns all products formed.
    NOTE: Enzymes must be valid keys of bio_alphabet.re_enzymes"""
    sites: dict = find_re_sites(input, *enzymes)

    circular = True if input.circular == BioProperty.CIRCULAR else False
    if not validate_sites(sites, len(input) - 1, circular): # NOTE: Must provide (len - 1) for valid re-indexing
        raise ReactionError("Cannot find valid binding site(s) with enzyme and input!")

    products: list[DNA] = []
    for enzyme in enzymes:
        for span in sites[enzyme]["Forward"]:
            if input.is_circular(): # Reindex one product
                product = input.reindex(span[0] + re_enzymes_offsets[enzyme][0], False)
                product.seq += product.seq[:re_enzymes_offsets[enzyme][1]]
                product.circular = BioProperty.LINEAR
                product.strandedness = BioProperty.CUT
                product.offsets = ((0, re_enzymes_offsets[enzyme][1]), (re_enzymes_offsets[enzyme][1], 0))
                products.append(product)
            else: # Make two products
                product_1 = input[:span[0] + re_enzymes_offsets[enzyme][0]]
                product_1.seq += input.seq[span[0]+re_enzymes_offsets[enzyme][0]:span[0]+re_enzymes_offsets[enzyme][0]+re_enzymes_offsets[enzyme][1]]
                product_1.strandedness = BioProperty.CUT
                product_1.offsets = ((0, re_enzymes_offsets[enzyme][1]), (0, 0))
                product_2 = input[span[0] + re_enzymes_offsets[enzyme][0]:]
                product_2.strandedness = BioProperty.CUT
                product_2.offsets = ((0, 0), (re_enzymes_offsets[enzyme][1], 0))
                products.extend([product_1, product_2])

        for span in sites[enzyme]["Reverse"]:
            if input.is_circular(): # Reindex one product
                product = input.reindex(span[0]+re_enzymes_offsets[enzyme][0], False)
                product.seq += product.seq[:re_enzymes_offsets[enzyme][1]]
                product.circular = BioProperty.LINEAR
                product.strandedness = BioProperty.CUT
                product.offsets = ((0, re_enzymes_offsets[enzyme][1]), (re_enzymes_offsets[enzyme][1], 0))
                products.append(product)
            else: # Make two products
                product_1 = input[:span[0]+re_enzymes_offsets[enzyme][0]]
                product_1.seq += input.seq[span[0]+re_enzymes_offsets[enzyme][0]:span[0]+re_enzymes_offsets[enzyme][0]+re_enzymes_offsets[enzyme][1]]
                product_1.strandedness = BioProperty.CUT
                product_1.offsets = ((0, re_enzymes_offsets[enzyme][1]), (0, 0))
                product_2 = input[span[0]+re_enzymes_offsets[enzyme][0]:]
                product_2.strandedness = BioProperty.CUT
                product_2.offsets = ((0, 0), (re_enzymes_offsets[enzyme][1], 0))
                products.extend([product_1, product_2])

    if len(products) == 1:
        return products[0]
    return products

def ligate(*inputs: DNA) -> DNA:
    """Ligates together the provided inputs."""
    for input in inputs:
        if input.is_circular():
            raise ReactionError("Cannot ligate circular input!")
    raise NotImplementedError

def gel_extract(*inputs: DNA, extraction_len: int) -> DNA | list[DNA]:
    raise NotImplementedError

def gibson(*inputs: DNA, homology_length: int = 15):
    """Ligates pieces together with some homology length."""
    raise NotImplementedError
