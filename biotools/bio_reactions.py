# Functions supporting reactions using DNA sequences
from dna import DNA
from oligo import Oligo
from bio_exceptions import ReactionError
from bio_annotation import BioOrientation
from bio_enums import *
from bio_utils import find_binding_site, rev_comp, find_re_sites, validate_sites, check_homology
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

def digest(input: DNA, enzymes: list[str], gel_extract: tuple[int, int] = None) -> DNA | list[DNA]:
    """Digests the input with the provided enzyme(s) and returns all products formed.
    Use gel_extract to grab products within a specific size range, inclusive.
    NOTE: Enzymes must be valid keys of bio_alphabet.re_enzymes"""
    sites: dict = find_re_sites(input, *enzymes)

    circular = True if input.circular == BioProperty.CIRCULAR else False
    if not (tuple_list := validate_sites(sites, len(input), circular)):
        raise ReactionError("Invalid enzyme parameters detected! Check cut sites do not overlap or cut other sites.")

    products: list[DNA] = []
    # Generate products based on spans
    product_num: int = 1
    for i, tup in enumerate(tuple_list):
        offset = tup[0]
        span = tup[1]
        enzyme = tup[2]
        try:
            if i == 0:
                if i == len(tuple_list) - 1: # First and last
                    if input.is_circular(): # Reindex one product
                        product = input.reindex(span[0] + offset[0], False)
                        product.seq += product.seq[:offset[1]] # One enzyme edge case, must add additional seq to look like product
                        product.circular = BioProperty.LINEAR
                        product.strandedness = BioProperty.CUT
                        if offset[-1] == BioOrientation.BOTTOM:
                            product.offsets = ((0, offset[1]), (offset[1], 0))
                        else:
                            product.offsets = ((offset[1], 0), (0, offset[1]))
                        product.name = input.name+f" Product {product_num}" # The whole original plasmid (no slice)
                        product.parent = input
                        product_num += 1
                        products.append(product)
                    else: # Make two products
                        product_1 = input[:span[0] + offset[0]]
                        product_1.seq += input.seq[span[0]+offset[0]:span[0]+offset[0]+offset[1]]
                        product_1.strandedness = BioProperty.CUT
                        if offset[-1] == BioOrientation.BOTTOM:
                            product_1.offsets = ((0, offset[1]), (0, 0))
                        else:
                            product_1.offsets = ((offset[1], 0), (0, 0))
                        product_1.name = input.name+f" Product {product_num} [0, {span[0] + offset[0]}]"
                        product_1.parent = input
                        product_num += 1
                        product_2 = input[span[0] + offset[0]:]
                        product_2.strandedness = BioProperty.CUT
                        if offset[-1] == BioOrientation.BOTTOM:
                            product_2.offsets = ((0, 0), (offset[1], 0))
                        else:
                            product_2.offsets = ((0, 0), (0, offset[1]))
                        product_2.name = input.name+f" Product {product_num} [{span[0] + offset[0]}, {len(input)}]"
                        product_2.parent = input
                        product_num += 1
                        products.extend([product_1, product_2])
                else: # Just first with at least one more
                    next_tup = tuple_list[i + 1]
                    next_span = next_tup[1]

                    # Generate offsets to use
                    us_slice = span[0] + offset[0] if tup[-1] == BioOrientation.FORWARD else span[1] - offset[0] - offset[1]
                    ds_slice = next_span[0] + next_tup[0][0] + next_tup[0][1] if next_tup[-1] == BioOrientation.FORWARD else (next_span[1] - next_tup[0][0])

                    if input.is_circular(): # Generate a product that uses the next span
                        product = input[us_slice:ds_slice]
                        product.circular = BioProperty.LINEAR
                        product.strandedness = BioProperty.CUT
                        if tuple_list[i + 1][0][-1] == BioOrientation.BOTTOM:
                            if offset[-1] == BioOrientation.BOTTOM:
                                product.offsets = ((0, tuple_list[i + 1][0][1]), (offset[1], 0))
                            else:
                                product.offsets = ((offset[1], tuple_list[i + 1][0][1]), (0, 0))
                        else:
                            if offset[-1] == BioOrientation.BOTTOM:
                                product.offsets = ((0, 0), (offset[1], tuple_list[i + 1][0][1]))
                            else:
                                product.offsets = ((offset[1], 0), (0, tuple_list[i + 1][0][1]))
                        product.name = input.name+f" Product {product_num} [{us_slice}, {ds_slice}]"
                        product.parent = input
                        product_num += 1
                        products.append(product)
                    else: # Generate two products
                        product_1 = input[:us_slice+offset[1]]
                        product_1.strandedness = BioProperty.CUT
                        if offset[-1] == BioOrientation.BOTTOM:
                            product_1.offsets = ((0, offset[1]), (0, 0))
                        else:
                            product_1.offsets = ((offset[1], 0), (0, 0))
                        product_1.name = input.name+f" Product {product_num} [0, {us_slice}]"
                        product_1.parent = input
                        product_num += 1
                        product_2 = input[us_slice:ds_slice]
                        product_2.strandedness = BioProperty.CUT
                        if tuple_list[i + 1][0][-1] == BioOrientation.BOTTOM:
                            if offset[-1] == BioOrientation.BOTTOM:
                                product_2.offsets = ((0, tuple_list[i + 1][0][1]), (offset[1], 0))
                            else:
                                product_2.offsets = ((offset[1], tuple_list[i + 1][0][1]), (0, 0))
                        else:
                            if offset[-1] == BioOrientation.BOTTOM:
                                product_2.offsets = ((0, 0), (offset[1], tuple_list[i + 1][0][1]))
                            else:
                                product_2.offsets = ((offset[1], 0), (0, tuple_list[i + 1][0][1]))
                        product_2.name = input.name+f" Product {product_num} [{us_slice}, {ds_slice}]"
                        product_2.parent = input
                        product_num += 1
                        products.extend([product_1, product_2])
            elif i == len(tuple_list) - 1: # Last
                if input.is_circular(): # Use first's span
                    first_tup = tuple_list[0]
                    first_span = first_tup[1]

                    # Generate offsets to use
                    us_slice = span[0] + offset[0] if tup[-1] == BioOrientation.FORWARD else span[1] - offset[0] - offset[1]
                    ds_slice = first_span[0] + first_tup[0][0] + first_tup[0][1] if first_tup[-1] == BioOrientation.FORWARD else (first_span[1] - first_tup[0][0])

                    product = input[us_slice:ds_slice]
                    product.circular = BioProperty.LINEAR
                    product.strandedness = BioProperty.CUT
                    if tuple_list[0][0][-1] == BioOrientation.BOTTOM:
                        if offset[-1] == BioOrientation.BOTTOM:
                            product.offsets = ((0, tuple_list[0][0][1]), (offset[1], 0))
                        else:
                            product.offsets = ((offset[1], tuple_list[0][0][1]), (0, 0))
                    else:
                        if offset[-1] == BioOrientation.BOTTOM:
                            product.offsets = ((0, 0), (offset[1], tuple_list[0][0][1]))
                        else:
                            product.offsets = ((offset[1], 0), (0, tuple_list[0][0][1]))
                    product.name = input.name+f" Product {product_num} [{us_slice}, {ds_slice}]"
                    product.parent = input
                    product_num += 1
                    products.append(product)

                else: # Go to end
                    product = input[span[0] + offset[0]:]
                    product.strandedness = BioProperty.CUT
                    if offset[-1] == BioOrientation.BOTTOM:
                        product.offsets = ((0, 0), (offset[1], 0))
                    else:
                        product.offsets = ((0, 0), (0, offset[1]))
                    product.name = input.name+f" Product {product_num} [{span[0] + offset[0]}, {len(input)}]"
                    product.parent = input
                    product_num += 1
                    products.append(product)
            else: # Neither first nor last
                next_tup = tuple_list[i + 1]
                next_span = next_tup[1]

                # Generate offsets to use
                us_slice = span[0] + offset[0] if tup[-1] == BioOrientation.FORWARD else span[1] - offset[0] - offset[1]
                ds_slice = next_span[0] + next_tup[0][0] + next_tup[0][1] if next_tup[-1] == BioOrientation.FORWARD else (next_span[1] - next_tup[0][0])

                product = input[us_slice:ds_slice]
                product.circular = BioProperty.LINEAR
                product.strandedness = BioProperty.CUT
                if tuple_list[i + 1][0][-1] == BioOrientation.BOTTOM:
                    if offset[-1] == BioOrientation.BOTTOM:
                        product.offsets = ((0, tuple_list[i + 1][0][1]), (offset[1], 0))
                    else:
                        product.offsets = ((offset[1], tuple_list[i + 1][0][1]), (0, 0))
                else:
                    if offset[-1] == BioOrientation.BOTTOM:
                        product.offsets = ((0, 0), (offset[1], tuple_list[i + 1][0][1]))
                    else:
                        product.offsets = ((offset[1], 0), (0, tuple_list[i + 1][0][1]))
                product.name = input.name+f" Product {product_num} [{us_slice}, {ds_slice}]"
                product.parent = input
                product_num += 1
                products.append(product)
        except IndexError:
            raise ReactionError(f"{enzyme} cut off of the input fragment. Check enzyme and cleavage location.")
        
    # Gel extraction
    if gel_extract:
        extracted_products: list[DNA] = [product for product in products if (product.length >= gel_extract[0] and product.length <= gel_extract[1])]
        if len(extracted_products) == 1: return extracted_products[0]
        return extracted_products

    if len(products) == 1: return products[0]
    return products

def ligate(*inputs: DNA) -> DNA:
    """Ligates together the provided inputs."""
    # TODO: How to handle blunt ligation? Maybe separate flag?
    for input in inputs:
        if input.is_circular():
            raise ReactionError("Cannot ligate circular input!")
        if not input.is_cut():
            raise ReactionError("Cannot ligate uncut input!")
    
    all_parts: list[DNA] = []
    for part in inputs:
        if part.is_circular():
            raise ReactionError(f"Unable to use circular products in a Gibson reaction! {part} is circular.")
        if part not in all_parts:  # In case we duplicate parts
            all_parts.append(part)
    all_parts = sorted(all_parts)  # Sort for deterministic fragment generation

    raise NotImplementedError

def gel_extract(*inputs: DNA, extraction_len_range: tuple[int,int]) -> DNA | list[DNA]:
    raise NotImplementedError

def gibson(*inputs: DNA, min_homology_len: int = 20, max_homology_len: int = 40, gel_extract: tuple[int, int] = None) -> DNA | list[DNA]:
    """Performs a Gibson cloning reaction on one or more products. The order of products does
    not matter, and homology will be determined. This will generate all possible products from the 
    input parts and return either a single product or list of products.
    """
    if min_homology_len < 20 or max_homology_len < 20:
        raise ValueError("Must check for at least 20 bases of homology between products.")

    all_parts: list[DNA] = []
    for part in inputs:
        if part.is_circular():
            raise ReactionError(f"Unable to use circular products in a Gibson reaction! {part} is circular.")
        if part not in all_parts:  # In case we duplicate parts
            all_parts.append(part)
    all_parts = sorted(all_parts)  # Sort for deterministic fragment generation

    parts_dict: dict[str, list[tuple[int, DNA]]] = {}
    for part in all_parts:
        parts_dict[part] = [[], []]
        for p in all_parts:
            if (h_idx := check_homology(part, p, min_homology_len, max_homology_len)) is not None:  # 5' checking
                parts_dict[part][0].append((p, h_idx))
            if (h_idx := check_homology(p, part, min_homology_len, max_homology_len)) is not None:  # 3' checking
                parts_dict[part][1].append((p, h_idx))

    # Check that everything matches at both ends
    for k in parts_dict.keys():
        if not parts_dict[k][0] or not parts_dict[k][1]:
            raise ReactionError(f"Input {k.name} did not have homology at one or both ends.")

    # Create the final Assembly
    final_assemblies: list[DNA] = []
    assembly_set_list: list[set[DNA]] = []
    for part in all_parts:
        final_name = [part.name]
        valid: bool = True
        used_assemblies: set[DNA] = {part}
        first_assem, assem = part, part
        final_assembly = first_assem
        final_idx: int = 0
        while not any(asmb[0] == first_assem for asmb in parts_dict[assem][1]):
            if parts_dict[assem][1][0][0] in used_assemblies:  # Cycle found
                valid = False
                break
            oth_assem = parts_dict[assem][1][0][0]
            final_name.append(oth_assem.name)
            i = parts_dict[assem][1][0][1]
            final_assembly = final_assembly.concatenate([oth_assem[i:]])
            assem = oth_assem
            used_assemblies.add(oth_assem)
            final_idx = i

        if valid:
            if used_assemblies not in assembly_set_list:  # Ensure configuration is unique
                assembly_set_list.append(used_assemblies)
                final_assembly = final_assembly[: len(final_assembly) - final_idx]
                final_assembly.name = f"{set(final_name)} Gibson Product"
                final_assembly.circular = BioProperty.CIRCULAR
                final_assemblies.append(final_assembly)

    # Gel extraction
    if gel_extract: 
        extracted_products: list[DNA] = [assembly for assembly in final_assemblies if (assembly.length >= gel_extract[0] and assembly.length <= gel_extract[1])]
        if len(extracted_products) == 1: return extracted_products[0]
        return extracted_products

    if len(final_assemblies) == 1:
        return final_assemblies[0]
    return final_assemblies

