"""bio_reactions.py
Functions supporting reactions using DNA sequences.

Available Functions:
1. amplify()
2. kld()
3. digest()
4. ligate()
5. gibson()
6. golden_gate()
7. gel_extract()
8. anneal_oligos()
"""
from biotools.dna import DNA
from biotools.oligo import Oligo
from biotools.bio_exceptions import ReactionError
from biotools.bio_annotation import BioOrientation
from biotools.bio_enums import BioProperty
from biotools.bio_utils import (
    find_binding_site, rev_comp, find_re_sites, 
    validate_sites, check_homology, contains_re_site
)
from biotools.bio_alphabet import RE_ENZYMES

### Constants ###
MIN_PLASMID_LEN: int = 740

def amplify(
    f_primer: list[Oligo], 
    r_primer: list[Oligo], 
    template: DNA, 
    min_binding_length: int = 15,
    gel_extraction: list[tuple[int, int]] | None = None
) -> DNA | list[DNA]:
    """Amplifies DNA using provided primers that must bind at least min_binding_length bases.

    Parameters
    ----------
    f_primer: list[Oligo]
        A list of Oligos to use as forward primers. All Oligos in this list will be
        paired with all Oligos in r_primer
    r_primer: list[Oligo]
        A list of Oligos to use as reverse primers. All Oligos in this list will be
        paired with all Oligos in f_primer
    template: DNA
        A DNA to use as template. Must contain binding sites of at least length
        min_binding_length for the reaction to work
    min_binding_length: int [DEFAULT = 15]
        The minimum binding length allowed for each primer
    gel_extraction: tuple[int, int] (Optional)
        The tuple of lengths to extract

    Returns
    -------
    DNA or a list of DNA of products from the primer of pairing(s)

    Raises
    ------
    `ReactionError` if the length of any primer is less than min_binding_length,
    if any primer doesn't bind to template, if forward primers are downstream of
    reverse primers on non-circular templates, or if any forward primer overlaps
    with any reverse primer
    """
    if any(len(f) < min_binding_length for f in f_primer) or \
       any(len(r) < min_binding_length for r in r_primer):
        raise ReactionError(f"Minimum length of primer must be at least {min_binding_length} bp!")

    f_binding_spans: list[tuple[int, int]] = []
    for f in f_primer:
        site = find_binding_site(primer=f, template=template, single=True)
        if site:
            f_binding_spans.append(site)
        else:
            raise ReactionError(f"No binding sites found for f primer {f.name}!")

    r_binding_spans: list[tuple[int, int]] = []
    for r in r_primer:
        site = find_binding_site(primer=r, template=template.rev_comp(), single=True)
        if site:
            # Convert to forward span indices
            new_r_span = (template.length - site[1], template.length - site[0])
            r_binding_spans.append(new_r_span)
        else:
            raise ReactionError(f"No binding sites found for r primer {r.name}!")

    # Check if all f priming sites are before r priming sites and not less than min_binding_length away
    if not template.is_circular():
        if not all(f_span[1] - r_span[0] <= -min_binding_length for f_span in f_binding_spans for r_span in r_binding_spans):
            raise ReactionError("Forward primers must be upstream of reverse primers and must not overlap!")
    else:
        if not all(f_span[1] < r_span[0] for f_span in f_binding_spans for r_span in r_binding_spans):  # TODO: Can span[0] and span[1] be the same?
            raise ReactionError("Forward primers must not overlap with reverse primers!")

    final_products: list[DNA] = []
    # Generate PCR product(s)
    for i, f_span in enumerate(f_binding_spans):
        for j, r_span in enumerate(r_binding_spans):
            product = template[f_span[0]:r_span[1]]
            if f_span[1] - f_span[0] < f_primer[i].length: # Need to take difference in primer
                product = DNA(seq=f_primer[i][:f_primer[i].length - (f_span[1] - f_span[0])].seq) + product
            if r_span[1] - r_span[0] < r_primer[j].length:
                product = product + DNA(seq=rev_comp(r_primer[j][:r_primer[j].length - (r_span[1] - r_span[0])].seq))
            product.name = template.name+f"_PCR_Product_{i+j+1}"
            final_products.append(product)

    if gel_extraction: return gel_extract(*final_products, extraction_len_range=gel_extraction)

    return final_products[0] if len(final_products) == 1 else final_products

def kld(input: DNA) -> DNA:
    """Performs KLD on a linear DNA fragment.
    For more information, see
    https://www.neb.com/en-ca/protocols/kld-enzyme-mix-reaction-protocol-m0554

    Parameters
    ----------
    input: DNA
        The input to circularize

    Returns
    -------
    Circularized DNA from the KLD reaction

    Raises
    ------
    `ReactionError` if input is shorter than MIN_PLASMID_LEN or if input
    is circular
    """
    if input.length < MIN_PLASMID_LEN:  # NOTE: 740 is the shortest theoretical length possible for a viable plasmid
        raise ReactionError("Cannot generate plasmids of length shorter than 740!")
    if input.is_circular():
        raise ReactionError("Cannot perform KLD on a circular input!")
    return DNA(
        input.seq,
        input.name,
        input.type,
        BioProperty.CIRCULAR,
        input.strandedness,
        input.annotations,
        input.parent
    )

def digest(
    input: DNA,
    enzymes: list[str],
    gel_extraction: list[tuple[int, int]] | None = None
) -> DNA | list[DNA]:
    """Digests the input with the provided enzyme(s) and returns all products formed.
    Use gel_extract to grab products within a specific size range, inclusive.
    NOTE: Enzymes must be valid keys of RE_ENZYMES

    Parameters
    ----------
    input: DNA
        The DNA input to digest with the specified enzymes
    enzymes: list[str]
        A list of strings of enzymes to cut input with. Must be
        a valid key in RE_ENZYMES
    gel_extraction: list[tuple[int, int]] | None (Optional)
        A list of tuples of ranges of lengths at which to extract products from the
        reaction

    Returns
    -------
    A single DNA product or a list of products if more than one was created

    Raises
    ------
    `ReactionError` if an IndexError is caught, suggesting that an enzyme cut off
    of the input fragment
    """
    sites: dict = find_re_sites(input, *enzymes)

    circular: bool = True if input.is_circular() else False
    tuple_list: list[tuple[int, int, int, tuple[int, int], int]] = validate_sites(sites, len(input), circular)

    products: list[DNA] = []
    # Generate products based on spans
    product_num: int = 1
    for i, tup in enumerate(tuple_list):
        offset = tup[0]
        span = tup[1]
        enzyme = tup[2]
        try:
            if i == 0:
                if i == len(tuple_list) - 1:  # First and last
                    if input.is_circular():  # Reindex one product
                        product = input.reindex(span[0] + offset[0], False)
                        product.seq += product.seq[:offset[1]]  # One enzyme edge case, must add additional seq to look like product
                        product.circular = BioProperty.LINEAR
                        product.strandedness = BioProperty.CUT
                        if offset[-1] == BioOrientation.BOTTOM:
                            product.offsets = ((0, offset[1]), (offset[1], 0))
                        else:
                            product.offsets = ((offset[1], 0), (0, offset[1]))
                        product.name = input.name+f" Product {product_num}"  # The whole original plasmid (no slice)
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
                    ds_slice = next_span[0] + next_tup[0][0] + next_tup[0][1] if next_tup[-1] == \
                               BioOrientation.FORWARD else (next_span[1] - next_tup[0][0])

                    if input.is_circular():  # Generate a product that uses the next span
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
                    else:  # Generate two products
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
            elif i == len(tuple_list) - 1:  # Last
                if input.is_circular():  # Use first's span
                    first_tup = tuple_list[0]
                    first_span = first_tup[1]

                    # Generate offsets to use
                    us_slice = span[0] + offset[0] if tup[-1] == BioOrientation.FORWARD else span[1] - offset[0] - offset[1]
                    ds_slice = first_span[0] + first_tup[0][0] + first_tup[0][1] if first_tup[-1] == \
                               BioOrientation.FORWARD else (first_span[1] - first_tup[0][0])

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
                    product = input[span[1] - offset[0] - offset[1] :]
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
                ds_slice = next_span[0] + next_tup[0][0] + next_tup[0][1] if next_tup[-1] == \
                           BioOrientation.FORWARD else (next_span[1] - next_tup[0][0])

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
            raise ReactionError(f"{enzyme} cut off of the input fragment. Check enzyme and cleavage locations.")

    # Gel extraction
    if gel_extraction: return gel_extract(*products, extraction_len_range=gel_extraction)
    return products[0] if len(products) == 1 else products

def ligate(
    inputs: list[DNA],
    gel_extraction: list[tuple[int, int]] | None = None,
    blunt: bool = False,
    allow_linear: bool = False
) -> DNA | list[DNA]:
    """Ligates the provided inputs together. Will generate self-ligating products if possible.

    Parameters
    ----------
    inputs: list[DNA]
        A list of DNA inputs to ligate together
    gel_extraction: list[tuple[int, int]] | None (Optional)
        A list of tuples of ranges of lengths at which to extract products from the
        reaction
    blunt: bool
        Whether or not to allow blunt ligations
    allow_linear: bool
        Whether or not to allow linear products

    Returns
    -------
    DNA if one product is generated or a list of DNA products if multiple are generated

    Raises
    ------
    `ReactionError` if any input is circular or uncut or if no parts ligate to one another
    """
    for input in inputs:
        if input.is_circular():
            raise ReactionError(f"Cannot ligate circular input! {part} is circular.")
        if not input.is_cut():
            raise ReactionError("Cannot ligate uncut input!")

    all_parts: list[DNA] = sorted(list(set(inputs)))  # Sort for deterministic fragment generation, remove dups

    parts_dict: dict[DNA, list[tuple[DNA, int]]] = {}
    for part in all_parts:
        for p in all_parts:
            if part not in parts_dict:
                parts_dict[part] = []
            if p.get_overhangs()[0] == rev_comp(part.get_overhangs()[1]):  # p's 5' checking to part's 3' rev comp
                if p.get_overhangs()[0] == "":
                    if blunt:
                        parts_dict[part].append((p, len(p.get_overhangs()[0])))
                else:  # Non-blunt
                    if part not in parts_dict:
                        parts_dict[part] = []
                    parts_dict[part].append((p, len(p.get_overhangs()[0])))

    if not parts_dict: raise ReactionError("No parts ligate to one another!")

    products: list[list[DNA]] = _generate_all_combinations(parts_dict)
    product_list: list[set[DNA]] = []
    sanitized_list: list[list[DNA]] = []
    # Generate sets to grab unique ones
    for product_check in products:
        if set(product_check) not in product_list:
            sanitized_list.append(product_check)
            product_list.append(set(product_check))

    # Create the final Assemblies
    final_assemblies: list[DNA] = []
    new_seq: DNA = None
    prev_part: DNA = None
    new_name: list[str] = []

    for product_list in sanitized_list:
        if product_list[0] != product_list[-1]:  # Linear product
            if allow_linear:  # Must be linear here
                for i, part in enumerate(product_list):
                    if i == 0:
                        new_seq = part
                        new_name = [part.name]
                    else:
                        slice_idx: int = [p[1] for p in parts_dict[prev_part] if p[0] == part][0]
                        if i == len(product_list) - 1:  # Final part, must add
                            new_name += [part.name]
                            new_seq.name = f"{new_name} Ligation Product"
                            new_seq = new_seq.concatenate([part[slice_idx :]])
                            new_seq.strandedness = BioProperty.CUT
                            new_seq.circular = BioProperty.LINEAR
                            final_assemblies.append(new_seq)
                        else:
                            new_name += [part.name]
                            new_seq = new_seq.concatenate([part[slice_idx :]])
                    prev_part = part
        else:  # Circular product
            for i, part in enumerate(product_list):
                if i == 0:
                    new_seq = part
                    new_name = [part.name]
                else:
                    slice_idx: int = [p[1] for p in parts_dict[prev_part] if p[0] == part][0]
                    if i == len(product_list) - 1:  # Final part, don't add because it repeats first part
                        new_seq.name = f"{new_name} Ligation Product"
                        new_seq.strandedness = BioProperty.DOUBLE_STRANDED
                        new_seq.circular = BioProperty.CIRCULAR
                        new_seq.offsets = ((0, 0), (0, 0))
                        new_seq = new_seq[: len(new_seq) - slice_idx]  # Trim to remove slicing repeat of overlap
                        final_assemblies.append(new_seq)
                    else:
                        new_name += [part.name]
                        new_seq = new_seq.concatenate([part[slice_idx :]])
                prev_part = part

    # Gel extraction
    if gel_extraction: return gel_extract(*final_assemblies, extraction_len_range=gel_extraction)
    return final_assemblies[0] if len(final_assemblies) == 1 else final_assemblies

def _generate_all_combinations(
    parts_dict: dict[DNA, list[tuple[DNA, int]], list[tuple[DNA, int]]],
    start: DNA = None,
    path: list[DNA] = [] 
) -> list[list[DNA]]:
    """Generates all possible combinations of ligation products,
    stopping when a cycle is found
    """

    if start is None:
        all_paths: list[list[DNA]] = []
        for key in parts_dict.keys():
            all_paths.extend(_generate_all_combinations(parts_dict, key))
        return all_paths

    current_path: list[DNA] = path + [start]
    paths: list[list[DNA]] = []
    children: list[DNA] = [dna for dna, _ in parts_dict[start]]

    if len(current_path) > 1:
        # Check for cycle starting with first part
        if children:
            for c in children:
                if current_path[0] == c:
                    paths.append(current_path + [c])  # Add to the graph so we know it was a cycle
        else:  # No children, so linear
            paths.append(current_path)
            return paths

    if children:
        for child in children:
            if child not in current_path:  # Check each child of the current node, adding to a graph
                paths.extend(_generate_all_combinations(parts_dict, child, current_path))

    return paths

def gibson(
    inputs: list[DNA],
    min_homology_len: int = 20,
    max_homology_len: int = 40,
    gel_extraction: list[tuple[int, int]] | None = None
) -> DNA | list[DNA]:
    """Performs a Gibson cloning reaction on one or more products. The order of products does
    not matter, and homology will be determined. This will generate all possible products from the
    input parts and return either a single product or list of products.

    Parameters
    ----------
    inputs: list[DNA]
        A list of DNA inputs to gibson together
    min_homology_len: int
        The minimum length of homology any two products can have betweeen one another
    max_homology_len: int
        The maximum length of homology any two products can have betweeen one another
    gel_extraction: list[tuple[int, int]] | None (Optional)
        A list of tuples of ranges of lengths at which to extract products from the
        reaction 

    Returns
    -------
    DNA or a list of DNA if more than one product was formed in the reaction

    Raises
    ------
    `ValueError` if either min_homology_len or max_homology_len is less than 20
    `ReactionError` if any part is shorter than min_homology_len, any part is
    circular, or any part is lacking homology at one or both ends
    """
    if min_homology_len < 20 or max_homology_len < 20:
        raise ValueError("Must check for at least 20 bases of homology between products.")

    if max_homology_len < min_homology_len:
        raise ValueError("max_homology_len cannot be less than min_homology_len!")

    all_parts: list[DNA] = []
    for part in inputs:
        if len(part) < min_homology_len:
            raise ReactionError(f"Part {part} is shorter than min_homology_len ({min_homology_len})!")
        if part.is_circular():
            raise ReactionError(f"Unable to use circular products in a Gibson reaction! {part} is circular.")
        if part not in all_parts:  # In case we duplicate parts
            all_parts.append(part)
    all_parts = sorted(all_parts)  # Sort for deterministic fragment generation

    # Generate linking dictionary
    parts_dict: dict[str, list[tuple[DNA, int]], tuple[DNA, int]] = {}
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
            if parts_dict[assem][1][0][0] in used_assemblies:  # Cycle found not involving first product
                valid = False
                break
            oth_assem = parts_dict[assem][1][0][0]
            final_name.append(oth_assem.name)
            final_idx = parts_dict[assem][1][0][1]
            final_assembly = final_assembly.concatenate([oth_assem[final_idx:]])
            assem = oth_assem
            used_assemblies.add(oth_assem)

        if valid:
            if used_assemblies not in assembly_set_list:  # Ensure configuration is unique
                assembly_set_list.append(used_assemblies)
                final_assembly = final_assembly[: len(final_assembly) - final_idx]
                final_assembly.name = f"{final_name} Gibson Product"
                final_assembly.circular = BioProperty.CIRCULAR
                final_assemblies.append(final_assembly)

    # Gel extraction
    if gel_extraction: return gel_extract(*final_assemblies, extraction_len_range=gel_extraction)
    return final_assemblies[0] if len(final_assemblies) == 1 else final_assemblies

def golden_gate(
    inputs: list[DNA],
    enzyme: str,
    gel_extraction: list[tuple[int, int]] | None = None
) -> DNA | list[DNA]:
    """Performs a Golden Gate Reaction, digesting and ligating the provided input with the provided
    enzyme.
    NOTE: This reaction never returns the original input even though it is technically possible,
    just unfavored.

    Parameters
    ----------
    inputs: DNA | list[DNA]
        DNA or a list of DNA to 
    enzyme: str
        The enzyme to digest the input(s) with
    gel_extraction: list[tuple[int, int]] (Optional)
        Gel extraction parameters if gel extraction of a specific length range or ranges
        is desired

    Returns
    -------
    One or more products from digesting and ligating via IIS enzymes

    Raises
    ------
    `ReactionError` if the enzyme is not type IIS
    """

    if enzyme not in list(RE_ENZYMES.keys())[0:4]:  # Check if enzyme is IIS
        raise ReactionError("Cannot perform Golden Gate reaction without type IIS enzyme!")

    digest_outputs: list[DNA] = []
    for inp in inputs:
        digest_outputs.extend(digest(input=inp, enzymes=[enzyme]))
    ligation_outputs: DNA | list[DNA] = ligate(inputs=digest_outputs, gel_extraction=gel_extraction)
    # Get rid of anything that was there originally or has the enzyme in it
    unique_outputs: list[DNA] = [lo for lo in ligation_outputs if all(lo != inp for inp in inputs)]
    # Golden gate doesn't return outputs with enzyme in them
    unique_outputs = [op for op in unique_outputs if not contains_re_site(op, RE_ENZYMES[enzyme])]
    if len(unique_outputs) == 0:
        raise ReactionError("No products formed from the golden gate reaction!")
    return unique_outputs[0] if len(unique_outputs) == 1 else unique_outputs

def gel_extract(
    *inputs: DNA,
    extraction_len_range: list[tuple[int,int]]
) -> DNA | list[DNA]:
    """From the provided inputs, extracts inputs within the specified extraction range, inclusive.

    Parameters
    ----------
    *inputs: DNA
        One or more DNA objects from which to extract a specific length range
    extraction_len_range: list[tuple[int,int]]
        A list of tuples of ranges to extract

    Returns
    -------
    DNA or a list of DNA within the bounds of extraction_len_range tuples, inclusive

    Examples
    --------
    >>> dna1.length
    550
    >>> dna2.length
    750
    >>> gel_extract(dna1, dna2, [(500,600), (1000,1100)])
    dna1 # We get back dna1 because it is between 500 and 600 but not
         # dna2 because it isn't in either range.
    """
    extracted_inputs: list[DNA] = [input for input in inputs if any(input.length >= elr[0] and input.length <= elr[1] \
                                                                    for elr in extraction_len_range)]
    return extracted_inputs[0] if len(extracted_inputs) == 1 else extracted_inputs

def anneal_oligos(f_oligo: Oligo, r_oligo: Oligo) -> DNA:
    """Anneals two oligos together, stitching them at their homology."""
    raise NotImplementedError
