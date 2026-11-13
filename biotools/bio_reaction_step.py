from __future__ import annotations

import biotools.bio_reactions as br
from biotools.bio_enums import BioReaction
from biotools.bio_exceptions import SimulationError
from biotools.dna import DNA
from biotools.oligo import Oligo
from typing import Callable, Any
import multiprocessing
import os
from itertools import product

# Callable function map
func_map: dict[BioReaction, Callable] = {
    BioReaction.AMPLIFY:    "amplify",
    BioReaction.ANNEAL:     "anneal_oligos",
    BioReaction.DIGEST:     "digest",
    BioReaction.GIBSON:     "gibson",
    BioReaction.GG:         "golden_gate",
    BioReaction.KLD:        "kld",
    BioReaction.LIGATE:     "ligate"
}

class BioReactionStep:
    """Class for representing reactions to be fed into a CompilationGraph.
    BioReactionSteps allow for multiprocessing of products with pools of sequences.
    BioReactionSteps will identify BioPools within the inputs and fork processes
    for each combination of pool.

    Parameters
    ----------
    type: BioReaction
        The type of reaction to simulate. One of the following:
            - Amplify       (PCR Reactions)
            - Anneal        (Oligo Annealing)
            - Digest        (Restriction Digest Reactions)
            - Gibson        (Gibson Assembly Reactions)
            - Golden Gate   (Golden Gate Assembly Reactions)
            - Input         (Individual Inputs)
            - KLD           (Kinase, Ligase, DpnI Reactions)
            - Ligate        (Ligation Reactions)
    name: str
        An identifying name for BioReactionStep.
        NOTE: In BioReactionGraphs, each Step must have a unique name
    kwargs: dict | DNA | None (Default: None)
        A dictionary of the keyword arguments used for the specific reaction

    Examples
    --------
    
    """

    def __init__(
        self,
        type: BioReaction,
        name: str,
        kwargs: dict[Any, Any] | DNA | None = None
    ):
        """Default constructor."""
        self.type = type
        self.name = name
        self.kwargs = kwargs
        self._has_pool = self._check_pools()

    def __eq__(self, oth: BioReactionStep) -> bool:
        return isinstance(oth, BioReactionStep) and self.name == oth.name and self.type == oth.type

    def __repr__(self) -> None:
        """Represented as BioReaction type and kwargs."""
        return f"{self.name}\n" \
                "-------\n" \
               f"KWARGS: {self._print_kwargs()}"

    def __hash__(self) -> int:
        return hash((self.name + self.type.name))

    def simulate(self, kwargs: dict[Any, Any] | None = None) -> DNA | list[Any]:
        """Simulates the reaction of self.type with self.kwargs.
        If there are pools present, will generate all combinations of inputs.

        Parameters
        ----------
        kwargs: dict[Any, Any] (Optional)
            A dictionary of arguments to pass into the reaction function of the Step

        Returns
        -------
        The product of the reaction, run as specified by self.kwargs or by the provided arguments
        """
        if self.type == BioReaction.INPUT:  # Input type returns its input
            return self.kwargs["input"]

        kargs: dict[str, Any] = {}
        if kwargs: kargs = kwargs
        else: kargs = self.kwargs
        if not kargs:
            raise SimulationError("Cannot simulate a reaction with no arguments provided!")

        self._has_pool = self._check_pools(kargs)  # Recheck for when we've added inputs in a graph

        # Generate arg list to create tuples for starmap
        arg_combos: list[tuple[Any]] = []
        func: Callable = getattr(br, func_map[self.type])

        if self._has_pool: # Perform pool workflow via multiprocessing
            match self.type:
                case BioReaction.AMPLIFY:
                    if kargs["template"].has_pool():
                        seqs: list[DNA] = kargs["template"].get_pools()
                    else:
                        seqs: list[DNA] = [kargs["template"]]

                    f_primers: list[Oligo] = []
                    for f_primer in kargs["f_primers"]:
                        if f_primer.has_pool():
                            f_primers.append(f_primer.get_pools())
                        else:
                            f_primers.append(f_primer)

                    r_primers: list[Oligo] = []
                    for r_primer in kargs["r_primers"]:
                        if r_primer.has_pool():
                            r_primers.append(r_primer.get_pools())
                        else:
                            r_primers.append(r_primer)

                    arg_combos = [
                        (f_primers, r_primers, seq, kargs["min_binding_length"], kargs["gel_extraction"]) \
                        for seq in seqs
                    ]
                case BioReaction.ANNEAL:
                    seqs: list[Any] = []
                    input_args: list[DNA | list[DNA]] = []

                    if kargs["f_oligo"].has_pool():
                        input_args.append(kargs["f_oligo"].get_pools())
                    else:
                        input_args.append([kargs["f_oligo"]])

                    if kargs["r_oligo"].has_pool():
                        input_args.append(kargs["r_oligo"].get_pools())
                    else:
                        input_args.append([kargs["r_oligo"]])

                    seqs: list[list[DNA | list[DNA]]] = product(*input_args)
                    arg_combos = [
                        (seqs[i][0], seqs[i][1]) for i in range(len(seqs))
                    ]
                case BioReaction.DIGEST:
                    # Check "input"
                    seqs = kargs["input"].get_pools()
                    arg_combos = [
                        (seq, kargs["enzymes"], kargs["gel_extraction"]) \
                        for seq in seqs
                    ]
                case BioReaction.GIBSON:
                    # Check "inputs"
                    input_args: list[DNA | list[DNA]] = []
                    for input in kargs["inputs"]:
                        if input.has_pool():
                            input_args.append(input.get_pools())
                        else:
                            input_args.append([input])

                    input_combos = product(*input_args)
                    arg_combos = [
                        (combo, kargs["min_homology_len"], kargs["max_homology_len"], kargs["gel_extraction"]) \
                        for combo in input_combos
                    ]
                case BioReaction.GG:
                    input_args: list[DNA | list[DNA]] = []
                    for input in kargs["inputs"]:
                        if input.has_pool():
                            input_args.append(input.get_pools())
                        else:
                            input_args.append([input])

                    input_combos = product(*input_args)
                    arg_combos = [
                        (combo, kargs["enzyme"], kargs["gel_extraction"]) \
                        for combo in input_combos
                    ]
                case BioReaction.KLD:
                    seqs: list[DNA] = kargs["input"].get_pools()
                    arg_combos = [
                        (seq,) for seq in seqs
                    ]
                case BioReaction.LIGATE:
                    # Check "inputs"
                    input_args: list[DNA | list[DNA]] = []
                    for input in kargs["inputs"]:
                        if input.has_pool():
                            input_args.append(input.get_pools())
                        else:
                            input_args.append([input])

                    input_combos = product(*input_args)
                    arg_combos = [
                        (combo, kargs["gel_extraction"], kargs["blunt"], kargs["allow_linear"]) \
                        for combo in input_combos
                    ]

            # Run and return combinations of args
            with multiprocessing.Pool(processes=os.cpu_count()-1) as pool:
                return pool.starmap(func, arg_combos)

        match self.type:  # Check if lists are present when they shouldn't be
            case BioReaction.ANNEAL:
                if isinstance(kargs["f_oligo"], list):
                    f_oligos = kargs["f_oligo"]
                else:
                    f_oligos = [kargs["f_oligo"]]
                if isinstance(kargs["r_oligo"], list):
                    r_oligos = kargs["r_oligo"]
                else:
                    r_oligos = [kargs["r_oligo"]]
                arg_combos = [
                    (f, r) for f in f_oligos for r in r_oligos
                ]
            case BioReaction.AMPLIFY:
                if isinstance(kargs["template"], list):
                    seqs: list[DNA] = kargs["template"]
                arg_combos = [
                    (kargs["f_primers"], kargs["r_primers"], seq, kargs["min_binding_length"], kargs["gel_extraction"]) \
                    for seq in seqs
                ]
            case BioReaction.DIGEST:
                # Check "input"
                if isinstance(kargs["input"], list):
                    seqs = kargs["input"]
                    arg_combos = [
                        (seq, kargs["enzymes"], kargs["gel_extraction"]) \
                        for seq in seqs
                    ]
            case BioReaction.GIBSON:
                # Check "inputs"
                input_args: list[DNA | list[DNA]] = []
                has_list: bool = False
                for input in kargs["inputs"]:
                    if isinstance(input, list):
                        input_args.append(input)
                        has_list = True
                    else:
                        input_args.append([input])
                if has_list:
                    input_combos = product(*input_args)
                    arg_combos = [
                        (combo, kargs["min_homology_len"], kargs["max_homology_len"], kargs["gel_extraction"]) \
                        for combo in input_combos
                    ]
            case BioReaction.GG:
                input_args: list[DNA | list[DNA]] = []
                has_list: bool = False
                for input in kargs["inputs"]:
                    if isinstance(input, list):
                        input_args.append(input)
                        has_list = True
                    else:
                        input_args.append([input])
                if has_list:
                    input_combos = product(*input_args)
                    arg_combos = [
                        (combo, kargs["enzyme"], kargs["gel_extraction"]) \
                        for combo in input_combos
                    ]
            case BioReaction.KLD:
                if isinstance(kargs["input"], list):
                    seqs: list[DNA] = kargs["input"]
                    arg_combos = [
                        (seq,) for seq in seqs
                    ]
            case BioReaction.LIGATE:
                # Check "inputs"
                input_args: list[DNA | list[DNA]] = []
                has_list: bool = False
                for input in kargs["inputs"]:
                    if isinstance(input, list):
                        input_args.append(input)
                        has_list = True
                    else:
                        input_args.append([input])

                if has_list:
                    input_combos = product(*input_args)
                    arg_combos = [
                        (combo, kargs["gel_extraction"], kargs["blunt"], kargs["allow_linear"]) \
                        for combo in input_combos 
                    ]

        if arg_combos:  # Multiprocessing workflow
            with multiprocessing.Pool(processes=os.cpu_count()-1) as pool:
                return pool.starmap(func, arg_combos)
        # We have no input lists where they usually aren't and no pools
        return func(**kargs)

    def _print_kwargs(self) -> str:
        """Generates printable dictionary in a pretty format for __repr__"""

        kwargs: str = "{\n"
        if self.kwargs:
            for bioreaction, func in self.kwargs.items():
                kwargs += f"         \"{bioreaction}\": {func}\n"
        kwargs += "        }\n"
        return kwargs

    def _check_pools(self, kwargs: dict[str, Any] | None = None) -> bool:
        """Returns True if pools exist in any input, False otherwise.
        If kwargs is provided, we check that instead of self.kwargs"""
        curr_kwargs: dict[str, Any]
        if not kwargs:
            if not self.kwargs:
                return False
            else:
                curr_kwargs = self.kwargs
        else:
            curr_kwargs = kwargs

        # We accept None for the inputs into the reactions and return False in that case
        match self.type:
            case BioReaction.ANNEAL:
                return curr_kwargs["f_primer"].has_pool() or \
                       curr_kwargs["r_primer"].has_pool()
            case BioReaction.AMPLIFY:
                if curr_kwargs["template"]:
                    if not curr_kwargs["template"].has_pool():
                        return any(p.has_pool() for p in curr_kwargs["f_primers"]) \
                            or any(p.has_pool() for p in curr_kwargs["r_primers"])
                    return True
                return any(p.has_pool() for p in curr_kwargs["f_primers"]) \
                    or any(p.has_pool() for p in curr_kwargs["r_primers"])
            case BioReaction.DIGEST | BioReaction.KLD:
                if not curr_kwargs["input"]: return False
                return curr_kwargs["input"].has_pool()
            case BioReaction.INPUT:
                return curr_kwargs["input"].has_pool()
            case BioReaction.GIBSON | BioReaction.GG:
                if not curr_kwargs["inputs"]:
                    return False
                for input in curr_kwargs["inputs"]:
                    if input.has_pool():
                        return True
            case BioReaction.LIGATE:
                if not curr_kwargs["inputs"]: return False
                for input in curr_kwargs["inputs"]:
                    if input.has_pool():
                        return True
        return False
