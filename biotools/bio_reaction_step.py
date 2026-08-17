from __future__ import annotations

from biotools.bio_reactions import (
    amplify, digest, ligate, gibson, anneal_oligos, kld,
    golden_gate
)
from biotools.bio_enums import BioReaction
from biotools.bio_exceptions import SimulationError
from biotools.dna import DNA
from typing import Callable, Any
import multiprocessing
import os
from itertools import product
from hashlib import sha256

# Callable function map 
func_map: dict[BioReaction, Callable] = {
    BioReaction.AMPLIFY:    amplify,
    BioReaction.ANNEAL:     anneal_oligos,
    BioReaction.DIGEST:     digest,
    BioReaction.GIBSON:     gibson,
    BioReaction.GG:         golden_gate,
    BioReaction.KLD:        kld,
    BioReaction.LIGATE:     ligate
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
        if not isinstance(oth, BioReactionStep):
            return False
        if self.name != oth.name:
            return False
        if self.kwargs != oth.kwargs:
            return False
        if self.type != oth.type:
            return False
        if self._has_pool != oth._has_pool:
            return False
        return True

    def __repr__(self) -> None:
        """Represented as BioReaction type and kwargs."""
        return f"{self.name}\n" \
                "-------\n" \
               f"KWARGS: {self._print_kwargs()}"

    def __hash__(self) -> int:
        return int(sha256(self.__repr__().encode()).hexdigest(), 16)

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

        self.kwargs: dict[str, Any] = {}
        if kwargs: self.kwargs = kwargs
        if not self.kwargs: 
            raise SimulationError("Cannot simulate a reaction with no arguments provided!")

        self._has_pool = self._check_pools()  # Recheck for when we've added inputs in a graph

        # Generate arg list to create tuples for starmap
        arg_combos: list[tuple[Any]] = []

        if self._has_pool: # Perform pool workflow via multiprocessing
            match self.type:
                case BioReaction.AMPLIFY:
                    if isinstance(self.kwargs["template"], DNA):
                        seqs: list[DNA] = self.kwargs["template"].get_pools()
                    elif isinstance(self.kwargs["template"], list):
                        seqs: list[DNA] = self.kwargs["template"]
                    arg_combos = [
                        (self.kwargs["f_primer"], self.kwargs["r_primer"], seqs[i], self.kwargs["min_binding_length"], self.kwargs["gel_extraction"]) \
                        for i in range(len(seqs))
                    ]
                case BioReaction.ANNEAL:
                    seqs: list[Any] = []
                    input_args: list[DNA | list[DNA]] = []

                    if self.kwargs["f_oligo"].has_pool():
                        input_args.append(self.kwargs["f_oligo"].get_pools())
                    else:
                        input_args.append([self.kwargs["f_oligo"]])

                    if self.kwargs["r_oligo"].has_pool():
                        input_args.append(self.kwargs["r_oligo"].get_pools())
                    else:
                        input_args.append([self.kwargs["r_oligo"]])

                    seqs: list[list[DNA | list[DNA]]] = list(product(*input_args))
                    arg_combos = [
                        (seqs[i][0], seqs[i][1]) for i in range(len(seqs))
                    ]
                case BioReaction.DIGEST:
                    # Check "input"
                    seqs = self.kwargs["input"].get_pools()
                    arg_combos = [
                        (seqs[i], self.kwargs["enzymes"], self.kwargs["gel_extraction"]) \
                        for i in range(len(seqs))
                    ]
                case BioReaction.GIBSON:
                    # Check "inputs"
                    input_args: list[DNA | list[DNA]] = []
                    for input in self.kwargs["inputs"]:
                        if input.has_pool():
                            input_args.append(input.get_pools())
                        else:
                            input_args.append([input])

                    input_combos = list(product(*input_args))
                    arg_combos = [
                        (input_combos[i], self.kwargs["min_homology_len"], self.kwargs["max_homology_len"], self.kwargs["gel_extraction"]) \
                        for i in range(len(input_combos))
                    ]
                case BioReaction.GG:
                    input_args: list[DNA | list[DNA]] = []
                    for input in self.kwargs["inputs"]:
                        if input.has_pool():
                            input_args.append(input.get_pools())
                        else:
                            input_args.append([input])

                    input_combos = list(product(*input_args))
                    arg_combos = [
                        (input_combos[i], self.kwargs["enzyme"], self.kwargs["gel_extraction"]) \
                        for i in range(len(input_combos)) 
                    ]
                case BioReaction.KLD:
                    seqs: list[DNA] = self.kwargs["input"].get_pools()
                    arg_combos = [
                        (seqs[i],) for i in range(len(seqs))
                    ]
                case BioReaction.LIGATE:
                    # Check "inputs"
                    input_args: list[DNA | list[DNA]] = []
                    for input in self.kwargs["inputs"]:
                        if input.has_pool():
                            input_args.append(input.get_pools())
                        else:
                            input_args.append([input])

                    input_combos = list(product(*input_args))
                    arg_combos = [
                        (input_combos[i], self.kwargs["gel_extraction"], self.kwargs["blunt"], self.kwargs["allow_linear"]) \
                        for i in range(len(input_combos)) 
                    ]

            # Run and return combinations of args
            with multiprocessing.Pool(processes=os.cpu_count()-1) as pool:
                return pool.starmap(func_map[self.type], arg_combos)

        match self.type:  # Check if lists are present when they shouldn't be.
            case BioReaction.AMPLIFY:  # NOTE: We don't consider ANNEAL here
                if isinstance(self.kwargs["template"], list):
                    seqs: list[DNA] = self.kwargs["template"]
                arg_combos = [
                    (self.kwargs["f_primer"], self.kwargs["r_primer"], seqs[i], self.kwargs["min_binding_length"], self.kwargs["gel_extraction"]) \
                    for i in range(len(seqs))
                ]
            case BioReaction.DIGEST:
                # Check "input"
                if isinstance(self.kwargs["input"], list):
                    seqs = self.kwargs["input"]
                    arg_combos = [
                        (seqs[i], self.kwargs["enzymes"], self.kwargs["gel_extraction"]) \
                        for i in range(len(seqs))
                    ]
            case BioReaction.GIBSON:
                # Check "inputs"
                input_args: list[DNA | list[DNA]] = []
                has_list: bool = False
                for input in self.kwargs["inputs"]:
                    if isinstance(input, list):
                        input_args.append(input)
                        has_list = True
                    else:
                        input_args.append([input])
                if has_list:
                    input_combos = list(product(*input_args))
                    arg_combos = [
                        (input_combos[i], self.kwargs["min_homology_len"], self.kwargs["max_homology_len"], self.kwargs["gel_extraction"]) \
                        for i in range(len(input_combos))
                    ]
            case BioReaction.GG:
                input_args: list[DNA | list[DNA]] = []
                has_list: bool = False
                for input in self.kwargs["inputs"]:
                    if isinstance(input, list):
                        input_args.append(input)
                        has_list = True
                    else:
                        input_args.append([input])
                if has_list:
                    input_combos = list(product(*input_args))
                    arg_combos = [
                        (input_combos[i], self.kwargs["enzyme"], self.kwargs["gel_extraction"]) \
                        for i in range(len(input_combos))
                    ]
            case BioReaction.KLD:
                if isinstance(self.kwargs["input"], list):
                    seqs: list[DNA] = self.kwargs["input"]
                    arg_combos = [
                        (seqs[i],) for i in range(len(seqs))
                    ]
            case BioReaction.LIGATE:
                # Check "inputs"
                input_args: list[DNA | list[DNA]] = []
                has_list: bool = False
                for input in self.kwargs["inputs"]:
                    if isinstance(input, list):
                        input_args.append(input)
                        has_list = True
                    else:
                        input_args.append([input])

                if has_list:
                    input_combos = list(product(*input_args))
                    arg_combos = [
                        (input_combos[i], self.kwargs["gel_extraction"], self.kwargs["blunt"], self.kwargs["allow_linear"]) \
                        for i in range(len(input_combos)) 
                    ]

        if arg_combos:  # Multiprocessing workflow
            with multiprocessing.Pool(processes=os.cpu_count()-1) as pool:
                return pool.starmap(func_map[self.type], arg_combos)
        # We have no input lists where they usually aren't and no pools
        return func_map[self.type](**self.kwargs)

    def _print_kwargs(self) -> str:
        """Generates printable dictionary in a pretty format for __repr__"""

        kwargs: str = "{\n"
        if self.kwargs:
            for bioreaction, func in self.kwargs.items():
                kwargs += f"         \"{bioreaction}\": {func}\n"
        kwargs += "        }\n"
        return kwargs

    def _check_pools(self) -> bool:
        """Returns True if pools exist in any input, False otherwise."""
        if not self.kwargs:
            return False

        # We accept None for the inputs into the reactions and return False
        if self.type == BioReaction.AMPLIFY:
            if not self.kwargs["template"]: return False
            return self.kwargs["template"].has_pool()
        elif self.type == BioReaction.DIGEST or self.type == BioReaction.KLD:
            if not self.kwargs["input"]: return False
            return self.kwargs["input"].has_pool()
        elif self.type == BioReaction.INPUT:
            return self.kwargs["input"].has_pool()
        elif self.type == BioReaction.GIBSON or self.type == BioReaction.GG:
            if not self.kwargs["inputs"]:
                return False
            for input in self.kwargs["inputs"]:
                if input.has_pool():
                    return True
        elif self.type == BioReaction.LIGATE:
            if not self.kwargs["inputs"]: return False
            for input in self.kwargs["inputs"]:
                if input.has_pool():
                    return True
        return False
