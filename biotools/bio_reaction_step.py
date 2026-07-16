from __future__ import annotations

from biotools.bio_reactions import (
    amplify, digest, ligate, gibson, anneal_oligos, kld
)
from biotools.bio_enums import BioReaction
from biotools.dna import DNA
from typing import Callable, Any
import multiprocessing
import os
from itertools import product

# Callable function map 
func_map: dict[BioReaction, Callable] = {
    BioReaction.AMPLIFY:    amplify,
    BioReaction.ANNEAL:     anneal_oligos,
    BioReaction.DIGEST:     digest,
    BioReaction.GIBSON:     gibson,
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
            - Amplify   (PCR Reactions)
            - Anneal    (Oligo Annealing)
            - Digest    (Restriction Digest Reactions)
            - Gibson    (Gibson Assembly Reactions)
            - KLD       (Kinase, Ligase, DpnI Reactions)
            - Ligate    (Ligation Reactions)
    kwargs: dict
        A dictionary of the keyword arguments used for the specific reaction
    name: str
        An identifying name for BioReactionStep.
        NOTE: In BioReactionGraphs, each Step must have a unique name
    """
    
    def __init__(
        self,
        type: BioReaction,
        kwargs: dict,
        name: str
    ):
        """Default constructor."""
        self.type = type
        self.kwargs = kwargs
        self.name = name
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
        return hash(self.__repr__())
    
    def simulate(self) -> DNA | list[Any]:
        """Simulates the reaction of self.type with self.kwargs.
        If there are pools present, will generate all combinations of inputs.
        
        Returns
        -------
        The product of the reaction, run as specified by self.kwargs
        """

        if self._has_pool: # Perform pool workflow via multiprocessing
            # Generate arg list to create tuples for starmap
            arg_combos: list[tuple[Any]] = []
            
            match self.type:
                case BioReaction.AMPLIFY:
                    seqs = self.kwargs["template"].get_pools()
                    arg_combos = [
                        (self.kwargs["f_primer"], self.kwargs["r_primer"], seqs[i], self.kwargs["min_binding_length"]) \
                        for i in range(len(seqs)) 
                    ]
                case BioReaction.ANNEAL:
                    seqs: list[Any] = []
                    input_args: list[DNA | list[DNA]] = []
                    if self.kwargs["f_primer"].has_pool():
                        input_args.append(self.kwargs["f_primer"].get_pools())
                    else:
                        input_args.append([self.kwargs["f_primer"]])
                    if self.kwargs["r_primer"].has_pool():
                        input_args.append(self.kwargs["r_primer"].get_pools())
                    else:
                        input_args.append([self.kwargs["f_primer"]])
                    
                    seqs = list(product(*input_args))

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
                        (input_combos[i], self.kwargs["gel_extraction"]) \
                        for i in range(len(input_combos)) 
                    ]

            # Run and return combinations of args
            with multiprocessing.Pool(processes=os.cpu_count()-1) as pool:
                return pool.starmap(func_map[self.type], arg_combos)

        return func_map[self.type](**self.kwargs)
    
    def _print_kwargs(self) -> str:
        """Generates printable dictionary in a pretty format for __repr__"""

        kwargs: str = "{\n"
        for bioreaction, func in self.kwargs.items():
            kwargs += f"         \"{bioreaction}\": {func}\n"
        kwargs += "        }\n"
        return kwargs
    
    def _check_pools(self) -> bool:
        """Returns True if pools exist in any input, False otherwise."""

        if self.type == BioReaction.AMPLIFY:
            return self.kwargs["template"].has_pool()
        elif self.type == BioReaction.DIGEST or self.type == BioReaction.KLD:
            return self.kwargs["input"].has_pool()
        elif self.type == BioReaction.GIBSON:
            for input in self.kwargs["inputs"]:
                if input.has_pool():
                    return True
                return False
        elif self.type == BioReaction.LIGATE:
            for input in self.kwargs["inputs"]:
                if input.has_pool():
                    return True
        return False
    