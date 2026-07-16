from __future__ import annotations

from graphviz import Digraph
import networkx as nx
import pydot
from biotools.bio_reaction_step import BioReactionStep
from biotools.bio_exceptions import (
    InvalidInstantiationException, ReactionError, SimulationError
)
from typing import Any

class BioReactionGraph:
    """Class for representing and executing reactions via specified, sequential steps.
    Compilation workflows can be generated and simulated to check for errors or
    inconsistencies.
    
    Parameters
    ----------
    *steps: tuple[BioReactionStep, BioReactionStep]
        One or more tuples of BioReactionSteps to generate edges between for the
        simulation graph
    name: str | None (Optional)
        An identifying name for a BioReactionGraph
    """
    
    def __init__(
        self,
        *steps: tuple[BioReactionStep, BioReactionStep],
        name: str | None = "BioReactionGraph"
    ):
        """Default constructor."""
        self.steps = [step for step in steps]
        self.name = name
        self.graph, self.graph_dict = self._generate_graph()

    def __repr__(self) -> str:
        raise NotImplementedError

    def show_graph(self) -> None:
        """Renders the visualization of the graph."""
        return self.graph

    def simulate(self, *steps: BioReactionStep | None) -> bool:
        """Simulates the steps of the BioReactionGraph.
        If one or more steps are provided, then only these steps will be simulated.

        Parameters
        ----------
        *steps: BioReactionStep | None (Optional)
            One or more BioReactionStep objects to simulate.
        
        Returns
        -------
        True if simulation was successful, raises an exception otherwise."""

        keys: list[Any] = list(self.graph_dict.keys())
        values: list[Any] = list(self.graph_dict.values())
        # Simulate steps if provided
        if steps:
            for step in steps:
                if step not in keys and step not in values:
                    raise ValueError(f"Could not find Step {step.name}!")
                else:
                    try:
                        step.simulate()
                        self.graph.node(step.name, fillcolor="#90EE90", style="filled")
                    except (ValueError, ReactionError) as e:
                        raise SimulationError(f"Error {type(e).__name__} caught during simulating step {step.name}!")
                    
        return True

    def _generate_graph(self) -> tuple[Digraph, dict[BioReactionGraph, set[BioReactionGraph]]]:
        """Generates a graph based on the input steps."""
        
        g = Digraph(format="svg")
        g.attr(label=self.name, labelloc="t", fontsize="15")

        graph_dict: dict[BioReactionGraph, set[BioReactionGraph]] = {}

        for step in self.steps:
            if step[0] not in graph_dict.keys():
                graph_dict[step[0]] = set()
            g.node(step[0].name, fillcolor="#FFCCCB", style="filled")
            g.node(step[1].name, fillcolor="#FFCCCB", style="filled")
            g.edge(step[0].name, step[1].name)
            graph_dict[step[0]].add(step[1])

        # Verify graph is a DAG
        if not nx.is_directed_acyclic_graph(nx.nx_pydot.from_pydot(pydot.graph_from_dot_data(g.source)[0])):
            raise InvalidInstantiationException(f"{self.name} is not a directed acyclic graph.\n" \
                                                "Please ensure no products eventually feed into themselves.")

        return g, graph_dict
