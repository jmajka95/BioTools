from __future__ import annotations

from collections import deque
from typing import Any

import networkx as nx
import pydot
from biotools.bio_enums import BioReaction
from biotools.bio_exceptions import (
    InvalidInstantiationException, ReactionError, SimulationError
)
from biotools.bio_reaction_step import BioReactionStep
from biotools.dna import DNA
from graphviz import Digraph
from collections import defaultdict
from IPython.display import display
import copy


class BioReactionGraph:
    """Class for representing and executing reactions via specified, sequential steps.
    Compilation workflows can be generated and simulated to check for errors or
    inconsistencies.

    Parameters
    ----------
    *edges: tuple[BioReactionStep, BioReactionStep]
        One or more tuples of BioReactionSteps to generate edges between for the
        simulation graph
    name: str (Default: "")
        An identifying name for a BioReactionGraph
    """

    def __init__(
        self,
        *edges: tuple[BioReactionStep, BioReactionStep],
        name: str = ""
    ):
        """Default constructor"""
        self.edges = [edge for edge in edges]
        self.steps = set(step for edge in self.edges for step in edge)
        self.name = name
        self.graph, self.graph_dict, self.nx_graph, self.step_dict = self._generate_graph()

    def __repr__(self) -> str:
        return f"{{ {self.name} | BioReactionGraph | [{len(self.steps)} Node(s)] }}"

    def set_name(self, name: str) -> None:
        """Sets the name of the graph
        
        Parameters
        ----------
        name: str
            The new name of the graph
        """
        self.graph.attr(label=name)

    def show_graph(self) -> None:
        """Displays the graph"""
        display(self.graph)

    def simulate(self, *steps: BioReactionStep | None) -> DNA | list[DNA]:
        """Simulates the steps of the BioReactionGraph.
        If one or more steps are provided, then only these steps will be simulated on
        their own

        Parameters
        ----------
        *steps: BioReactionStep | None (Optional)
            One or more BioReactionSteps to simulate.

        Returns
        -------
        Products from the reaction simulation if the simulation was successful,
        raises a SimulationError otherwise.
        """

        output: DNA | list[DNA] = None
        keys: list[Any] = list(self.graph_dict.keys())
        values: list[Any] = list(self.graph_dict.values())
        all_values: list[Any] = [v for val in values for v in val]
        step_dict = copy.deepcopy(self.step_dict)  # Copy to avoid issues running simulate() many times
        if steps:  # Simulate steps if provided
            for step in steps:
                if step not in keys and step not in all_values:  # Check the value(s) of each set
                    raise ValueError(f"Could not find Step '{step.name}'!")
                else:
                    try:
                        step.simulate()
                        if step.type != BioReaction.INPUT:
                            self.graph.node(step.name, fillcolor="#90EE90", style="filled")
                    except (ValueError, ReactionError) as e:
                        raise SimulationError(f"The following {type(e).__name__} was caught during simulating step {step.name}:\n{e}")
        else:  # Simulate steps first to last, passing outputs into the next step
            # Grab all nodes with no parents
            initialization_list: list[BioReactionStep] = []
            for node in self.graph_dict.keys():
                if not list(self.nx_graph.predecessors(node.name)):
                    initialization_list.append(node)
            reaction_queue: deque = deque(initialization_list)

            curr_step: BioReactionStep
            output: DNA | list[DNA]
            while reaction_queue:
                curr_step = reaction_queue.popleft()
                try:
                    if curr_step in initialization_list:  # Simulate every node w/ no predecessor
                        output = curr_step.simulate()
                    else:
                        output = curr_step.simulate(kwargs=step_dict[curr_step])
                except (ValueError, ReactionError) as e:
                        raise SimulationError(f"The following {type(e).__name__} was caught during simulating step {step.name}:\n{e}")
                if curr_step.type != BioReaction.INPUT:
                    self.graph.node(curr_step.name, fillcolor="#90EE90", style="filled")
                for step in self.graph_dict[curr_step]:
                    if step not in reaction_queue:
                        reaction_queue.append(step)  # step is child, so we can just check/extend their list
                    if step_dict[step]:  # Add to the input
                        step_dict[step] = self._format_input(step, output, step_dict[step])
                    else:  # Create the input
                        step_dict[step] = self._format_input(step, output)  # Check reaction type to generate correct input

        return output

    def _format_input(
        self,
        step: BioReactionStep,
        output: DNA | list[DNA],
        kwargs: dict[BioReactionStep, dict[str, Any]] | None = None
    ) -> dict[str, Any]:
        """Formats input for step for proper simulation"""
        if kwargs:  # We need to extend the input arg for the reaction
            match step.type:
                case BioReaction.AMPLIFY:
                    raise ValueError("Cannot provide more than one input for amplification reaction!")  # TODO: Update this
                case BioReaction.ANNEAL:
                    raise ValueError("Cannot provide more than one input for amplification reaction!")  # TODO: Update this
                case BioReaction.DIGEST | BioReaction.KLD | \
                BioReaction.INPUT:
                    kwargs["input"].append(output)
                case BioReaction.GIBSON | BioReaction.GG | \
                BioReaction.LIGATE:
                    kwargs["inputs"].append(output)
            return kwargs
        match step.type:
            case BioReaction.AMPLIFY:
                return {
                    "f_primer":         step.kwargs["f_primer"],
                    "r_primer":         step.kwargs["r_primer"],
                    "template":         output,
                    "min_binding_len":  step.kwargs["min_binding_length"],
                    "gel_extraction":   step.kwargs["gel_extraction"]
                }
            case BioReaction.ANNEAL:
                return {
                    "f_oligo":          step.kwargs["f_oligo"],
                    "r_oligo":          step.kwargs["r_oligo"]
                }
            case BioReaction.DIGEST:
                return {
                    "input":            output,
                    "enzymes":          step.kwargs["enzymes"],
                    "gel_extraction":   step.kwargs["gel_extraction"]
                }
            case BioReaction.GIBSON:
                return {
                    "inputs":           [output],
                    "min_homology_len": step.kwargs["min_homology_len"],
                    "max_homology_len": step.kwargs["max_homology_len"],
                    "gel_extraction":   step.kwargs["gel_extraction"]
                }
            case BioReaction.GG:
                return {
                    "inputs":           [output],
                    "enzyme":           step.kwargs["enzyme"],
                    "gel_extraction":   step.kwargs["gel_extraction"]
                }
            case BioReaction.KLD:
                return {
                    "input":            output
                }
            case BioReaction.LIGATE:
                return {
                    "inputs":           [output],
                    "gel_extraction":   step.kwargs["gel_extraction"]
                }

    def _generate_graph(
        self
    ) -> tuple[Digraph, dict[BioReactionGraph, set[BioReactionStep]], nx.DiGraph, dict[BioReactionStep, dict[str, Any]]]:
        """Generates a graph based on the input steps"""

        g: Digraph = Digraph(format="svg")
        g.attr(label=self.name, labelloc="t", fontsize="15")

        graph_dict: dict[BioReactionGraph, set[BioReactionStep]] = defaultdict(set)
        step_dict: dict[BioReactionStep, dict[str, Any]] = {}

        for step in self.edges:  # Create graphviz nodes and edges
            if step[0].type == BioReaction.INPUT:
                g.node(step[0].name, fillcolor="#AFD6F9", style="filled")
            else:
                g.node(step[0].name, fillcolor="#FFCCCB", style="filled")
            if step[1].type == BioReaction.INPUT:
                g.node(step[1].name, fillcolor="#AFD6F9", style="filled")
            else:
                g.node(step[1].name, fillcolor="#FFCCCB", style="filled")
            g.edge(step[0].name, step[1].name)
            graph_dict[step[0]].add(step[1])

            if step[0] not in step_dict:  # TODO: Need to make step_dict every time we run simulate otherwise the dict saves stuff and adds more each time
                step_dict[step[0]] = {}
            if step[1] not in step_dict:
                step_dict[step[1]] = {}

        # Verify graph is a DAG
        nx_graph: nx.DiGraph = nx.nx_pydot.from_pydot(pydot.graph_from_dot_data(g.source)[0])
        if not nx.is_directed_acyclic_graph(nx_graph):
            raise InvalidInstantiationException(f"{self.name} is not a directed acyclic graph.\n" \
                                                "Please ensure no products eventually feed into their own reactions.")

        return g, graph_dict, nx_graph, step_dict

    def remove_steps(self, *steps: BioReactionStep) -> None:
        """Removes one or more steps from the BioReactionGraph

        Parameters
        ----------
        steps: BioReactionStep
            One or more steps to remove from the graph

        Returns
        -------
        None
        """
        for step in steps:
            if not step in self.steps:
                raise ValueError(f"Step {step.name} not found!")
            for edge in self.edges:
                if step in edge:
                    self.edges.remove(edge)
        self.graph, self.graph_dict, self.nx_graph, self.step_dict = self._generate_graph()

    def remove_edges(self, *edges: tuple[BioReactionStep, BioReactionStep]) -> None:
        """Removes one or more steps from the BioReactionGraph

        Parameters
        ----------
        edges: BioReactionStep
            One or more edges to remove from the graph

        Returns
        -------
        None
        """
        for edge in edges:
            if not edge in self.edges:
                raise ValueError(f"Edge {edge} not found!")
            self.edges.remove(edge)
        self.graph, self.graph_dict, self.nx_graph, self.step_dict = self._generate_graph()

    def add_edges(self, *edges: tuple[BioReactionStep, BioReactionStep]) -> None:
        """Add one or more edges to the graph

        Parameters
        edges: BioReactionStep
            One or more edges to add to the graph

        Returns
        -------
        None
        """
        if all(edge in self.edges for edge in edges): return

        for edge in edges:
            if not edge in self.edges:
                self.edges.append(edge)

        self.graph, self.graph_dict, self.nx_graph, self.step_dict = self._generate_graph()
