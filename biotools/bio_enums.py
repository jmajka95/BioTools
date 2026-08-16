"""Enums for common biological identifying information. 
"""

from enum import Enum

class BioMolecule(Enum):
    """Class for representing biomolecules."""
    DNA =               "DNA"
    OLIGO =             "Oligo"
    RNA =               "RNA"
    PROTEIN =           "Protein"

class BioProperty(Enum):
    """Class for representing biological properties like strandedness, circularity, etc."""
    CIRCULAR =          "Circular"
    LINEAR =            "Linear"
    DOUBLE_STRANDED =   "Double-Stranded"
    SINGLE_STRANDED =   "Single-Stranded"
    CUT =               "Cut"  # For resembling products of reactions

class BioOrientation(Enum):
    """Class for representing orientation properties."""
    FORWARD =           "Forward"
    REVERSE =           "Reverse"
    TOP =               "Top"
    BOTTOM =            "Bottom"

class BioFile(Enum):
    """Class for biological file names."""
    FASTA =             "Fasta"
    GENBANK =           "Genbank"
    FASTQ =             "Fastq"

class BioReaction(Enum):
    """Class representing cloning reactions used in BioReactionGraphs."""
    AMPLIFY =           "Amplify"
    ANNEAL =            "Anneal"
    DIGEST =            "Digest"
    GIBSON =            "Gibson"
    GG =                "Golden Gate"
    INPUT =             "Input"
    KLD =               "Kld"
    LIGATE =            "Ligate"
