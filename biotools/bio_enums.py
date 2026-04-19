# Enums for common indentifiers
from enum import Enum

class BioMolecule(Enum):
    """Class for representing biomolecules."""
    DNA = "DNA"
    OLIGO = "Oligo"
    RNA = "RNA"
    PROTEIN = "Protein"

class BioProperty(Enum):
    """Class for representing biological properties like strandedness, circularity, etc."""
    CIRCULAR = "Circular"
    LINEAR = "Linear"
    DOUBLE_STRANDED = "Double-Stranded"
    SINGLE_STRANDED = "Single-Stranded"
    CUT = "Cut" # For resembling products of reactions

class BioOrientation(Enum):
    """Class for representing orientation properties."""
    FORWARD = "Forward"
    REVERSE = "Reverse"

class BioFile(Enum):
    """Class for biological file names."""
    FASTA = "Fasta"
    GENBANK = "Genbank"
    FASTQ = "Fastq"
