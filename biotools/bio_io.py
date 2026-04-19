from dna import DNA
from bio_annotation import BioAnnotation
from bio_enums import BioOrientation, BioProperty
import re

# TODO: Add visualization of stuff?
# TODO: SQLite database locally saves stuff?
# TODO: Would be cool to have like Fasta object that save important information saved in .fasta files

class BioFileParser():
    """Class for biological file parsing."""
    
    # TODO: Init needed for any reason?
    # If ^ isn't needed, then why make a class?
    # There has to be something to add in init...
    # def __init__(
    #     self,
    # ):
    #     pass

    def parse_fasta(self, file_path: str) -> list[DNA]:
        """Parses a .fasta file, saving all sequences as DNA objects.
        Assumes sequences are double-stranded DNA."""
        
        # TODO: This assumes that they will be DNA sequences. Any way to infer this? Could also have parameter
        # TODO: Add FileParsingError exception and add as appropriate

        name: str = ""
        seq: str
        dna_list: list[DNA] = []
        save: bool = False

        with open(file_path, "r") as file:
            for line in file:
                if line.startswith(">"):
                    seq = "" # Reset sequence
                    m = re.match(r"^>(.*)", line.rstrip("\n"))
                    if m:
                        name = m.group(1)
                    else:
                        name = ""
                    save = True
                elif save:
                    addition = line.rstrip("\n")
                    if (addition != ""):
                        seq += addition
                    else:
                        dna_list.append(DNA(seq, name))                
                        save = False
            if save:
                dna_list.append(DNA(seq, name))  
                
        return dna_list
    
    def parse_genbank(self, file_path: str) -> list[DNA]:
        """Parses a .gbk file, saving the sequence as a DNA object."""
        
        sequences: list[DNA] = []
        annotations: list[BioAnnotation] = []
        seq: str = ""
        name: str = ""
        annot_name: str = ""
        add_annotations: bool = False
        add_sequence: bool = False
        orientation: BioOrientation = BioOrientation.FORWARD
        circular: BioProperty = BioProperty.LINEAR
        span: tuple[int, int] = (0, 0)

        # Multiple sequences can be separated by //
        with open(file_path, "r") as file:
            for line in file:
                if line.startswith("//"):
                    sequences.append(DNA(seq, name, circular=circular, annotations=annotations))
                    seq = ""
                    name = ""
                    annotations = []
                    add_sequence = False
                    circular = BioProperty.LINEAR
                else:
                    # Find name
                    if line.startswith("LOCUS"):
                        m = re.match(r"^LOCUS\s+([A-Za-z0-9_]+)\s+", line)
                        if m:
                            name = m.group(1)
                        if "circular" in line:
                            circular = BioProperty.CIRCULAR
                    elif line.startswith("FEATURES"): # Annotations creation workflow
                        add_annotations = True
                    elif line.startswith("ORIGIN"):
                        add_annotations = False
                        add_sequence = True
                    elif add_annotations:
                        if "complement" in line:
                            m = re.search(r"complement\((\d+)\.\.(\d+)\)", line)
                            if m:
                                span = int(m.group(1)), int(m.group(2))
                            orientation = BioOrientation.REVERSE
                        elif ".." in line:
                            m = re.search(r"\s*(\d+)\.\.(\d+)", line)
                            if m:
                                span = int(m.group(1)), int(m.group(2))
                        if "label" in line:
                            m = re.search(r"\s*/label=\"(.*)\"", line)
                            if m:
                                annot_name = m.group(1)
                            
                            annotations.append(BioAnnotation(span, annot_name, orientation))
                            annot_name = ""
                            span = (0,0)
                            orientation = BioOrientation.FORWARD
                    elif add_sequence:
                        m = re.search(r"\d+\s+([atgcnwsmkrybvhdATGCNWSMKRYBVHD\s]+)", line) # TODO: Allow degenerate bases?
                        new_seq = re.sub(r"\s+", "", m.group(1))
                        seq += new_seq
        return sequences

    def to_fasta(self, *seqs: DNA | list[DNA], filename: str = "test.fasta") -> None:
        """Generates a fasta file of the provided sequences"""
        if not filename.endswith(".fasta"):
            filename += ".fasta"
        with open(filename, "w") as f:
            for seq in seqs:
                f.write(f">{seq.name}\n")
                f.write(f"{seq.seq}\n")
                f.write("\n")
