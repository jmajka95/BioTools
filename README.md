# BioTools
BioTools is a diverse, lightweight computational toolkit for biology, spanning sequence file conversion and parsing, virtual cloning, and more.\
\
**NOTE:** This is **not** a vibe-coded repository.

What can you do with BioTools?
1. Represent DNA, RNA, and Proteins
2. Parse FASTA and Genbank files
3. Virtually simulate reactions and high-throughput cloning workflows

Examples:
1. Visualization
```python
from biotools.dna import DNA
from biotools.rna import RNA
from biotools.protein import Protein

# This is how you create DNA, RNA, and Protein
dna =       DNA(seq="ATGC", name="foo")
rna =       RNA(seq="AUGC", name="bar")
protein =   Protein(seq="HEYYY", name="cześć")

dna.sequence() # To visualize the sequence, 
               # also works with RNA and Protein
```

2. File Parsing
```python
from biotools.bio_io import BioFileParser as BFP

genbank =   BFP().parse_genbank("path/to/genbank.gb")
# Multiple sequences can be parsed in one Genbank file
# Circularity is also inferred
fasta =     BFP().parse_fasta("path/to/fasta.fa", circular=False)
# Multiple sequences can be read in a file
# Linearity is assumed unless otherwise stated
```

3. Reaction Simulation
```python
from biotools.bio_reactions import digest, ligate

plasmid_1 = BFP().parse_genbank("path/to/genbank")
plasmid_2 = BFP().parse_genbank("path/to/genbank")
digest_1 =  digest(plasmid_1, ["HindIII", "XbaI"], (500, 550))
digest_2 =  digest(plasmid_2, ["HindIII", "XbaI"], (1500, 1600))
# We can digest both plasmids with HindIII and XbaI, and then
# extract products between 500-550 bp and 1500-1600 bp, inclusive

ligation = ligate(digest_1, digest_2)
# Cloning has been simulated and confirmed to work if the reaction
# returns products before you order expensive parts
```

3. Pooled Cloning
```python
from biotools.bio_reaction_step import BioReactionStep as BRS
from biotools.bio_annotation import Block
from biotools.bio_enums import BioOrientation, BioReaction
from biotools.bio_pool import BioPool

plasmid =   BFP().parse_genbank("path/to/genbank")
bp =        BioPool(BFP().parse_fasta("path/to/library/seqs"))
b =         Block((100, 150), "block", BioOrientation.FORWARD, bp)
# Parse DNA sequences into a pool, then make an annotation
plasmid.add_annotations(b)
# Adds the Block from slice 100 to 150 of plasmid

digest_reaction = BRS(
    BioReaction.DIGEST,
    "digest",
    {"input": plasmid, "enzymes": ["HindIII"], "gel_extraction": None}
)
digest_reaction.simulate()
# Simulates a digest reaction, generating outputs for all
# possible combinations of plasmid sequence from the BioPool
```

4. Multi-Step Simulations
```python
from biotools.oligo import Oligo
from biotools.bio_reaction_graph import BioReactionGraph as BRG
from biotools.bio_reactions import pcr

f_primer = plasmid.to_oligo((50, 75))
r_primer = plasmid.to_oligo((175, 200), reverse=True)
# Generate a primers from plasmid to use in a PCR reaction

pcr_step = BRS(
    BioReaction.AMPLIFY,
    "pcr_step",
    {
        "f_primer": f_primer,
        "r_primer": r_primer,
        "template": plasmid,
        "min_binding_length": 15, 
        "gel_extraction": None
    }
)
# We make a PCR step to start the graph

digest_step = BRS(
    BioReaction.DIGEST,
    "digest_step",
    {"input": None, "enzymes": ["BsaI"], "gel_extraction": None}
)
# This digest step happens after PCR. "input" is None because the 
# outputs from the PCR step will become the new inputs

reaction_graph = BRG((pcr_step, digest_step), name="test_simulation")
# Each tuple is an edge. Multiple edges can be added
reaction_graph.show_graph()
# This will produce a graphviz representation of the graph

reaction_graph.simulate()
# This simualates each step, passing products from previous steps to the next ones
# Final outputs will be returned if the simulation was a success

reaction_graph.show_graph()
# Completed steps become green, while uncompleted ones remain red
```

Why use BioTools?
1. BioTools was created with computational biologists in mind. It aims to be as lightweight, flexible, and simple as possible.
2. BioTools is designed to simulate complex, high-throughput cloning. As DNA synthesis costs continue to fall,
cheap, large-scale libraries are more attainable than ever. BioTools was built with this style of cloning in mind.
3. BioTools contains comprehensive tutorials in Python notebooks directly within the repository. You should never
have to dig through complicated, boring manuals.
