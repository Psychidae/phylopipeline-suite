import os

# --- PATH CONSTANTS ---
# Assuming structure: root/modules/constants.py
MODULES_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(MODULES_DIR)
TOOLS_DIR = os.path.join(ROOT_DIR, "tools")

IUPAC_CODES = {
    'A': 'Adenine', 'C': 'Cytosine', 'G': 'Guanine', 'T': 'Thymine', '-': 'Gap',
    'R': 'A or G', 'Y': 'C or T', 'W': 'A or T', 'S': 'G or C',
    'K': 'G or T', 'M': 'A or C', 'B': 'not A', 'D': 'not C',
    'H': 'not G', 'V': 'not T', 'N': 'Any', 'X': 'Any'
}

GENETIC_CODES = {
    "Standard (1)": 1,
    "Vertebrate Mitochondrial (2)": 2,
    "Yeast Mitochondrial (3)": 3,
    "Mold/Protozoan Mitochondrial (4)": 4,
    "Invertebrate Mitochondrial (5)": 5,
    "Ciliate/Dasycladacean/Hexamita Nuclear (6)": 6,
    "Echinoderm/Flatworm Mitochondrial (9)": 9,
    "Euplotid Nuclear (10)": 10,
    "Bacterial/Archaeal/Plant Plastid (11)": 11,
}
