from Bio import Entrez

Entrez.email = "dkwabiad@mit.edu"
Entrez.tool = "tRNA_Synthetase_Export.py"

MAX_RETURNED = 50 #number of organisms 
PHYLOGENETIC_CLOSENESS = 4 #For tRNAdb1: From 1-4, 1 being close and 4 being really far
OUTPUT_DIRECTORY = r"C:\Users\Bernard Kwabi-Addo\Documents\tRNA-aaRS Pipeline\trna-aars-pipeline\src\output" + "/"


AA_SET = {("Alanine", "Alanyl",), ("Arginine", "Arginyl",), ("Asparagine", "Asparaginyl",), ("Aspartic acid", "Aspartate", "Aspartyl",), 
                       ("Cysteine", "Cysteinyl",), ("Glutamic acid", "Glutamate", "Glutamyl",), ("Glutamine", "Glutaminyl",), ("Glycine", "Glycyl",), 
                       ("Histidine", "Histidyl",), ("Isoleucine", "Isoleucyl",), ("Leucine", "Leucyl", ), ("Lysine", "Lysyl",), ("Methionine", "Methionyl", ), 
                       ("Phenylalanine", "Phenylalanyl",), ("Proline", "Prolyl",), ("Serine", "Seryl",), ("Threonine", "Threonyl",), ("Tryptophan", "Tryptophanyl",), 
                       ("Tyrosine", "Tyrosyl", ), ("Valine", "Valyl",), ("Phosphoserine", "Phosphoseryl",), ("Pyrrolysine", "Pyrrolysyl",)}


AA_DICT = {'Ala': ('Alanine', 'Alanyl'), 'Arg': ('Arginine', 'Arginyl'), 
                    'Asn': ('Asparagine', 'Asparaginyl'), 'Asp': ('Aspartate', 'Aspartyl'), 
                    'Cys': ('Cysteine', 'Cysteinyl'), 'Glu': ('Glutamate', 'Glutamyl'), 
                    'Gln': ('Glutamine', 'Glutaminyl'), 'Gly': ('Glycine', 'Glycyl'), 'His': ('Histidine', 'Histidyl'), 
                    'Ile': ('Isoleucine', 'Isoleucyl'), 'Leu': ('Leucine', 'Leucyl'), 'Lys': ('Lysine', 'Lysyl'), 
                    'Met': ('Methionine', 'Methionyl'), 'Phe': ('Phenylalanine', 'Phenylalanyl'), 'Pro': ('Proline', 'Prolyl'), 
                    'Ser': ('Serine', 'Seryl'), 'Thr': ('Threonine', 'Threonyl'), 'Trp': ('Tryptophan', 'Tryptophanyl'), 
                    'Tyr': ('Tyrosine', 'Tyrosyl'), 'Val': ('Valine', 'Valyl'), 'Sep': ('Phosphoserine', 'Phosphoseryl'), 'Pyl': ('Pyrrolysine', 'Pyrrolysyl')}

TRNADB1_ALL_SEQ = r"datasets\full_tRNAdb.fst" #Replace with your path name to All Sequences
ECOLI_K12_PATH = r"datasets\tRNAdb1_E_Coli_K12.fst" #From tRNA Database 1

## From tRNAdb1
ENTERO_PATH = r"datasets\tRNAdb1_Enterobacteriales.fst"
GAMMA_PATH = r"datasets\tRNAdb1_Gammaproteobacteria.fst"
PROTEO_PATH = r"datasets\tRNAdb1_Proteobacteria.fst"
BACTERIA_PATH = r"datasets\tRNAdb1_Bacteria.fst"
ARCHEA_PATH = r"datasets\tRNAdb1_Archea.fst"
EUKARYA_PATH = r"datasets\tRNAdb1_Eukaryota.fst"
VIRUSES_PATH = r"datasets\tRNAdb1_Viruses.fst"


## From tRNAdb2
VIRUSES_PATH_2 = r"datasets\trna_sequence_virus_1.fasta"
PLANT_PATH_2 = r"datasets\trna_sequence_plant_1.fasta"
PHAGE_PATH_2 = r"datasets\trna_sequence_phage_1.fasta"
BACTERIA_PATH_2 = r"datasets\trna_sequence_cmp_bac_1.fasta"
ARCHEA_PATH_2 = r"datasets\trna_sequence_cmp_arc_1.fasta"
FUNGI_PATH_2 = r"datasets\trna_sequence_fungi_1.fasta"
PLASMID_PATH_2 = r"datasets\trna_sequence_plasmid_1.fasta"
ENV_PATH_2 = r"datasets\trna_sequence_env_1.fasta"
CHLORO_PATH_2 = r"datasets\trna_sequence_chloro_1.fasta"






