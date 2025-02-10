from constants.global_constants import *
from common_functions import *

#Define output file names
TRNA_DATASET_NAME = "fungi_trna_dataset.csv"
ALIGNMENT_OUTPUT_NAME = "alignment_values_fungi.csv"
AARS_FILE_NAME = "aars_output.csv"
AARS_TRNA_FILE_NAME = "associated_aars_trna.csv"

no_search_results = []

## Get tRNAs
tRNA_df = get_tRNAdb_data_2(FUNGI_PATH_2) # change path name to path of interest
tRNA_df["Organism Name"] = tRNA_df["Species"].apply(lambda s: s[:s.find('(') -1] if s.find('(') != -1 else s ) #Gets rid of the () parts of species (not incldued in ncbi database)

full_trna_df = tRNA_df
full_trna_df["Valid tRNA Sequence"] = full_trna_df["tRNA Sequence"].apply(is_valid_dna_sequence)
tRNA_df = tRNA_df.loc[full_trna_df["Valid tRNA Sequence"]]
tRNA_df = tRNA_df.drop(columns=["Valid tRNA Sequence"])

#Save tRNA data as a CSV
tRNA_df.to_csv(OUTPUT_DIRECTORY + TRNA_DATASET_NAME, index=False)

#Conduct alignment scores
get_tRNA_type(tRNA_df) #Makes a new column in the dataframe that is a combination of the amino acid and anticodon, e.g. Ala(GGC)
tRNA_df["Sequence"] = tRNA_df["tRNA Sequence"]

k12_df = get_tRNAdb_data(ECOLI_K12_PATH)
get_tRNA_type(k12_df)

k12_df["tRNA Sequence"] = k12_df["tRNA Sequence"].apply(lambda x: remove_position_chars(x)) #remove the positional structure for the trna sequneces

alignment_df = align_dataframes(tRNA_df, k12_df)
alignment_output_df = alignment_df.sort_values(by=['Amino Acid Name', 'AlignmentScore'], ascending=True)

alignment_output_df = alignment_output_df.drop_duplicates()
alignment_output_df.to_csv(OUTPUT_DIRECTORY + ALIGNMENT_OUTPUT_NAME, index=False)

# Search for synthetases that correspond to the DNA 
has_tax_id = (alignment_output_df["Taxonomic ID"] != "") | (alignment_output_df["Taxonomic ID"] != None)


org_name_set = set(alignment_output_df[has_tax_id]["Organism Name"])
total = len(org_name_set)
print(total)

total_data_list = []
# org_name_set = {"Escherichia_coli_K12"}
# org_name_set = {"Candidatus Nanopusillus acidilobi 7A"}
# org_name_set = {"Vulcanisaeta moutnovskia 768-28"}
count = 0

for name in org_name_set:
    count += 1
    txid = species_tax_id[name]
    
    print(f"Checking Organism {name}, organism {count} of {total} and txid {txid}")
    
    
    if ('(' in txid): # this is a tuple of values, meaning that it specifies a taxonomy id and strain
        txid = convert_string_to_tuple(txid)
        strain = txid[1]
        search_term = f"txid{txid[0]}[organism:noexp] AND {strain}[strain] AND (trna ligase OR trna synthetase)"
        # print(search_term)
    else:
        search_term = f"txid{txid}[organism:noexp] AND (trna ligase OR trna synthetase)"
        
    # if isinstance(species_tax_id[name], tuple): #specifies the strain as well
    #     strain = species_tax_id[name][1]
    #     org_name = name[:-len(strain)]
    #     search_term = f"{name[0]}[organism] and {name[1]}[strain] and chloroplast and (trna ligase or trna synthetase)"
    # else: #just the taxon id 
    #     search_term = f"{name}[organism] and chloroplast and (trna ligase or trna synthetase)" #search query for tRNA database
    
    
    
    max_returned = MAX_RETURNED # maximum number of sequences returned
    
    total_data_list += get_extracted_aaRS_from_NCBI(name, search_term, max_returned, no_search_results)
    
    # try:
    #     total_data_list += get_extracted_aaRS_from_NCBI(name, search_term, max_returned)
    # except: 
    #     print(f'Organism {name} had an error occur')
    #     with_errors.append(name)
    #     continue

aaRS_df = pd.DataFrame(total_data_list)
print("Organisms with no search results:", no_search_results)



#Replace each amino acid with its three letter code to be consistent
aaRS_df["Amino Acid Name"].replace(AA_DICT, inplace=True)

#remove sequences that don't have 'tRNA ligase' or 'tRNA synthetase' in their title
df = aaRS_df
df["aaRS Description"] = df["aaRS Description"].str.lower()
mask = (df['aaRS Description'].str.contains("trna ligase") | df["aaRS Description"].str.contains("trna synthetase")) & ~(df['aaRS Description'].str.contains("partial"))  \
    &  ~(df['aaRS Description'].str.contains("chain")) &  ~(df['aaRS Description'].str.contains("domain")) & ~(df['aaRS Description'].str.contains("related")) \
    & ~(df['aaRS Description'].str.contains("family")) & ~(df['aaRS Description'].str.contains("fragment")) & ~(df['aaRS Description'].str.contains("peptide")) \
    & ~(df['aaRS Description'].str.contains("binding")) & ~(df['aaRS Description'].str.contains("core"))
    
filtered_df = df[mask]

# make sure the organisms that are in the result match the organisms that were searched
filtered_df["Organism Name"] = filtered_df["Organism Name"].str.replace('_', ' ')
# tRNA_df["Organism Name"] = tRNA_df["Organism Name"].str.replace('_', ' ')
# mask = filtered_df.apply(lambda row: row['Organism Name'].lower() in row['Actual Organism'].lower(), axis=1)
# mask = filtered_df.apply(lambda row: row['Actual Organism'].lower().startswith(row['Organism Name'].lower()), axis=1)
# mask = filtered_df.apply(lambda row: row['Organism Name'].lower() == row['Actual Organism'].lower(), axis=1)

# mask = mask.reset_index(drop=True)
final_df = filtered_df[mask]
final_df["Probable?"] = "probable" in final_df["aaRS Description"]
final_df["Putative?"] = "putative" in final_df["aaRS Description"]


#merging tRNA, aaRS 

alignment_output_df.rename(columns={'Amino Acid': 'Amino Acid Name'}, inplace=True)

output_df = associate_tRNA_aaRS_df(alignment_output_df, final_df)
l = output_df.columns

# output_df = output_df[["Organism Name", "Actual Organism", "aaRS Description", "AlignmentScore",  "Amino Acid Name", "Anticodon", "tRNA Sequence", "aaRS Sequence", "Probable?", "Putative?", "Origin", "Accession ID", "SequenceID", "E.Coli K12 Cognate tRNA", "E.Coli K12 Cognate SequenceID", "Taxonomic ID"]]
#Decided not to put tRNA name in the final output -- just a combo of Amino Acid name and Anticodon

output_df.drop_duplicates(inplace=True)

final_df.to_csv(OUTPUT_DIRECTORY + AARS_FILE_NAME, index=False)
output_df.to_csv(OUTPUT_DIRECTORY + AARS_TRNA_FILE_NAME, index=False)