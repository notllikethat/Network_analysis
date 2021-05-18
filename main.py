import lipid_preprocessing as lp
import pandas as pd

# Reading data
TL_4sp = pd.read_csv('data/newTL_4sp.txt', sep='\t')

# Renaming columns
TL_4sp.columns = ["Lipid_ID", "Lipid_class", "Bulk_structure",
                  "Detailed_structure", "Adduct", "Index_Othermode", "MZ", "RT"]

# Replacing adducts with supported form
TL_4sp.loc[TL_4sp["Adduct"] == "M-H+FA", "Adduct"] = "M+HCOO"

# Adding info on modes of the features
pos_neg_data = pd.read_csv('data/pos_neg_mode.csv', sep='\t')
TL_4sp = pd.merge(TL_4sp, pos_neg_data, how="left")

''' Input should be pd.DataFrame, formatted like this:
    --- Lipid_ID, MZ, Mode ---
    where mode should be either "pos" or "neg"
    
    If there are adducts, add the parameter "adducts".
    "adducts" should be pd.DataFrame, formatted like this:
    --- Lipid_ID, Adduct ---
'''

# Adjusting input
data = TL_4sp[["Lipid_ID", "MZ", "Mode"]]
adducts = TL_4sp[["Lipid_ID", "Adduct"]]

# # Annotating data with LIPYD
# annotated_data = lp.AnnotateDataWithLipyd(data=data, adducts=adducts)
# print(annotated_data.head())
# print(len(annotated_data))
#
#
# # Saving intermediate results
# annotated_data.to_csv(r'results/lipyd_annotated_data.csv', index=False)

'''
    If you want to check generated annotation, you should provide parameter "curation"
    "curation" should be pd.DataFrame, formatted like this:
    --- Lipid_ID, Curation ---
    And look similar to this: SM 38:1 (LipidMaps style)
    Otherwise the function will return all possible annotations for lipids
'''

# # Reading intermediate results
# annotated_data = pd.read_csv('results/lipyd_annotated_data.csv', sep=',')
# # print(annotated_data.head())
#
#
# # Adjusting input
# curation = other[["Lipid_ID", "Bulk_structure"]]
# curation.columns = ["Lipid_ID", "Curation"]


# # For test use only
# annotated_data = annotated_data[0:10]
#
# # Getting representatives
# representatives, repr_to_swl = lp.GetRepresentatives(data=annotated_data, curation=curation, check_annotation=True)
# print(repr_to_swl.head())
#
# # Saving intermediate results
# representatives.to_csv(r'results/representatives.csv', index=False)
# repr_to_swl.to_csv(r'results/representatives_to_swl.csv', index=False)


'''
    Parameter "representatives" should be pd.DataFrame, formatted like this:
    --- Lipid_ID, SwissLipids_ID ---
    Parameter "repr_to_swl" should be pd.DataFrame, formatted like this:
    --- SwissLipids_ID, Representative_ID --- 
'''

# # Reading intermediate results
# representatives = pd.read_csv('results/representatives.csv', sep=',')
# repr_to_swl = pd.read_csv('results/representatives_to_swl.csv', sep=',')
# other = pd.read_csv('results/lipyd_other.csv', sep=',')
# # print(representatives.head())
#
# # Adjusting input
# representatives = representatives[["Lipid_ID", "SwissLipids_ID"]]
#
# # Go back from representatives to lipids
# annotated_data = lp.RepresentativesToLipids(representatives, repr_to_swl)
# # print(annotated_data.head())
#
#
# # Saving intermediate results
# annotated_data.to_csv(r'results/annotated_data_representatives.csv', index=False)


'''
    Parameter "annotated_data" should be pd.DataFrame, formatted like this:
    --- Lipid_ID, SwissLipids_ID , Representative_ID ---
'''

# # Reading intermediate results
# annotated_data = pd.read_csv('results/annotated_data_representatives.csv', sep=',')
#
#
# # Mapping to graph
# annotated_data = lp.MappingToGraph(annotated_data)
# annotated_data.head()
#
# # Saving intermediate results
# annotated_data.to_csv(r'results/annotated_data_chebi.csv', index=False)


'''
    Merging everything back
'''

# Reading intermediate results
annotated_data = pd.read_csv('results/annotated_data_chebi.csv', sep=',')

# Reading data
TL_4sp = pd.read_csv('data/newTL_4sp.txt', sep='\t')

# Renaming columns
TL_4sp.columns = ["Lipid_ID", "Lipid_class", "Bulk_structure",
                  "Detailed_structure", "Adduct", "Index_Othermode", "MZ", "RT"]


# Merging everything back
result = pd.merge(annotated_data, TL_4sp, how="left")

result = result[["Lipid_ID", "ChEBI_ID", "SwissLipids_ID", "Level",
                 "Initial_SwissLipids_ID", "Depth", "Representative_ID",
                 "Lipid_class", "Bulk_structure",
                 "Detailed_structure", "Adduct", "Index_Othermode", "MZ", "RT"]]

result.columns = ['index_newtable_4sp', 'Lipid class', 'Bulk structure',
                  'Detailed structure', 'adduct', 'index_othermode_newtable_4sp',
                  'newtable_4sp_mz', 'newtable_4sp_rt',
                  "ChEBI_ID", "SwissLipids_ID", "Level",
                  "Initial_SwissLipids_ID", "Depth", "Representative_ID"]


# Saving final results
annotated_data.to_csv(r'results/annotated_data_chebi.csv', index=False)
