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
other = TL_4sp[["Lipid_ID", "Lipid_class", "Bulk_structure",
                "Detailed_structure", "Index_Othermode", "RT"]]


# # Annotating data with LIPYD
# annotated_data = lp.AnnotateDataWithLipyd(data=data, adducts=adducts)
# print(annotated_data.head())
# print(len(annotated_data))
#
#
# # Saving intermediate results
# annotated_data.to_csv(r'results/lipyd_annotated_data.csv', index=False)
# other.to_csv(r'results/lipyd_other.csv', index=False)


# Reading intermediate results
annotated_data = pd.read_csv('results/lipyd_annotated_data.csv', sep=',')
other = pd.read_csv('results/lipyd_other.csv', sep=',')
print(annotated_data.head())

'''
    If you want to check generated annotation, you should provide parameter "curation"
    "curation" should be pd.DataFrame, formatted like this:
    --- Lipid_ID, Curation ---
    And look similar to this: SM 38:1 (LipidMaps style)
    Otherwise the function will return all possible annotations for lipids
'''

# Adjusting input
curation = other[["Lipid_ID", "Bulk_structure"]]

# Getting representatives
representatives = lp.GetRepresentatives(data=annotated_data, curation=curation, check_annotation=True)
