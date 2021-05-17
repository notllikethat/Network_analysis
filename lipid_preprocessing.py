import pandas as pd
import numpy as np
import csv
import re
import requests
import time
from lipyd.lipproc import *
from lipyd import name
from lipyd import moldb


def AddAdducts(data, adducts):
    data = pd.merge(data, adducts, how="left")
    return data


def SeparatePosNegModes(data):
    positive_data = data[data.Mode == "pos"]
    negative_data = data[data.Mode == "neg"]

    positive_data = positive_data.reset_index(drop=True)
    negative_data = negative_data.reset_index(drop=True)

    positive_data = positive_data[["Lipid_ID", "MZ"]]
    negative_data = negative_data[["Lipid_ID", "MZ"]]

    return positive_data, negative_data


def GetAnnotation(pos_neg_data, db, mode, adducts=None):

    if adducts is not None:
        pos_neg_data = AddAdducts(pos_neg_data, adducts)
    annotation = pd.DataFrame(columns=["Lipid_ID", "SwissLipids_ID",
                                       "Formula",
                                       "Modification", "MZ"])

    for lipid_id in pos_neg_data.Lipid_ID:
        mass = float(pos_neg_data[pos_neg_data.Lipid_ID == lipid_id].MZ)
        result = db.adduct_lookup(mass, ionmode=mode)

        for modification in result.keys():
            r = result[modification]

            for Lipid_Record in r[1]:
                if adducts is not None:
                    # Checking adduct
                    tmp = pos_neg_data[pos_neg_data.Lipid_ID == lipid_id]
                    tmp = tmp.reset_index(drop=True)

                    if str(tmp.Adduct.iloc[0]) in modification:
                        if Lipid_Record.lab.db == "SwissLipids":
                            annotation = annotation.append({"Lipid_ID": lipid_id,
                                                            "SwissLipids_ID": Lipid_Record.lab.db_id,
                                                            "Formula": Lipid_Record.lab.formula,
                                                            "Modification": modification,
                                                            "MZ": mass},
                                                           ignore_index=True)
                else:
                    if Lipid_Record.lab.db == "SwissLipids":
                        annotation = annotation.append({"Lipid_ID": lipid_id,
                                                        "SwissLipids_ID": Lipid_Record.lab.db_id,
                                                        "Formula": Lipid_Record.lab.formula,
                                                        "Modification": modification,
                                                        "MZ": mass},
                                                       ignore_index=True)

    return annotation


def AnnotateDataWithLipyd(data, adducts=None):
    positive_data, negative_data = SeparatePosNegModes(data)

    print("...Composing MoleculeDatabaseAggregator....")
    db = moldb.MoleculeDatabaseAggregator(resources={
        'SwissLipids': (moldb.SwissLipids, {'levels': {'Species', 'Isomeric subspecies',
                                                       'Molecular subspecies',
                                                       'Structural subspecies'}}),
        'LipidMaps': (moldb.LipidMaps, {})
    })

    print("...Annotating positive lipids...")
    positive_data = GetAnnotation(positive_data, db=db, mode='pos', adducts=adducts)

    print("...Annotating negative lipids...")
    negative_data = GetAnnotation(negative_data, db=db, mode='neg', adducts=adducts)

    annotated_data = pd.concat([positive_data, negative_data])

    return annotated_data
