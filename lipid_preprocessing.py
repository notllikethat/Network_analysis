import matplotlib
import pandas as pd
import numpy as np
import csv
import re
import requests
from lipyd.lipproc import *
from lipyd import name
from lipyd import moldb
from networkx import *


''' 
    This part is dedicated to creating the annotation of the data via lipyd
'''

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


def AddAdducts(data, adducts):
    data = pd.merge(data, adducts, how="left")
    return data


''' 
    This part is dedicated to getting representatives
'''


def GetRepresentatives(data, curation=None, check_annotation=False):
    levels_data = pd.read_csv("data/levels_data.csv", sep=",")
    levels_data.columns = ["SwissLipids_ID", "Level"]

    # Adding levels to data
    swl_ids = AddingLevelsToData(data, levels_data=levels_data)

    # Creating graph
    g = CreateAnnotatedGraph(levels_data=levels_data)

    # Getting representatives
    representatives, more_species, no_species, getting_error = RunDFS(swl_ids, graph=g)

    # Compress results:
    if len(getting_error) > 0:
        print("Something went wrong. We got KeyError")

    elif len(more_species) > 0:
        print("Something went wrong. We got too many species")

    elif len(no_species) > 0:
        representatives_no_species = GetRepresentativesNoSpecies(no_species=no_species,
                                                                 graph=g,
                                                                 levels=levels_data)

    representatives = pd.concat([representatives, representatives_no_species])
    representatives = representatives.drop_duplicates()
    representatives = representatives.reset_index(drop=True)

    # Merging results with data
    representatives = pd.merge(representatives, data, how="left")
    representatives = representatives[["Representative_ID", "Lipid_ID"]]
    representatives = representatives.drop_duplicates()
    representatives = representatives.reset_index(drop=True)

    representatives = GetRepresentativesClasses(representatives)

    # Auto curation check
    if check_annotation and curation is not None:
        representatives = CheckAnnotation(representatives, curation)

    return representatives


def CheckAnnotation(representatives, curation):
    representatives = pd.merge(representatives, curation, how="left")

    pre_columns = list(representatives.columns)

    repl = representatives.drop_duplicates(subset=['Class_SwissLipids'], keep="first")
    repl['Abbreviation_SwissLipids'] = repl['Abbreviation_SwissLipids'].str.replace('\(d', '(')
    repl['Abbreviation_SwissLipids'] = repl['Abbreviation_SwissLipids'].str.replace(' \(O-', '_O(')
    repl['Abbreviation_SwissLipids'] = repl['Abbreviation_SwissLipids'].str.replace('\(O-', '_O(')
    repl['Abbreviation_SwissLipids'] = repl['Abbreviation_SwissLipids'].str.replace(' \(P-', '_P(')
    repl['Abbreviation_SwissLipids'] = repl['Abbreviation_SwissLipids'].str.replace('\(P-', '_P(')

    tmp = repl["Abbreviation_SwissLipids"].str.split("(", n=1, expand=True)
    repl["Class_Abbr_SL"] = tmp[0]
    repl = repl[["Class_SwissLipids", "Class_Abbr_SL"]]

    anno = pd.read_csv("Databases_data/full_classes_data_LM_SWL.csv", sep=",")
    anno["Annotation"] = anno["Abbreviation_LipidMaps"]
    anno = anno[["SwissLipids_ID", "Annotation", "Class_SwissLipids", "Abbreviation_SwissLipids",
                 "LipidMaps_ID", "Class_LipidMaps", "Abbreviation_LipidMaps", "SwissLipids_name"]]
    anno = anno.reset_index(drop=True)

    lm_classes = list(anno.Class_LipidMaps.drop_duplicates().values)
    lm_classes.sort()
    if "-" in lm_classes:
        lm_classes = lm_classes[1:]

    cls_df = pd.DataFrame(columns=["Class_LipidMaps", "Class_Abbr_LM"])
    for cl in lm_classes:
        abb = re.findall(r'\[(.+)\]', str(cl))
        cls_df = cls_df.append({"Class_LipidMaps": cl,
                                "Class_Abbr_LM": abb[0]},
                               ignore_index=True)

    representatives = pd.merge(representatives, repl, how="left")
    representatives = pd.merge(representatives, cls_df, how="left")
    representatives.loc[representatives["Class_Abbr_SL"].isnull(), "Class_Abbr_SL"] = "-"
    representatives.loc[representatives["Class_Abbr_LM"].isnull(), "Class_Abbr_LM"] = "-"
    representatives.loc[representatives["Class_Abbr_SL"] == "DHDG", "Class_Abbr_SL"] = "FA"
    representatives.loc[representatives["Class_SwissLipids"] == "Fatty acid methyl esters",
                        "Class_Abbr_SL"] = "FA"

    representatives = representatives[pre_columns]

    representatives = representatives[['SwissLipids_ID', 'Lipid_ID', 'Annotation', 'Curation',
                                       'Class_SwissLipids', 'Abbreviation_SwissLipids', 'LipidMaps_ID',
                                       'Class_LipidMaps', 'Abbreviation_LipidMaps', 'SwissLipids_name']]

    representatives.loc[representatives["Annotation"] == representatives["SwissLipids_name"], "Annotation"] = "-"
    representatives['Annotation'] = representatives['Annotation'].str.replace(' \(O-', '_O(')
    representatives['Annotation'] = representatives['Annotation'].str.replace('\(O-', '_O(')
    representatives['Annotation'] = representatives['Annotation'].str.replace(' \(P-', '_P(')
    representatives['Annotation'] = representatives['Annotation'].str.replace('\(P-', '_P(')
    representatives['Annotation'] = representatives['Annotation'].str.replace(' \(', '(')

    representatives.loc[~representatives["Annotation"].str.contains("NAPE"),
                        "Annotation"] = representatives['Annotation'].str.replace('\)', '')
    representatives.loc[~representatives["Annotation"].str.contains("NAPE"),
                        "Annotation"] = representatives['Annotation'].str.replace('\(', ' ')

    representatives.loc[representatives["Annotation"] == "-",
                        "Annotation"] = representatives["SwissLipids_name"]

    # Real check
    representatives["OK"] = "-"
    representatives.loc[representatives["Annotation"] == representatives["Bulk structure"], "OK"] = "+"

    representatives = representatives[representatives.OK == "+"]

    representatives = representatives[['Lipid_ID', 'SwissLipids_ID']]
    return representatives


def GetRepresentativesClasses(representatives):
    annotation = pd.read_csv("data/full_classes_data_LM_SWL.csv", sep=",")
    annotation["Annotation"] = annotation["Abbreviation_LipidMaps"]
    annotation = annotation[["SwissLipids_ID", "Annotation", "Class_SwissLipids", "Abbreviation_SwissLipids",
                             "LipidMaps_ID", "Class_LipidMaps", "Abbreviation_LipidMaps", "SwissLipids_name"]]
    annotation = annotation[annotation.SwissLipids_ID != "-"]
    annotation = annotation.reset_index(drop=True)

    annotation.loc[annotation["Annotation"] == "-", "Annotation"] = annotation["Abbreviation_SwissLipids"]
    annotation.loc[annotation["Annotation"] == "-", "Annotation"] = annotation["SwissLipids_name"]
    annotation = annotation.drop_duplicates()
    annotation = annotation.reset_index(drop=True)

    representatives = pd.merge(representatives, annotation, how="left")

    return representatives


def GetRepresentativesNoSpecies(no_species, graph, levels):
    representatives_no_species = pd.DataFrame(columns=["SwissLipids_ID", "Representative_ID"])

    for lipid in no_species:
        dfs_swl_ids = GetParentsWithLevels(graph=graph, levels=levels, source=lipid)
        df = dfs_swl_ids[(dfs_swl_ids.Level != "Class") &
                         (dfs_swl_ids.Level != "Category") &
                         (dfs_swl_ids.Level != "-")]
        if len(df) == 1:
            df = df.reset_index(drop=True)
            repr_id = df.SwissLipids_ID.iloc[0]
            representatives_no_species = representatives_no_species.append({"SwissLipids_ID": lipid,
                                                                            "Representative_ID": repr_id},
                                                                           ignore_index=True)
        else:
            representatives_no_species = representatives_no_species.append({"SwissLipids_ID": lipid,
                                                                            "Representative_ID": lipid},
                                                                           ignore_index=True)
    return representatives_no_species


def GetParentsWithLevels(graph, levels, source):
    dfs_output = list(nx.dfs_preorder_nodes(graph, source=source))
    dfs_swl_ids = pd.DataFrame({"SwissLipids_ID": dfs_output})
    dfs_swl_ids = pd.merge(dfs_swl_ids, levels, how="left",
                           left_on="SwissLipids_ID", right_on="SwissLipids_ID")
    dfs_swl_ids.loc[dfs_swl_ids["Level"].isnull(), "Level"] = "-"
    return dfs_swl_ids


def RunDFS(swl_ids, graph):
    more_species = []
    no_species = []
    getting_error = []
    representatives = pd.DataFrame(columns=["SwissLipids_ID", "Representative_ID"])

    for i in range(len(swl_ids)):
        if swl_ids.SwissLipids_ID[i] in graph.nodes():
            try:
                SG = [n for n in nx.dfs_preorder_nodes(graph, source=swl_ids.SwissLipids_ID[i])
                      if graph.nodes[n]["Level"] == "Species"]

            except KeyError:
                getting_error.append(swl_ids.SwissLipids_ID[i])

            if len(SG) > 1:
                more_species.append(swl_ids.SwissLipids_ID[i])

            elif len(SG) < 1:
                no_species.append(swl_ids.SwissLipids_ID[i])

            elif len(SG) == 1:
                repr_id = SG[0]
                representatives = representatives.append({"SwissLipids_ID": swl_ids.SwissLipids_ID[i],
                                                          "Representative_ID": repr_id},
                                                         ignore_index=True)
        else:
            no_species.append(swl_ids.SwissLipids_ID[i])

    return representatives, more_species, no_species, getting_error


def CreateAnnotatedGraph(levels_data):
    file = "data/graph.txt"
    n, m, edges_list = GetData(file)

    # Creating graph
    print("\n..Creating Graph..")
    g = MakeGraph(n, edges_list)

    # Annotating graph
    g_ann = pd.DataFrame({"SwissLipids_ID": list(g.nodes())})
    levels_data = pd.merge(g_ann, levels_data, how="left",
                           left_on="SwissLipids_ID", right_on="SwissLipids_ID")
    levels_data.loc[levels_data["Level"].isnull(), "Level"] = "-"
    g = AnnotateGraph(g, levels=levels_data)

    return g


def AddingLevelsToData(data, levels_data):
    swl_ids = data.SwissLipids_ID.drop_duplicates().values
    swl_ids = pd.merge(swl_ids, levels_data, how="left",
                       left_on="SwissLipids_ID", right_on="SwissLipids_ID")
    swl_ids.loc[swl_ids["Level"].isnull(), "Level"] = "-"
    return swl_ids


def AnnotateGraph(g, levels):
    node_attr = levels.set_index('SwissLipids_ID').to_dict('index')
    nx.set_node_attributes(g, node_attr)
    return g


def MakeGraph(n, edges_list):
    g = nx.DiGraph()
    for edge in edges_list:
        try:
            g.add_edge(edge[0], edge[1])
        except IndexError:
            print(edge)
    return g


def ReadFile(file):
    with open(file, 'r') as f:
        n = int(f.readline())
        m = int(f.readline())
        arr = ReadList(f.readline())
        keys = ReadList(f.readline())
    return n, m, arr, keys


def ReadList(s):
    return list(map(str, s.split()))


def GetData(file):
    edges_list = []

    with open(file, "r") as f:
        n, m = map(int, f.readline().split())

        for i in range(m):
            edges_list.append(ReadList(f.readline()))

    return n, m, edges_list


''' 
    This part is dedicated to getting representatives
'''

















