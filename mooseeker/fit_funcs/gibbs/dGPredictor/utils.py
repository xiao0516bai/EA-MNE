import streamlit as st
import pandas as pd
import numpy as np
import re
from PIL import Image
import webbrowser
import json
import pickle
import sys
import joblib
import pdb

# sys.path.append('./CC/')
# import chemaxon
# from chemaxon import *
# from compound import Compound
# from compound_cacher import CompoundCacher
# sys.path.append('./CC/')

sys.path.append('.')

from dGPredictor.CC.chemaxon import *
from dGPredictor.CC.compound import Compound
from dGPredictor.CC.compound_cacher import CompoundCacher
from rdkit.Chem import rdChemReactions as Reactions
from rdkit.Chem import Draw
from rdkit import Chem

import os
current_dir = os.path.realpath(os.curdir)#获取当前目录的绝对路径



# Bulk Estimation of dG for a list of KEGG Reactions


def load_smiles():
    db = pd.read_csv('./dGPredictor/data/cache_compounds_20160818.csv', # KEGG 化合物的汇集
                     index_col='compound_id')
    db_smiles = db['smiles_pH7'].to_dict()
    return db_smiles


def load_molsig_rad1():
    molecular_signature_r1 = json.load(open('./dGPredictor/data/decompose_vector_ac.json')) ## KEGG 化合物的分解
    return molecular_signature_r1


def load_molsig_rad2():
    molecular_signature_r2 = json.load(
        open('./dGPredictor/data/decompose_vector_ac_r2_py3_indent_modified_manual.json')) ##
    return molecular_signature_r2


def load_model():
    filename = './dGPredictor/model/M12_model_BR.pkl'
    loaded_model = joblib.load(open(filename, 'rb'))
    return loaded_model


def load_compound_cache():
    ccache = CompoundCacher()
    return ccache



def parse_reaction_formula_side(s):
    """
        Parses the side formula, e.g. '2 C00001 + C00002 + 3 C00003'
        Ignores stoichiometry.

        Returns:
            The set of CIDs.
    """
    if s.strip() == "null":
        return {}

    compound_bag = {}
    for member in re.split('\s+\+\s+', s):
        tokens = member.split(None, 1)
        if len(tokens) == 0:
            continue
        if len(tokens) == 1:
            amount = 1
            key = member
        else:
            amount = float(tokens[0])
            key = tokens[1]

        compound_bag[key] = compound_bag.get(key, 0) + amount

    return compound_bag


def parse_formula(formula, arrow='<=>', rid=None):
    """
    Parses a two-sided formula such as: 2 C00001 => C00002 + C00003

    Return:
        The set of substrates, products and the direction of the reaction
    """

    tokens = formula.split(arrow)
    if len(tokens) < 2:
        print(('Reaction does not contain the arrow sign (%s): %s'
               % (arrow, formula)))
    if len(tokens) > 2:
        print(('Reaction contains more than one arrow sign (%s): %s'
               % (arrow, formula)))

    left = tokens[0].strip()
    right = tokens[1].strip()

    sparse_reactioin = {}

    for cid, count in parse_reaction_formula_side(left).items():
        sparse_reactioin[cid] = sparse_reactioin.get(cid, 0) - count

    for cid, count in parse_reaction_formula_side(right).items():
        sparse_reactioin[cid] = sparse_reactioin.get(cid, 0) + count

    return sparse_reactioin



def get_rule(rxn_dict, molsig1, molsig2, novel_decomposed1, novel_decomposed2):
    if novel_decomposed1 != None:
        for cid in novel_decomposed1:
            molsig1[cid] = novel_decomposed1[cid]

    if novel_decomposed2 != None:
        for cid in novel_decomposed2:
            molsig2[cid] = novel_decomposed2[cid]

    molsigna_df1 = pd.DataFrame.from_dict(molsig1).fillna(0)

    molsigna_df1 = pd.DataFrame.from_dict(molsig1).fillna(0)
    all_mets1 = molsigna_df1.columns.tolist()
    all_mets1.append("C00080")
    all_mets1.append("C00282")

    molsigna_df2 = pd.DataFrame.from_dict(molsig2).fillna(0)
    all_mets2 = molsigna_df2.columns.tolist()
    all_mets2.append("C00080")
    all_mets2.append("C00282")

    moieties_r1 = open('./dGPredictor/data/group_names_r1.txt')
    moieties_r2 = open('./dGPredictor/data/group_names_r2_py3_modified_manual.txt')
    moie_r1 = moieties_r1.read().splitlines()
    moie_r2 = moieties_r2.read().splitlines()

    molsigna_df1 = molsigna_df1.reindex(moie_r1)
    molsigna_df2 = molsigna_df2.reindex(moie_r2)

    rule_df1 = pd.DataFrame(index=molsigna_df1.index)
    rule_df2 = pd.DataFrame(index=molsigna_df2.index)
    # for rid, value in reaction_dict.items():
    #     # skip the reactions with missing metabolites
    #     mets = value.keys()
    #     flag = False
    #     for met in mets:
    #         if met not in all_mets:
    #             flag = True
    #             break
    #     if flag: continue

    rule_df1['change'] = 0
    for met, stoic in rxn_dict.items():
        if met == "C00080" or met == "C00282":
            continue  # hydogen is zero
        rule_df1['change'] += molsigna_df1[met] * stoic

    rule_df2['change'] = 0
    for met, stoic in rxn_dict.items():
        if met == "C00080" or met == "C00282":
            continue  # hydogen is zero
        rule_df2['change'] += molsigna_df2[met] * stoic

    rule_vec1 = rule_df1.to_numpy().T
    rule_vec2 = rule_df2.to_numpy().T

    m1, n1 = rule_vec1.shape
    m2, n2 = rule_vec2.shape

    zeros1 = np.zeros((m1, 44))
    zeros2 = np.zeros((m2, 44))
    X1 = np.concatenate((rule_vec1, zeros1), 1)
    X2 = np.concatenate((rule_vec2, zeros2), 1)

    rule_comb = np.concatenate((X1, X2), 1)

    # rule_df_final = {}
    # rule_df_final['rad1'] = rule_df1
    # rule_df_final['rad2'] = rule_df2
    return rule_comb, rule_df1, rule_df2




def get_ddG0(rxn_dict, pH, I, novel_mets):
    ccache = CompoundCacher()
    # ddG0 = get_transform_ddG0(rxn_dict, ccache, pH, I, T)
    T = 298.15
    ddG0_forward = 0
    for compound_id, coeff in rxn_dict.items():
        if novel_mets != None and compound_id in novel_mets:
            comp = novel_mets[compound_id]
        else:
            comp = ccache.get_compound(compound_id)
        ddG0_forward += coeff * comp.transform_pH7(pH, I, T)

    return ddG0_forward


def get_dG0(rxn_dict, rid, pH, I, loaded_model, molsig_r1, molsig_r2,
            novel_decomposed_r1, novel_decomposed_r2, novel_mets):
    # rule_df = get_rxn_rule(rid)
    rule_comb, rule_df1, rule_df2 = get_rule(
        rxn_dict, molsig_r1, molsig_r2, novel_decomposed_r1, novel_decomposed_r2)

    X = rule_comb #先把方程式改成自己设定的规则，然后进行预测

    ymean, ystd = loaded_model.predict(X, return_std=True)

    return ymean[0] + get_ddG0(rxn_dict, pH, I, novel_mets), ystd[0], rule_df1, rule_df2



def get_mu_std(rxn_dict_list, pH=7, I=0.1):
    db_smiles = load_smiles()
    molsig_r1 = load_molsig_rad1()
    molsig_r2 = load_molsig_rad2()

    loaded_model = load_model()
    ccache = load_compound_cache()

    result_dict = {}

    for keys in rxn_dict_list:
        kegg_rxn_string = rxn_dict_list[keys]
        kegg_rxn_dict = parse_formula(kegg_rxn_string)
        mu, std, rule_df1, rule_df2 = get_dG0(kegg_rxn_dict, keys, pH, I, loaded_model, molsig_r1, molsig_r2,
                                              [], [], [])
        result_dict[keys] = [mu, std]

    return result_dict




# KEGG_rxn_list = {"R00010": "C01083 + C00001 <=> 2 C00031",
#                  "R00303": "C00092 + C00001 <=> C00031 + C00009",
#                  "R00304": "C00103 + C00001 <=> C00031 + C00009",
#                  "R07294": "C15524 + C00001 <=> C02137 + C00010",
#                  "R01252": "C00148 + C00026 + C00007 <=> C01157 + C00042 + C00011",
#                  "R00406": "C00091 + C00149 <=> C00042 + C04348"
#                  }



# print(get_mu_std(KEGG_rxn_list))



# pH = 7  # any number between 0-14
# I = 0.1  # min_value=0.0, max_value=0.5)
#
# db_smiles = load_smiles()
# molsig_r1 = load_molsig_rad1()
# molsig_r2 = load_molsig_rad2()

# loaded_model = load_model()
# ccache = load_compound_cache()
#
#
# for keys in KEGG_rxn_list:
#     kegg_rxn_string = KEGG_rxn_list[keys]
#     kegg_rxn_dict = parse_formula(kegg_rxn_string)
#     mu, std, rule_df1, rule_df2 = get_dG0(kegg_rxn_dict, keys, pH, I, loaded_model, molsig_r1, molsig_r2, [], [], [])
#     print(keys)
#     print(kegg_rxn_string)
#     print("dG = %.2f ± %.2f kJ/mol" % (mu, std))