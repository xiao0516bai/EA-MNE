"""
Description:
Autor: caoyh
Date: 2022-09-16 17:02:26
LastEditTime: 2022-11-20 21:25:23
"""

from mooseeker.kegg_helper.pykegg import split_equation
from mooseeker.fit_funcs.gibbs.mdf import *
from mooseeker.utils import read_json


def get_gibbs_mlp(cfg, rxn_list):
    # dG_dict = read_json(cfg['file_path']['dG_dict'])
    dG_dict = read_json(cfg['file_path']['dG_dict_mlp'])
    rxn_dict = read_json(cfg['file_path']['all_rxn_dict'])
    cpd_dict = read_json(cfg['file_path']['all_cpd_dict'])
    cons_list = pd.read_csv(cfg["file_path"]['gibb_concentration'], sep='\t', header=None).to_dict('list')

    reaction = [] # Rxxxxx 反应方程式
    drGs = [] # Rxxxxx dG_Mean
    con = []
    cpb_list = set() # 存储这条链表中所有的化合物

    for rxn in rxn_list:
        reaction += [rxn + '\t' + rxn_dict[rxn]['equation']] # Rxxxxx 反应方程式
        try: 
            # 如果没有这个值 或者没有这个反应，则认为为0
            drGs += [rxn + '\t' + str(dG_dict[rxn]['dG_Mean'])] # Rxxxxx dG_Mean 根据反应id去dG.json中查值
        except Exception as e: # dG_dict 没有 R13013
            print(e)
            drGs += [rxn + '\t' + str(0)] # dG.json中没有这个反应的话 将此反应的吉布斯自由能值视为0

        s_list, p_list = split_equation(rxn_dict[rxn]['equation'])

        for c in s_list+p_list:
            try:
                cpb_list.add(c) # 存储这条链表中所有的化合物
            except Exception as e:
                print(e)
                continue

    for cpb in list(cpb_list): # 遍历所有的化合物
        if cpb in cons_list[0]:
            idx = cons_list[0].index(cpb)
            # min = cons_list[1][idx]
            # max = cons_list[2][idx]
            con += [cpb + '\t' + str(cons_list[1][idx]) + '\t' + str(cons_list[2][idx])] # Cxxxxx cons_list[1][idx] cons_list[2][idx]

    reactions_text = "\n".join(reaction) + "\n"
    drGs_text = "\n".join(drGs) + "\n"
    cons_text = '\n'.join(con)


    # Load stoichiometric matrix, drGs and constraints
    S = read_reactions(reactions_text)
    drGs = read_reaction_drGs(drGs_text)
    constraints = read_constraints(cons_text)
    # Calculate MDF inputs
    c = mdf_c(S)
    A = mdf_A(S)
    b = mdf_b(S, drGs, constraints)
    # Perform MDF 求一条路径的吉布斯自由能值
    mdf_result = mdf(c, A, b)

    return mdf_result.fun














