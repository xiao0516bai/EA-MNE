import os
import json
from pathlib import Path
import cobra
# import cobra.test
from cobra.io import read_sbml_model
from cobra import Reaction, Metabolite
# import pubchempy as pcp
import sys
import numpy as np
# 调用里面的split_equation
from mooseeker.kegg_helper.pykegg import *
import re
# SingleReaction 和 SingleLinkList类的声明 get_config read_json
from mooseeker.utils import *

#############################################################################
## 路径扩展为网络相关函数
#############################################################################

# 初始化酶池 为每个反应随机选择酶

def init_enzyme_pool(chrom,rxn_dict):
    reactions_list = chrom.get_rxn_list()
    enzyme_pool = set()
    for i in range(chrom.length()):
        reaction = reactions_list[i]
        ## 以reaction为键值 去all_rxn_dict中查找酶
        enzymes = rxn_dict[reaction]['enzyme']
        # 酶可能是多个，随机选择一个酶，也可能没有
        if enzymes:
            enzyme = enzymes[np.random.randint(len(enzymes))]
            # enzyme = enzymes[0]
            enzyme_pool.add(enzyme)
    # print("选择酶：",enzyme_pool)
    return enzyme_pool

# 初始化底物池
def init_substrate_pool_from_path(chrom,all_pool):
    reactions_dict = chrom.get_rxn_dict()
    reactions_list = chrom.get_rxn_list()
    compounds_list = chrom.get_cpd_list()
    # print("reactions_list:",reactions_list)
    # print("compounds_list:",compounds_list)
    # print("equations:")
    # for reaction in reactions_list:
    #     print(rxn_dict[reaction]['equation'])

    substrate_list = [] # 初始化底物列表 需要确保初始化的底物列表中没有重复的

    # 先把反应的所有的反应物都加底物列表中
    for i,reaction in enumerate(reactions_list):
        # 划分出来反应物 split_equation函数在pykegg里面定义了 已经引用了
        lefts,rights = split_equation(reactions_dict[reaction])
        ## 判断遍历的这个底物s 是在lefts里面还是 rights里面
        ## 在哪个里面 哪边就是反应物
        s = compounds_list[i]
        if s in lefts:
            reactants = lefts
        elif s in rights:
            reactants = rights
        else:
            print(s,"没有在方程式里面")
        for r in reactants:
            if r not in substrate_list:
                substrate_list.append(r)
    for i in range(1,len(compounds_list)-1):
        substrate_list.remove(compounds_list[i])

    # 初始化底物池
    substrate_pool = {}
    all_pool_cpds = list(all_pool.keys()) 
    for c in substrate_list: 
        if c not in all_pool_cpds: 
            substrate_pool[c] = None
        else:
            temp = list(all_pool[c]) 
            substrate_pool[c] = temp 
           
    return substrate_pool

def init_total_rxn(chrom,rxn_dict):
    total_rxn_list = []
    total_rxn_dict = {}
    reactions_list = chrom.get_rxn_list()
    # reactions_dict = chrom.get_rxn_dict()
    compounds_list = chrom.get_cpd_list()
    for i, key in enumerate(reactions_list):
        total_rxn_dict[key] = rxn_dict[key]

        equation = total_rxn_dict[key]['equation']
        lefts, rights = split_equation(equation)
        s = compounds_list[i]
        if s in lefts:
            total_rxn_dict[key]['direction'] = 0
        elif s in rights:
            total_rxn_dict[key]['direction'] = 1
        else:
            print(s, "没有在方程式里面")
    for rxn in chrom.get_rxn_list():
        total_rxn_list.append(rxn)
    return total_rxn_list, total_rxn_dict

def reactants_are_all_in_substrate_pool(reactants,substrate_pool):
    count = 0
    for reactant in reactants:
        if reactant not in substrate_pool.keys():
            break
        else:
            count+=1
    if(count==len(reactants)):
        return True
    else:
        return False

def add_products_to_substrates_pool(products,need_to_be_added_products,substrate_pool,all_pool_cpds,all_pool):
    for product in products:
        if product not in substrate_pool.keys():
            if product not in need_to_be_added_products.keys():
                if product not in all_pool_cpds: 
                    need_to_be_added_products[product] = None
                else:
                    temp = list(all_pool[product])
                    need_to_be_added_products[product] = temp 
    return need_to_be_added_products

def add_reaction_to_total(reaction,direction,total_rxn_list,total_rxn_dict):
  
    if reaction['kegg_id'] in total_rxn_list:
        pass
    else:
        total_rxn_list.append(reaction['kegg_id'])
        total_rxn_dict[reaction['kegg_id']] = reaction
        total_rxn_dict[reaction['kegg_id']]['direction'] = direction
    return total_rxn_list,total_rxn_dict

def start_having_reactions_in_substrate_pool(substrate_pool,all_pool,enzyme_pool,total_rxn_list,total_rxn_dict):
    ### 遍历底物池
    all_pool_cpds = list(all_pool.keys())
    need_to_be_added_products = {}
    for s in substrate_pool.keys():
#         print("--------------------------------")
#         print("遍历化合物",s,"的所有反应")
        if substrate_pool[s] is None:
            continue
        else:
            singleReactions = substrate_pool[s]# singleReactions是一个列表。包含了某个底物的所有反应
#             print(s,"包含的所有反应个数：",len(singleReactions))
            for singleReaction in singleReactions:
                reaction = singleReaction.reaction 
                enzymes = reaction['enzyme']
    #             print(reaction)
                if enzymes is not None:
                    for e in enzymes:
                        # 有催化反应的酶
                        if e in enzyme_pool:
                            # 判断底物池中是否存在所有的反应物
                            lefts,rights = split_equation(reaction['equation'])
                            if s in lefts:
                                reactants = lefts
                                products = rights
                                direction = 0
                            elif s in rights:
                                reactants = rights
                                products = lefts
                                direction = 1
                            else:
                                print(s,"没有在方程式里面")
                            if reactants_are_all_in_substrate_pool(reactants,substrate_pool):
    #                             print(reaction)
                                total_rxn_list,total_rxn_dict = add_reaction_to_total(reaction,direction,total_rxn_list,total_rxn_dict)
                               
                                need_to_be_added_products = add_products_to_substrates_pool(products,need_to_be_added_products,substrate_pool,all_pool_cpds,all_pool)
                            break; 
                else:
                    # 如果没有酶，跳过这条反应
                    # print(f"反应 {reaction['equation']} 没有酶或酶不在酶池中，跳过此反应。")
                    continue
    return need_to_be_added_products,total_rxn_list,total_rxn_dict,substrate_pool

def react_thoroughly(substrate_pool,all_pool,enzyme_pool,total_rxn_list,total_rxn_dict):
    while True:
        need_to_be_added_products = {}
        need_to_be_added_products,total_rxn_list,total_rxn_dict,substrate_pool = start_having_reactions_in_substrate_pool(substrate_pool,all_pool,enzyme_pool,total_rxn_list,total_rxn_dict)
        if len(need_to_be_added_products)!= 0:
            substrate_pool.update(need_to_be_added_products)
        else:
            return total_rxn_list,total_rxn_dict,substrate_pool

def traverse_reactions_in_path(chrom,substrate_pool,all_pool,enzyme_pool,total_rxn_list,total_rxn_dict):
    reactions_dict = chrom.get_rxn_dict()
    reactions_list = chrom.get_rxn_list()
    compounds_list = chrom.get_cpd_list()
    all_pool_cpds = list(all_pool.keys())
    for i,reaction in enumerate(reactions_dict):
        equation = reactions_dict[reaction]
        lefts,rights = split_equation(equation)
        s = compounds_list[i]
        if s in lefts:
            reactants = lefts
            products = rights
        elif s in rights:
            reactants = rights
            products = lefts
        else:
            print(s,"没有在方程式里面")
        added_products = {}

        for product in products:
            added_products = add_products_to_substrates_pool(products,added_products,substrate_pool,all_pool_cpds,all_pool)

        substrate_pool.update(added_products)

        total_rxn_list,total_rxn_dict,substrate_pool = react_thoroughly(substrate_pool,all_pool,enzyme_pool,total_rxn_list,total_rxn_dict)
    return total_rxn_list,total_rxn_dict,substrate_pool

def get_network(cfg, chrom,log,i):

    rxn_dict = read_json(cfg['file_path']['all_rxn_dict'])
    # 获取所有化合物
    cpd_dict = read_json(cfg['file_path']['all_cpd_dict'])
    all_pool = np.load(cfg['file_path']['mypool'], allow_pickle=True).item()


    # rxn_list = chrom.get_rxn_list()
    # cpd_list = chrom.get_cpd_list()
    # print("get_network_yield -> rxn_list:",rxn_list)
    # print("get_network_yield -> cpd_list:", cpd_list)
    # 初始化酶池
    log.logger.info('--- %d expanding the path into a network' %i)
    # log.logger.info('--- (fit_funcs) get enzyme_pool ---')
    # enzyme_pool = init_enzyme_pool(chrom,rxn_dict)[i]  # 初始化酶池
    enzyme_pool = init_enzyme_pool_list(chrom,rxn_dict)[i]  # 初始化酶池
    # log.logger.info(enzyme_pool)

    # 初始化底物池
    # log.logger.info('--- (fit_funcs) get substrate_pool from path---')
    substrate_pool = init_substrate_pool_from_path(chrom,all_pool)  # 初始化底物池 不包含主路径中除了起始化合物之外的化合物

    # 初始化总反应池
    # log.logger.info('--- (fit_funcs) init_total_rxn ---')
    total_rxn_list, total_rxn_dict = init_total_rxn(chrom,rxn_dict)

    # 初始底物池反应发生彻底
    # log.logger.info('--- (fit_funcs) react_thoroughly ---')
    total_rxn_list,total_rxn_dict,substrate_pool = react_thoroughly(substrate_pool,all_pool,enzyme_pool,total_rxn_list,total_rxn_dict)  # 初始底物池中的反应彻底发生完

    # 遍历路径上的反应 不断添加反应所需的反应物 扩展底物池 使底物池反应发生彻底
    # log.logger.info('--- (fit_funcs) traverse_reactions_in_path ---')
    total_rxn_list,total_rxn_dict,substrate_pool = traverse_reactions_in_path(chrom,substrate_pool,all_pool,enzyme_pool,total_rxn_list,total_rxn_dict)  # 依次把主路径上的反应的生成物添加到底物池中，然后使其彻底反应
    # log.logger.info(total_rxn_list)

    ####################  计算网络的理论产量值 ####################
    # 先默认返回10 方便调试
    return total_rxn_list, total_rxn_dict, substrate_pool,enzyme_pool

# 初始化酶池
def init_enzyme_pool_list(chrom,rxn_dict):
    # 反应列表
    reactions_list = chrom.get_rxn_list()
    return generate_enzyme_combinations(reactions_list,rxn_dict, current_combination=[], index=0, result=[])

def generate_enzyme_combinations(reactions_list,rxn_dict, current_combination=[], index=0, result=[]):
    if index == len(reactions_list):
        if current_combination:
            result.append(current_combination.copy())
        return
    reaction = reactions_list[index]
    enzymes = rxn_dict[reaction]['enzyme']
    if enzymes is not None:
        for enzyme in enzymes:
            current_combination.append(enzyme)
            generate_enzyme_combinations(reactions_list,rxn_dict, current_combination, index + 1, result)
            current_combination.pop()
    else:
        generate_enzyme_combinations(reactions_list, rxn_dict, current_combination, index + 1, result)
    return result

# 并且记录当前的酶选择 需要写到文档里面 （最终肯定要记录的，因为如果对酶筛选的话，需要知道酶信息）
def get_network_totalReactions(cfg: object, chrom: object, log: object) -> object:
    rxn_list = chrom.chrom.get_rxn_list()
    cpd_list = chrom.chrom.get_cpd_list()
    log.logger.info("get_network_totalReactions -> original-rxn_list:")
    log.logger.info(rxn_list)
    log.logger.info("get_network_totalReactions -> original-cpd_list:")
    log.logger.info(cpd_list)

    rxn_dict = read_json(cfg['file_path']['all_rxn_dict'])
    cpd_dict = read_json(cfg['file_path']['all_cpd_dict'])
    enzyme_list = init_enzyme_pool_list(chrom.chrom, rxn_dict)
    count = sys.maxsize  # 设定网络中反应总数最大值100000
    best_enzyme_pool = []
    best_total_rxn_list = []

    best_total_rxn_dict = []

    # 最好的底物池 之前没有保留
    best_substrate_pool = []

    for i in range(len(enzyme_list)):
        total_rxn_list, total_rxn_dict, substrate_pool, enzyme_pool = get_network(cfg, chrom.chrom,log, i)
        if len(total_rxn_list) < count:
            count = len(total_rxn_list)
            best_enzyme_pool = enzyme_pool
            best_total_rxn_list = total_rxn_list
            best_total_rxn_dict = total_rxn_dict
            best_substrate_pool = substrate_pool
    # log 输出
    log.logger.info("For this path,best_enzyme_pool:")
    log.logger.info(best_enzyme_pool)
    log.logger.info("Min-count: %d" %count)
    # log.logger.info("-"*20)
    # 返回适应度值为，网络中最小的反应数
    # 这个酶池需要保存的 在chromesome中添加了一个 best_enzyme_pool 属性
    # chrom.set_best_network(best_enzyme_pool,best_total_rxn_list,count)
    print("get_network_totalReactions--->")
    print("best_enzyme_pool:")
    print(best_enzyme_pool)
    chrom.best_enzyme_pool = best_enzyme_pool
    print("chrom.best_enzyme_pool:")
    print(chrom.best_enzyme_pool)
    chrom.total_reactions = count
    chrom.best_total_rxn_list = best_total_rxn_list
    chrom.best_total_rxn_dict = best_total_rxn_dict
    chrom.best_substrate_pool = best_substrate_pool

    # # 所有的化合物
    # total_compounds = []
    # for key in substrate_pool.keys():
    #     total_compounds.append(key)

    # 所有的化合物
    total_compounds = []
    for key in best_substrate_pool.keys():
        total_compounds.append(key)

    chrom.total_compounds = total_compounds

    # 记录这个网络扩展后的总反应 方便之后画图
    return count

## 毒性评估
# def get_toxicity(cfg,chrom,log):
#
#     if chrom.total_compounds == []:
#         count = get_network_totalReactions(cfg,chrom,log)
#         # 重点不是得到count 而是得到最好的网络 然后再进行毒性评估
#         # 计算所有化合物的毒性。  负值有毒，正值无毒。
#         # 因为 适应度值需要越小越好。所以对于有毒的负值，最后要进行取反。然后累加。
#     total_compounds = chrom.total_compounds
#     cpd_dict = read_json(cfg['file_path']['all_cpd_dict'])
#
#     toxicity_sum = 0
#
#     for c in total_compounds:
#         print(cpd_dict[c]['toxicity'])
#         if (cpd_dict[c]['toxicity'] < 0):
#             toxicity_sum += cpd_dict[c]['toxicity']
#
#     return -toxicity_sum

def get_toxicity(cfg,chrom,log,way=0):

    # way
    # 0 代表是评估网络，默认是网络
    # 1 代表评估是路径

    if way==0:
        print('---- 计算网络中化合物的毒性 ----')
        # 网络，需要扩展网络
        if chrom.total_compounds == []:
            count = get_network_totalReactions(cfg,chrom,log)
            # 重点不是得到count 而是得到最好的网络 然后再进行毒性评估
            # 计算所有化合物的毒性。  负值有毒，正值无毒。
            # 因为 适应度值需要越小越好。所以对于有毒的负值，最后要进行取反。然后累加。
        total_compounds = chrom.total_compounds
    else:
        print('---- 计算路径中化合物的毒性 ----')
        # 路径，可以直接获取化合物列表
        total_compounds = chrom.chrom.get_cpd_list()

    cpd_dict = read_json(cfg['file_path']['all_cpd_dict'])

    toxicity_sum = 0

    for c in total_compounds:
        print(c,"---->",cpd_dict[c]['toxicity'])
        if cpd_dict[c]['toxicity'] is None:
            cpd_dict[c]['toxicity'] = 0
        if (cpd_dict[c]['toxicity'] < 0):
            toxicity_sum += cpd_dict[c]['toxicity']

    return -toxicity_sum


####################### 网络评估f1，f2，f3所需的所有函数 ##########################################
# 定义多叉树的树节点
class MyTreeNode:
    def __init__(self, compound):
        self.compound = compound  # 节点的值
        self.children = []  # 子节点列表

    def add_child(self, child_node):
        """添加一个子节点，如果该子节点不在children中"""
        if child_node not in self.children:
            self.children.append(child_node)


# # 表示每个树节点，放在字典列表中
# def get_node_dict(final_compounds, all_pool, total_rxn_dict, myCompounds, main_reactions, debug=False):
#     # 去除辅因子的化合物
#     final_compounds_no_confactors = [c for c in final_compounds if c in all_pool]
#     node_dict = {}

#     # 记录所有参与主路径反应的化合物，但不在主路径化合物列表中
#     involved_in_main_reactions = set()
#     for r in main_reactions:
#         reaction = total_rxn_dict[r]
#         direction = reaction["direction"]
#         reactants = reaction["reactants"]
#         products = reaction["products"]

#         # 双向检查以确保化合物完整性
#         for c in reactants + products:
#             if c not in myCompounds:
#                 involved_in_main_reactions.add(c)

#     # 初始化节点字典
#     for c in final_compounds_no_confactors:
#         node_dict[c] = MyTreeNode(c)

#     # 构建节点关系
#     for c in final_compounds_no_confactors:
#         if c in involved_in_main_reactions:
#             continue

#         cNode = node_dict[c]
#         for r in all_pool.get(c, []):
#             r_id = r.reaction.get("kegg_id")
#             if r_id in total_rxn_dict.keys():
#                 nextCompounds = []
#                 # 根据反应方向确定下一个化合物
#                 if total_rxn_dict[r_id]["direction"] == 0 and c not in total_rxn_dict[r_id]["products"]:
#                     nextCompounds = total_rxn_dict[r_id]["products"]
#                 elif total_rxn_dict[r_id]["direction"] == 1 and c not in total_rxn_dict[r_id]["reactants"]:
#                     nextCompounds = total_rxn_dict[r_id]["reactants"]

#                 # 添加子节点
#                 for nextC in nextCompounds:
#                     if (
#                         nextC in final_compounds_no_confactors
#                         and nextC != c
#                         and nextC not in myCompounds
#                         and nextC not in involved_in_main_reactions
#                     ):
#                         cNode.add_child(node_dict[nextC])  # 添加 MyTreeNode 对象

#     if debug:
#         print("构建的节点字典:")
#         for k, v in node_dict.items():
#             children = [child.compound for child in v.children]
#             print(f"{k}: 子节点 -> {children}")

#     return node_dict
class MyTreeNode:
    def __init__(self, compound):
        self.compound = compound  # 节点的值
        self.children = []  # 子节点列表

    def add_child(self, child_node):
        """添加一个子节点，如果该子节点不在children中"""
        if child_node not in self.children:
            self.children.append(child_node)


# 表示每个树节点，放在字典列表中
def get_node_dict(final_compounds, all_pool, total_rxn_dict, myCompounds, main_reactions, debug=False):
    # 去除辅因子的化合物
    final_compounds_no_confactors = [c for c in final_compounds if c in all_pool]
    node_dict = {}

    # 记录所有参与主路径反应的化合物，但不在主路径化合物列表中
    involved_in_main_reactions = set()
    for r in main_reactions:
        reaction = total_rxn_dict[r]
        direction = reaction["direction"]
        reactants = reaction["reactants"]
        products = reaction["products"]

        # 双向检查以确保化合物完整性
        for c in reactants + products:
            if c not in myCompounds:
                involved_in_main_reactions.add(c)

    # 初始化节点字典
    for c in final_compounds_no_confactors:
        if c not in node_dict:  # 检查化合物是否已存在，避免重复创建
            node_dict[c] = MyTreeNode(c)

    # 构建节点关系
    for c in final_compounds_no_confactors:
        if c in involved_in_main_reactions:
            continue

        cNode = node_dict[c]
        for r in all_pool.get(c, []):
            r_id = r.reaction.get("kegg_id")
            if r_id in total_rxn_dict.keys():
                nextCompounds = []
                # 根据反应方向确定下一个化合物
                if total_rxn_dict[r_id]["direction"] == 0 and c not in total_rxn_dict[r_id]["products"]:
                    nextCompounds = total_rxn_dict[r_id]["products"]
                elif total_rxn_dict[r_id]["direction"] == 1 and c not in total_rxn_dict[r_id]["reactants"]:
                    nextCompounds = total_rxn_dict[r_id]["reactants"]

                # 添加子节点
                for nextC in nextCompounds:
                    # 检查 nextC 是否是 final_compounds_no_confactors，并且 nextC 不在 myCompounds 中
                    if (
                        nextC in final_compounds_no_confactors
                        and nextC != c
                        and nextC not in myCompounds
                        and nextC not in involved_in_main_reactions
                    ):
                        # 确保 nextC 已经在 node_dict 中
                        if nextC not in node_dict:
                            node_dict[nextC] = MyTreeNode(nextC)
                        cNode.add_child(node_dict[nextC])  # 添加 MyTreeNode 对象

    # 打印调试信息
    if debug:
        print("构建的节点字典:")
        for k, v in node_dict.items():
            children = [child.compound for child in v.children]
            print(f"{k}: 子节点 -> {children}")

    return node_dict


# # 深度优先算法递归计算每个根节点的分支权重加和
# def compute_branch(compound, node_dict, myCompounds, line, visited, debug=False):
#     if compound in visited:
#         return 0  # 防止重复计算
#     visited.add(compound)

#     # 当前节点的权重
#     line_score = 1 / line
#     branch_fitness = line_score * len(node_dict[compound].children)

#     if debug:
#         print(f"节点 {compound} 在层级 {line} 的权重为 {line_score:.2f}，分支数为 {len(node_dict[compound].children)}，分支和为 {branch_fitness:.2f}")
#         print(f"分支节点: {[child.compound for child in node_dict[compound].children]}")

#     # 递归计算子节点
#     for child in node_dict[compound].children:
#         if child.compound not in myCompounds:
#             branch_fitness += compute_branch(child.compound, node_dict, myCompounds, line + 1, visited, debug)

#     return branch_fitness


# # 遍历每个主路径上的化合物得到整条路径的分支权重和 f2
# def get_f2(myCompounds, node_dict, debug=False):
#     branch_fitness = 0
#     visited = set()
#     for i, c in enumerate(myCompounds, start=1):
#         score = compute_branch(c, node_dict, myCompounds, line=1, visited=visited, debug=debug)
#         if debug:
#             print(f"主路径节点 {c} 的分支和: {score:.2f}")
#         branch_fitness += score / i

#     return branch_fitness

def compute_branch(compound, node_dict, myCompounds, line, visited_nodes_global, level_fitness, debug=True):
    
    # 检查 compound 是否已经被访问过
    if compound in visited_nodes_global:
        print(f"跳过节点 {compound}（已经访问过），权重设为 0")  # 添加日志帮助调试
        return 0  # 防止重复计算，直接返回 0 权重

    # 将当前节点标记为已访问
    visited_nodes_global.add(compound)

    # 当前层的分数按照 1 / line 计算
    line_score = 1 / line
    yield_branch_fitness = 0  # 初始化当前节点的分支权重

    # 计算当前节点的分支权重
    for child in node_dict[compound].children:
        # 如果子节点已经被访问过，跳过这个子节点的分支计算
        if child.compound in visited_nodes_global:
            if debug:
                print(f"跳过子节点 {child.compound}（已经访问过）")
            continue

        if child.compound in myCompounds:
            # 如果子节点在主路径化合物中，将分支权重设为负值
            yield_branch_fitness -= line_score
        else:
            # 子节点不在主路径中，按正常逻辑累加
            yield_branch_fitness += line_score

    # 更新当前层的分支总和
    if line not in level_fitness:
        level_fitness[line] = 0
    level_fitness[line] += yield_branch_fitness

    # 打印当前节点的分支总和和子节点
    child_names = [child.compound for child in node_dict[compound].children]
    if debug:
        print(f"节点 {compound} 的分支总和为 {yield_branch_fitness:.2f}")
        print(f"节点 {compound} 的子节点: {child_names}")

    # 遍历所有子节点，跳过主路径节点
    for child in node_dict[compound].children:
        if child.compound in myCompounds:
            continue
        # 递归计算子节点的分支
        yield_branch_fitness += compute_branch(child.compound, node_dict, myCompounds, line + 1, visited_nodes_global, level_fitness, debug)

    return yield_branch_fitness


# 动态生成主路径层级的权重
def generate_main_path_weights(num_levels):
    level_weights = [1 / i for i in range(1, num_levels + 1)]
    print("主路径层级的权重:", level_weights)
    return level_weights


def get_f2(myCompounds, node_dict, debug=True):
    level_weights = generate_main_path_weights(len(myCompounds))
    visited_nodes_global = set()  # 创建一个全局访问集合
    level_fitness = {}  # 创建字典记录每层的分支总和
    yield_branch_fitness = 0

    # 使用 reversed 使得从最后一个元素开始遍历 myCompounds
    for i, compound in enumerate(reversed(myCompounds), start=1):
        current_weight = level_weights[-i]  # 对应的层级权重
        score = compute_branch(compound, node_dict, myCompounds, line=1, visited_nodes_global=visited_nodes_global, level_fitness=level_fitness)
        print(f"主路径节点 {compound} ：{score:.2f}，层级权重：{current_weight:.2f}")
        yield_branch_fitness += score * current_weight

    # 打印每一层的分支总和
    for level in sorted(level_fitness.keys()):
        print(f"层级 {level} 的分支总和：{level_fitness[level]:.2f}")

    print(f"总体分支评估f2: {yield_branch_fitness}")

    return yield_branch_fitness



# 计算指向主路径上化合物的反应的个数 f3
def get_f3(myCompounds, node_dict):
    l = len(myCompounds) - 1
    num = 0
    for c in node_dict:
        cNode = node_dict[c]
        for child in cNode.children:
            if child.compound in myCompounds:
                num += 1
    f3 = num
    return -f3


# 网络中的总反应个数 f1
def get_f1(total_rxn_dict):
    return len(total_rxn_dict)


# 代谢网络分支量评估的总结果
def get_network_branch_fitness(cfg, chrom, log, debug=True):
    log.logger.info("开始代谢网络分支评估:")

    # 主路径上的化合物和反应列表
    myCompounds = chrom.chrom.get_cpd_list()
    main_reactions = chrom.chrom.get_rxn_list()

    if not chrom.total_compounds:
        count = get_network_totalReactions(cfg, chrom, log)

    final_compounds = chrom.total_compounds
    total_rxn_dict = chrom.best_total_rxn_dict
    all_pool = np.load(cfg["file_path"]["mypool"], allow_pickle=True).item()

    w1, w2, w3 = 0.5, 0.3, 0.2

    # 构建节点字典
    # node_dict = get_node_dict(final_compounds, all_pool, total_rxn_dict, myCompounds, main_reactions, debug)
    node_dict = get_node_dict(chrom.total_compounds, all_pool, total_rxn_dict, myCompounds, main_reactions)

    f1 = get_f1(total_rxn_dict)
    f2 = get_f2(myCompounds, node_dict, debug)
    f3 = get_f3(myCompounds, node_dict)
    
    print(f" 代谢网络分支量评估f1: {f1}, f2: {f2}, f3: {f3}")

    return w1 * f1 + w2 * f2 + w3 * f3



##############################计算理论产率的权重值###################################################

def compute_yield_branch(compound, node_dict, myCompounds, line, visited_nodes_global, level_fitness):
    """
    计算以 compound 为根的分支的 fitness。
    如果分支的子节点在主路径化合物列表内，将该分支的权重设为负值。
    """
    # 检查 compound 是否已经被访问过
    if compound in visited_nodes_global:
        print(f"跳过节点 {compound}（已经访问过），权重设为 0")  # 添加日志帮助调试
        return 0  # 防止重复计算，直接返回 0 权重

    # 将当前节点标记为已访问
    visited_nodes_global.add(compound)

    # 当前层的分数按照 1 / line 计算
    line_score = 1 / (2 ** line)
    yield_branch_fitness = 0  # 初始化当前节点的分支权重

    # 计算当前节点的分支权重
    for child in node_dict[compound].children:
        # 如果子节点已经被访问过，跳过这个子节点的分支计算
        if child.compound in visited_nodes_global:
            print(f"跳过子节点 {child.compound}（已经访问过）")  # 添加日志帮助调试
            continue

        if child.compound in myCompounds:
            # 如果子节点在主路径化合物中，将分支权重设为负值
            yield_branch_fitness -= line_score
        else:
            # 子节点不在主路径中，按正常逻辑累加
            yield_branch_fitness += line_score

    # 更新当前层的分支总和
    if line not in level_fitness:
        level_fitness[line] = 0
    level_fitness[line] += yield_branch_fitness

    # 打印当前节点的分支总和和子节点
    child_names = [child.compound for child in node_dict[compound].children]
    print(f"节点 {compound} 的分支总和为 {yield_branch_fitness:.2f}")
    print(f"节点 {compound} 的子节点: {child_names}")

    # 遍历所有子节点，跳过主路径节点
    for child in node_dict[compound].children:
        if child.compound in myCompounds:
            continue
        # 递归计算子节点的分支
        yield_branch_fitness += compute_yield_branch(child.compound, node_dict, myCompounds, line + 1, visited_nodes_global, level_fitness)

    return yield_branch_fitness



# 动态生成主路径层级的权重
def generate_main_path_yield_weights(num_levels):
    level_weights_yield = [1 / i for i in range(1, num_levels + 1)]
    print("主路径层级的权重:", level_weights_yield)
    return level_weights_yield

def get_yield_weight(myCompounds, node_dict):
    level_weights_yield = generate_main_path_yield_weights(len(myCompounds))
    visited_nodes_global = set()  # 创建一个全局访问集合
    level_fitness = {}  # 创建字典记录每层的分支总和
    yield_branch_fitness = 0

    # 使用 reversed 使得从最后一个元素开始遍历 myCompounds
    for i, compound in enumerate(reversed(myCompounds), start=1):
        current_weight = level_weights_yield[-i]  # 对应的层级权重
        score = compute_yield_branch(compound, node_dict, myCompounds, line=1, visited_nodes_global=visited_nodes_global, level_fitness=level_fitness)
        print(f"主路径节点 {compound} 的分支理论产率：{score:.2f}，层级权重：{current_weight:.2f}")
        yield_branch_fitness += score * current_weight

    # 打印每一层的分支总和
    for level in sorted(level_fitness.keys()):
        print(f"层级 {level} 的分支总和：{level_fitness[level]:.2f}")

    print(f"总体分支理论损耗权重: {yield_branch_fitness}")

    return yield_branch_fitness

def add_cpd2rxn(model, rxn, total_rxn_dict, cpd_dict):

    r = total_rxn_dict[rxn]
    # parser_equation_coef,split_equation 调用的是 师姐在pykegg中写的函数
    rxn_coff = parser_equation_coef(r['equation']) # 得到 化合物和 化学计量系数 字典
    reactants, products = split_equation(r['equation'])
    
    names = locals()
    
    meta_dict = {}
    for cpd in reactants:
        cpd = re.search(r'.*(C\d{5}).*', cpd).group(1)
        # c = MyCompound(cpd)
        c = cpd_dict[cpd]
        
        if c['bigg_id']!=None:
            for bid in c['bigg_id']:
                if Metabolite(bid+'_c') in model.metabolites:
                        # meta_dict[model.metabolites.get_by_id(bid+'_c')] = -1
                        meta_dict[model.metabolites.get_by_id(bid+'_c')] = rxn_coff[cpd]
                        break
                else: ## 这个必须得加 不加就会出错
                    c_name = bid + '_c'
                    names[c_name] = Metabolite(id = c_name, 
                                                compartment='c', 
                                                name=c_name, 
                                                formula=c['formula']) 
                    # meta_dict[names.get(c_name)] = -1 
                    meta_dict[names.get(c_name)] = rxn_coff[cpd]
                    break

        else:
            # print('Something is error. Mannully add compound %s to metabolites ' %(cpd))
            c_name = c['kegg_id'] + '_' + str(c['pubchem'])
            names[c_name] = Metabolite(id = cpd + '_c', 
                                        compartment='c', 
                                        name=c_name, 
                                        formula=c['formula']) 
            # meta_dict[names.get(c_name)] = -1 
            meta_dict[names.get(c_name)] = rxn_coff[cpd]
      
    for cpd in products:
        cpd = re.search(r'.*(C\d{5}).*', cpd).group(1)
        # c = MyCompound(cpd)
        c = cpd_dict[cpd]

        if c['bigg_id']!=None:
            for bid in c['bigg_id']:
                if Metabolite(bid+'_c') in model.metabolites:
                        # meta_dict[model.metabolites.get_by_id(bid+'_c')] = 1
                        meta_dict[model.metabolites.get_by_id(bid+'_c')] = rxn_coff[cpd]
                        break
                else:
                    c_name = bid + '_c'
                    names[c_name] = Metabolite(id = c_name, 
                                                compartment='c', 
                                                name=c_name, 
                                                formula=c['formula']) 
                    # meta_dict[names.get(c_name)] = -1 
                    meta_dict[names.get(c_name)] = rxn_coff[cpd]
                    break

        else:
            # print('Something is error. Mannully add compound %s to metabolites ' %(cpd))
            c_name = c['kegg_id'] + '_' + str(c['pubchem'])
            names[c_name] = Metabolite(id = cpd + '_c', 
                                        compartment='c', 
                                        name=c_name, 
                                        formula=c['formula']) 
            # meta_dict[names.get(c_name)] = 1  
            meta_dict[names.get(c_name)] = rxn_coff[cpd]

    return meta_dict

def get_yield(cfg, chrom, product, log, way=0):
    """
    使用 chrom 和 log 输入计算理论产率的主函数。

    参数：
    - cfg: 配置文件，包括路径信息。
    - chrom: 包含主路径数据的 Chromosome 对象。
    - product: 目标产物 KEGG 编号。
    - log: 日志记录对象。
    - way: 指标选项，0 表示乘以 (1 - yield_weight)，1 表示直接返回原始值。

    返回：
    - maximum_yield: 计算得到的理论产率。
    """
    log.logger.info("开始计算理论产率...")

    # 加载必要数据
    cpd_dict = read_json(cfg['file_path']['all_cpd_dict'])
    all_pool = np.load(cfg["file_path"]["mypool"], allow_pickle=True).item()
    total_rxn_dict = chrom.best_total_rxn_dict

    myCompounds = chrom.chrom.get_cpd_list()
    main_reactions = chrom.chrom.get_rxn_list()
    node_dict = get_node_dict(chrom.total_compounds, all_pool, total_rxn_dict, myCompounds, main_reactions)

    original = read_sbml_model(cfg['file_path']['host_sbml'])  # 读取模型文件
    model = original.copy()

    # 设置模型的培养基条件
    medium = model.medium
    medium["EX_o2_e"] = 20.0  # 设置葡萄糖和氧气的吸收速率为 20
    medium["EX_glc__D_e"] = 20.0
    model.medium = medium

    # 计算生物量
    wt_growth = model.optimize()
    max_growth = wt_growth.objective_value
    min_growth = 0.8 * max_growth  # 最小碳通量为野生型生长速率的 80%

    # 获取理论产率权重
    yield_weight = get_yield_weight(myCompounds, node_dict)

    # 处理反应列表
    for rxn in main_reactions:
        if rxn not in total_rxn_dict.keys():
            log.logger.error(f"反应 {rxn} 不在 total_rxn_dict 中。")
            return 0

        r = total_rxn_dict[rxn]
        direction = r['direction']

        if r['bigg_id'] is not None:
            found_in_model = False
            for rid in r['bigg_id']:
                if rid in model.reactions:  # 检查反应是否已存在于模型中
                    log.logger.info(f"反应已存在：{rxn} {rid}")
                    found_in_model = True
            if not found_in_model:
                # 添加新反应
                rxn_obj = cobra.Reaction(r['bigg_id'][0])
                rxn_obj.lower_bound = 0 if direction == 0 else -1000
                rxn_obj.upper_bound = 1000 if direction == 0 else 0
                rxn_obj.add_metabolites(add_cpd2rxn(model, rxn, total_rxn_dict, cpd_dict))
                model.add_reactions([rxn_obj])
                log.logger.info(f"添加反应：{rxn} {r['bigg_id'][0]}")
        else:
            # 添加新反应
            rxn_obj = cobra.Reaction(rxn)
            rxn_obj.lower_bound = 0 if direction == 0 else -1000
            rxn_obj.upper_bound = 1000 if direction == 0 else 0
            rxn_obj.add_metabolites(add_cpd2rxn(model, rxn, total_rxn_dict, cpd_dict))
            model.add_reactions([rxn_obj])
            log.logger.info(f"添加反应：{rxn}")

    # # 获取目标产物的 bigg_id
    # target_bigg_id = None
    # if product in cpd_dict:
    #     product_info = cpd_dict[product]
    #     if product_info['bigg_id'] is not None:
    #         target_bigg_id = product_info['bigg_id'][0]

    # if target_bigg_id is None:
    #     log.logger.error(f"未找到 {product} 对应的 bigg_id。")
    #     return -999

    # try:
    #     target_metabolite = model.metabolites.get_by_id(target_bigg_id + "_c")
    #     target_metabolite_id_e = target_metabolite.id.replace("_c", "_e")
    # except KeyError:
    #     log.logger.error(f"未在模型中找到目标产物：{target_bigg_id}")
    #     return -999

    # # 创建交换反应
    # exchange_rxn = cobra.Reaction(f"EX_{target_metabolite_id_e}")
    # exchange_rxn.lower_bound = -1000
    # exchange_rxn.upper_bound = 1000
    # exchange_rxn.add_metabolites({target_metabolite: -1})
    target_bigg_id = None
    if product in cpd_dict:
        product_info = cpd_dict[product]
        if product_info['bigg_id'] is not None:
            target_bigg_id = product_info['bigg_id'][0]  # 使用第一个bigg_id

    if target_bigg_id is None:
        print(f"未找到 {product} 对应的 bigg_id。")
        return original, model, -999

    try:
        target_metabolite = model.metabolites.get_by_id(target_bigg_id + "_c")
    except KeyError:
        print(f"未在模型中找到目标产物：{target_bigg_id}")
        return original, model, -999
    # 创建交换反应（目标产物的出口反应）
    exchange_rxn = cobra.Reaction(f"EX_{target_metabolite.id}")
    exchange_rxn.lower_bound = -1000  # 允许目标产物的排出
    exchange_rxn.upper_bound = 1000
    exchange_rxn.add_metabolites({target_metabolite: -1})  # 目标产物离开细胞
    model.add_reactions([exchange_rxn])

    # 将交换反应设为目标反应
    model.objective = exchange_rxn
    print(model.objective)

    # 优化模型，计算理论产量
    solution = model.optimize()
    max_biomass = solution.objective_value
    log.logger.info(f"最大生物量：{max_biomass}")

    if max_biomass < min_growth:
        log.logger.warning("This pathway is not feasible!")
        return -999
    else:
        if model.reactions.get_by_id('EX_glc__D_e').flux == 0:
            log.logger.warning("This pathway cannot get EX_glc__D_e, Maximum theoretical yield = 0")
            return 0
        else:
            maximum_yield = max_biomass / (-1 * model.reactions.get_by_id('EX_glc__D_e').flux)
            if way == 0:
                maximum_yield *= (1 - yield_weight)
                if maximum_yield < 0:
                    log.logger.warning("Calculated yield is negative, returning -999")
                    return -999 
            log.logger.info(f'理论产率 = {maximum_yield} mmol-Dgly/mmol-glc')
            return maximum_yield
def get_path_yield(cfg, chrom, product, log, way=0):
    """
    直接基于主路径计算理论产率，不进行网络扩展。
    
    参数：
    - cfg: 配置文件，包括路径信息。
    - chrom: 包含主路径数据的 Chromosome 对象。
    - product: 目标产物 KEGG 编号。
    - log: 日志记录对象。
    - way: 指标选项（暂未使用，保留接口）。
    
    返回：
    - maximum_yield: 计算得到的理论产率。
    """
    log.logger.info("开始计算理论产率（不进行网络扩展）...")

    # 加载必要数据
    cpd_dict = read_json(cfg['file_path']['all_cpd_dict'])
    original = read_sbml_model(cfg['file_path']['host_sbml'])  # 读取模型文件
    model = original.copy()

    # 设置模型的培养基条件
    medium = model.medium
    medium["EX_o2_e"] = 20.0  # 设置氧气吸收速率为 20
    medium["EX_glc__D_e"] = 20.0  # 设置葡萄糖吸收速率为 20
    model.medium = medium

    # 计算野生型生物量
    wt_growth = model.optimize()
    max_growth = wt_growth.objective_value
    min_growth = 0.8 * max_growth  # 最小生长速率为野生型的 80%

    # 获取主路径的反应列表
    main_reactions = chrom.chrom.get_rxn_list()
    total_rxn_dict = read_json(cfg['file_path']['all_rxn_dict'])  # 加载所有反应字典

    # 将主路径反应添加到模型中
    for rxn in main_reactions:
        if rxn not in total_rxn_dict:
            log.logger.error(f"反应 {rxn} 不在 total_rxn_dict 中。")
            return -999

        r = total_rxn_dict[rxn]
        if 'direction' not in r:
            log.logger.warning(f"反应 {rxn} 缺少 'direction' 键，使用默认值 0（不可逆）")
            direction = 1
        else:
            direction = r['direction']

        if r['bigg_id'] is not None:
            # 检查反应是否已存在于模型中
            found_in_model = False
            for rid in r['bigg_id']:
                if rid in model.reactions:
                    log.logger.info(f"反应已存在：{rxn} {rid}")
                    found_in_model = True
                    break

            if not found_in_model:
                # 添加新反应
                rxn_obj = cobra.Reaction(r['bigg_id'][0])
                rxn_obj.lower_bound = 0 if direction == 0 else -1000
                rxn_obj.upper_bound = 1000 if direction == 0 else 0
                rxn_obj.add_metabolites(add_cpd2rxn(model, rxn, total_rxn_dict, cpd_dict))
                model.add_reactions([rxn_obj])
                log.logger.info(f"添加反应：{rxn} {r['bigg_id'][0]}")
        else:
            # 添加新反应（无 BiGG ID）
            rxn_obj = cobra.Reaction(rxn)
            rxn_obj.lower_bound = 0 if direction == 0 else -1000
            rxn_obj.upper_bound = 1000 if direction == 0 else 0
            rxn_obj.add_metabolites(add_cpd2rxn(model, rxn, total_rxn_dict, cpd_dict))
            model.add_reactions([rxn_obj])
            log.logger.info(f"添加反应：{rxn}")

    # 设置目标产物交换反应
    if product not in cpd_dict:
        log.logger.error(f"未找到 {product} 对应的化合物信息。")
        return -999

    product_info = cpd_dict[product]
    if product_info['bigg_id'] is None:
        log.logger.error(f"目标产物 {product} 无对应的 BiGG ID。")
        return -999

    target_bigg_id = product_info['bigg_id'][0]
    try:
        target_metabolite = model.metabolites.get_by_id(target_bigg_id + "_c")
    except KeyError:
        log.logger.error(f"未在模型中找到目标产物：{target_bigg_id}")
        return -999

    # 创建交换反应
    exchange_rxn = cobra.Reaction(f"EX_{target_metabolite.id}")
    exchange_rxn.lower_bound = -1000  # 允许产物排出
    exchange_rxn.upper_bound = 1000
    exchange_rxn.add_metabolites({target_metabolite: -1})
    model.add_reactions([exchange_rxn])

    # 优化模型，计算理论产量
    model.objective = exchange_rxn
    solution = model.optimize()
    max_biomass = solution.objective_value

    if max_biomass < min_growth:
        log.logger.warning("路径不可行（生物量低于阈值）！")
        return -999
    else:
        if model.reactions.get_by_id('EX_glc__D_e').flux == 0:
            log.logger.warning("葡萄糖吸收通量为 0，理论产率 = 0")
            return 0
        else:
            maximum_yield = max_biomass / (-1 * model.reactions.get_by_id('EX_glc__D_e').flux)
            log.logger.info(f'理论产率 = {maximum_yield} mmol-产物/mmol-葡萄糖')
            return maximum_yield