import os
import sys
import csv
import numpy as np
import time
import argparse
import random
import shutil
import datetime

# 保留现有的导入内容和其他函数
from mooseeker.model.BiOperators import *
from mooseeker.utils import *
from mooseeker.fit_funcs import get_gibbs_mlp, get_yield,get_network_totalReactions,get_toxicity,get_network_branch_fitness
from mooseeker.utils import get_config
from mooseeker.utils import SingleLinkList
from mooseeker.utils import SingleReaction
from mooseeker.log import Logger
def parser():
    params = {
        "task_id": "013",
        "algorithm": "single",
        "NrMax": 15,
        "ob_substrate": "C00103",
        "ob_product": "C00631",
        "abundant": ["C00001", "C00002", "C00080"],
        "database": ["KEGG", "BiGG"],
        "host": ["ecoli"],
    }
    args = argparse.Namespace(**params)
    args.save_path = ".."
    args.project = '_'.join((args.task_id, args.ob_substrate, args.ob_product))
    args.save_path = '/'.join((args.save_path, args.project)) + '/'
    if not os.path.exists(args.save_path):
        os.mkdir(args.save_path)
    return args

args = parser()
cfg = get_config()
unique_paths = set()
# log = Logger(cfg['file_path']['log_dir'] + 'SFLA.log', fmt='%(message)s')

pop_size = 15
max_iter = 20
CROSS_RATE = 0.8
MUTATION_RATE = 0.4  # 根据路径特殊性变异概率大一点
weight1 = 0.3
weight2 = 0.2
weight3 = 0.3
weight4 = 0.2
# 获取当前时间戳
timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
project = '_'.join((
    'SFLA', 
    str(args.ob_product), 
    str(weight1), 
    str(weight2), 
    str(weight3), 
    str(weight4), 
    timestamp
))
result_dir = cfg['file_path']['result_dir'] + project + '/'

if not os.path.exists(result_dir):
    os.mkdir(result_dir)

log_file_path = result_dir + 'SFLA.log'
log = Logger(log_file_path, fmt='%(message)s')

if os.path.exists(result_dir + 'all_sfla.txt'):
    os.remove(result_dir + 'all_sfla.txt')

def evaluate(chrom):
    weight1 = 0.3
    weight2 = 0.2
    weight3 = 0.3
    weight4 = 0.2
    product = 'C00631'

    length = chrom.chrom.length()
    gibbs = get_gibbs_mlp(cfg,chrom.chrom.get_rxn_list())
    total_reactions = get_network_totalReactions(cfg, chrom, log)
    net_tox = get_toxicity(cfg, chrom, log, 0)
    net_count = get_network_branch_fitness(cfg, chrom, log)
    network_yield = get_yield(cfg, chrom, product, log, 0)

    fit = float(weight1) * net_count + float(weight2) * net_tox - float(weight3) * network_yield + float(weight4) * gibbs


    log.logger.info(f'Length: {length}\nNet_count: {net_count}\nNet_Tox: {net_tox}\nNetwork_yield: {network_yield}\nGibbs: {gibbs}\nFitness: {fit}\n')

    csv_file_path = result_dir + 'evaluation_data.csv'
    if not os.path.exists(csv_file_path):
        with open(csv_file_path, 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["Fitness", "Total Reactions", "Length", "Net_count", "Net_Tox", "Network_yield", "Gibbs"])

    with open(csv_file_path, 'a', newline='') as file:
        writer = csv.writer(file)
        writer.writerow([fit, total_reactions, length, net_count, net_tox, network_yield, gibbs])

    with open(result_dir + 'all_sfla.txt', 'a') as f:
        f.write(f"[{length}, {fit}]\n")
        reactions_list = chrom.chrom.get_rxn_list()
        reactions_str = ','.join(reactions_list)
        f.write(f"{reactions_str}\n")

        travel_list = chrom.chrom.travel()
        compounds = [line[1] for line in travel_list]
        if travel_list:
            compounds.append(travel_list[-1][2])
        compounds_str = ','.join(compounds)
        f.write(f"{compounds_str}\n")

        f.write('==' * 10 + '\n')

    unique_path_str = f"Reactions: {reactions_str} | Compounds: {compounds_str}"
    unique_paths.add(unique_path_str)

    chrom.set_fit(fit)
    return fit, length


# 整合变异函数
def mut(Chrome):
    chrom = copy.deepcopy(Chrome)
    cond = False  # 用来表示是否变异成功

    while True:
        m = random.randint(0, chrom.chrom.length() - 2)
        cur = chrom.chrom._head
        for i in range(m):
            cur = cur.next

        s = cur.P

        # 生成突变个体
        mut_chrom = Chromesome(args, cfg, log, initial_substrate=s)
        if mut_chrom.get_a_chrom():
            cur.next = mut_chrom.chrom._head

        chrom = trim_chrom(chrom)

        if chrom.chrom.length() <= 20:
            cond = True  # 变异成功
            break

    chrom.update_pool_after_change()
    return cond, chrom  # 返回是否成功和新的个体


# 修剪函数
def trim_chrom(input_chrom):
    chrom = copy.deepcopy(input_chrom)
    head = chrom.chrom._head
    if chrom.chrom.length() > 2:
        dummy = SingleReaction(S=None, P=None, reaction={})
        dummy.next = head
        pre = dummy
        left = pre.next
        while left is not None:
            right = left.next
            while right is not None:
                if left.S == right.S:
                    pre.next = right
                    if pre.S is None:
                        chrom.chrom._head = pre.next
                    left = pre.next
                right = right.next
            pre = pre.next
            left = left.next
            pre.next = left
        return chrom
    else:
        return input_chrom

def get_chrom_path(chrom):
    """从SingleLinkList提取路径"""
    current = chrom._head  # 从链表的头节点开始
    path = []
    while current is not None:
        path.append(current)  # 将节点添加到路径
        current = current.next  # 移动到下一个节点
    return path

def local_search(frog, best_frog, worst_frog):
    # 计算步长，步长决定修改路径的程度
    step_size = np.random.uniform(0, 1) * (worst_frog.chrom.length() - frog.chrom.length())

    # 获取青蛙个体的路径
    frog_path = get_chrom_path(frog.chrom)
    best_frog_path = get_chrom_path(best_frog.chrom)

    new_frog = copy.deepcopy(frog)

    # 更新路径：调整某些节点以使其向最优路径靠近
    for i in range(int(step_size)):  # 根据步长调整路径中的某些节点
        if i < len(frog_path) and i < len(best_frog_path):
            frog_path[i] = best_frog_path[i]  # 向最优路径靠拢

    # 返回调整后的青蛙个体
    return new_frog


def global_mix(frogs, best_frog=None):
    # 如果有最优个体，先将其加入新的种群中
    if best_frog:
        frogs = [best_frog] + frogs
    
    random.shuffle(frogs)
    return frogs

# 初始化种群
def initialize_population(pop_size, cfg, log):
    chrom_list = []  # 改为普通列表
    count = 0
    log.logger.info("-----Sampling-----")
    while count < pop_size:
        chrom = Chromesome(args, cfg, log)
        if chrom.get_a_chrom():
            log.logger.info('--This is the %dth chromsome---' % (count + 1))
            chrom_list.append(chrom)  # 直接添加 chrom 对象
            count += 1
    log.logger.info("-----Sampling Done-----")
    return chrom_list

import numpy as np

def select(pop, fitness):
    # 计算适应度总和
    total_fitness = sum(fitness)
    
    # 确保适应度值为正，避免零或负数
    probabilities = [1 / (fit if fit > 0 else 1e-10) for fit in fitness]
    
    # 计算总概率并归一化
    total_probability = sum(probabilities)
    if total_probability == 0:
        # 如果概率和为零，采用均匀分布
        probabilities = [1 / len(fitness)] * len(fitness)
    else:
        probabilities = [prob / total_probability for prob in probabilities]

    # 使用numpy进行随机选择
    selected_index = np.random.choice(np.arange(len(fitness)), size=pop_size, replace=True, p=probabilities)
    
    # 根据选择的索引更新种群
    pop = [pop[idx] for idx in selected_index]
    return pop

# 精英选择策略
def elite_select(all_pop, all_fitness, all_length):
    best_index = np.argmin(all_fitness)  # 选择适应度最小的个体（假设适应度越小越好）
    elite_individual = all_pop[best_index]
    elite_fitness = all_fitness[best_index]
    elite_length = all_length[best_index]

    # 保证最优个体一定进入下一代
    selected_pop = [elite_individual]
    selected_fitness = [elite_fitness]
    selected_length = [elite_length]

    # 确保适应度值为正，避免零或负数
    probabilities = [1 / (fit if fit > 0 else 1e-10) for fit in all_fitness]
    
    # 计算总概率并归一化
    total_probability = sum(probabilities)
    if total_probability == 0:
        # 如果概率和为零，采用均匀分布
        probabilities = [1 / len(all_fitness)] * len(all_fitness)
    else:
        probabilities = [prob / total_probability for prob in probabilities]

    # 选择剩余个体
    selected_index = np.random.choice(np.arange(len(all_fitness)), size=pop_size - 1, replace=True, p=probabilities)

    # 加入选中的个体
    selected_pop.extend([all_pop[idx] for idx in selected_index])
    selected_fitness.extend([all_fitness[idx] for idx in selected_index])
    selected_length.extend([all_length[idx] for idx in selected_index])

    return selected_pop, selected_fitness, selected_length




# 保存每代个体信息
def save_every_generation(index, pop, fitness, chrom_length):
    with open(result_dir + 'every_generation_childs.txt', 'a') as f:
        f.write(f"===== Generation {index} =====\n")
        for i in range(len(fitness)):
            length = chrom_length[i]
            fit = fitness[i]
            f.write(str([length, fit]) + '\n')
            chrom = pop[i]
            travel_list = chrom.chrom.travel()
            for line in travel_list:
                for l in line:
                    f.write(l)
                    f.write('\t')
                f.write('\n')


def sfla_algorithm():
    # 初始化种群
    chrom_list = initialize_population(pop_size, cfg, log)
    fitness = []
    chrom_length = []

    # 计算初始个体的适应度值，并修剪个体
    for chrom in chrom_list:
        # **先进行剪枝操作**
        chrom = trim_chrom(chrom)

        # **保存剪枝后的化合物路径**
        fit, length = evaluate(chrom)
        fitness.append(fit)
        chrom_length.append(length)

    # 保存初始最优个体信息
    best_chrom_idx = np.argmin(fitness)
    best_fitness = min(fitness)
    with open(result_dir + "best.csv", "w", newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Generation", "Best Fitness"])
        csvwriter.writerow([0, best_fitness])

    # 迭代过程
    for i in range(max_iter):
        log.logger.info(f"----- Generation {i + 1} -----")

        # 对种群进行局部搜索
        for group in range(0, len(chrom_list), 5):  # 分成子群体
            group_frogs = chrom_list[group:group + 5]
            best_frog = min(group_frogs, key=lambda frog: frog.fit)
            worst_frog = max(group_frogs, key=lambda frog: frog.fit)

            for frog in group_frogs:
                # 局部搜索
                new_frog = local_search(frog, best_frog, worst_frog)
                
                # **再次进行剪枝操作**：确保局部搜索后的个体经过剪枝
                new_frog = trim_chrom(new_frog)
                
                # **评估和保存剪枝后的路径信息**
                new_frog_fit, _ = evaluate(new_frog)
                if new_frog_fit < worst_frog.fit:
                    # 找到worst_frog在chrom_list中的位置并替换
                    worst_frog_idx = next((i for i, chrom in enumerate(chrom_list) if id(chrom) == id(worst_frog)), None)
                    if worst_frog_idx is not None:
                        chrom_list[worst_frog_idx] = new_frog
                    else:
                        log.logger.info("Warning: worst_frog not found in chrom_list")

                # 变异操作
                if np.random.rand() < MUTATION_RATE:
                    cond, mut_cross = mut(frog)
                    frog_idx = next((i for i, chrom in enumerate(chrom_list) if id(chrom) == id(frog)), None)
                    if frog_idx is not None:
                        chrom_list[frog_idx] = mut_cross
                    else:
                        log.logger.info("Warning: frog not found in chrom_list")

        # 全局混合
        chrom_list = global_mix(chrom_list)

        # 更新适应度并保存最优解
        fitness = [evaluate(chrom)[0] for chrom in chrom_list]
        best_fitness = min(fitness)
        with open(result_dir + "best.csv", "a", newline='') as csvfile:
            csvwriter = csv.writer(csvfile)
            csvwriter.writerow([i + 1, best_fitness])

        # 保存每代最优个体详细信息
        best_chrom_idx = np.argmin(fitness)
        best_chrom = chrom_list[best_chrom_idx]
        save_every_generation(i + 1, chrom_list, fitness, chrom_length)
        with open(result_dir + 'best_childs.txt', 'a') as f:
            f.write(f"Generation {i+1}\n")
            # 路径长度，网络中总反应个数，适应度值（其实现在就是网络中总反应个数）
            length = best_chrom.chrom.length()
            total_reactions = best_chrom.total_reactions
            fit = best_chrom.fit

            f.write(f"Length: {length}, Total Reactions: {total_reactions}, Fitness: {fit}\n")
            
            # 写入路径中的反应及化合物
            travel_list = best_chrom.chrom.travel()
            f.write("Reactions and Compounds:\n")
            for line in travel_list:
                f.write("\t".join(line) + "\n")
            
            # 写入最好的酶组合
            f.write(f"Best Enzyme Pool: {best_chrom.best_enzyme_pool}\n")
            # 写入当前酶组合产生的总反应
            f.write(f"Best Total Reactions: {best_chrom.best_total_rxn_list}\n")
            f.write("=" * 20 + "\n")

        # 每代结束后选择最优个体进入下一代
        chrom_list, fitness, chrom_length = elite_select(chrom_list, fitness, chrom_length)

    # 算法结束后输出所有不重复路径到 final.txt
    with open(result_dir + 'final.txt', 'w') as f:
        for path_str in unique_paths:
            f.write(path_str + '\n')

start_time = time.time()  
sfla_algorithm()
end_time = time.time()
log.logger.info('--Execution Time: %f h---' % ((end_time - start_time) / 3600))




