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

# 吉布斯自由能和理论产量都要考虑

moo_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
print("moo_dir:", moo_dir)

sys.path.append(moo_dir)  # 添加project路径
sys.path.append(moo_dir + "/mooseeker/model")  # 添加project路径下的mooseeker module路径
sys.path.append(moo_dir + "/mooseeker/utils")  # 添加project路径下的mooseeker module路径

# 这个是我自己加的 不然运行nohup命令的时候运行不了
sys.path.append("..")

from mooseeker.model.BiOperators import *
from mooseeker.utils import *

def parser():
    params = {
        "task_id": "013",
        "algorithm": "single",
        "NrMax": 15,
        "ob_substrate": "C00082",
        "ob_product": "C00755",
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
# log = Logger(cfg['file_path']['log_dir'] + 'PSO.log', fmt='%(message)s')


# 获取当前时间戳
timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
project = '_'.join(('PSO', args.ob_product, timestamp))
result_dir = cfg['file_path']['result_dir'] + project + '/'

if not os.path.exists(result_dir):
    os.mkdir(result_dir)

log_file_path = result_dir + 'PSO.log'
log = Logger(log_file_path, fmt='%(message)s')

if os.path.exists(result_dir + 'all_pso.txt'):
    os.remove(result_dir + 'all_PSO.txt')

def evaluate(chrom):
    weight1 = 0.3
    weight2 = 0.2
    weight3 = 0.3
    weight4 = 0.2
    product = 'C00755'

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

    with open(result_dir + 'all_pso.txt', 'a') as f:
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
    return fit

# 修剪函数
def trim_chrom(input_chrom):
    chrom = copy.deepcopy(input_chrom)
    if chrom.chrom is None or chrom.chrom._head is None:  # 检查是否有效
        log.logger.warning("Invalid chromosome, skipping trimming.")
        return input_chrom  # 如果无效，返回原始染色体

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


def cross(chromA, chromB):
    """cross chromA and chromB
    Args:
        chromA (class Chromsomes):
        chromB (class Chromsomes):
    Returns:
        newA (class Chromsomes):
        newB (class Chromsomes):
    """
    newA, newB = copy.deepcopy(chromA), copy.deepcopy(chromB)
    curA = newA.chrom._head
    while (curA != None and curA.next != None):
        curB = newB.chrom._head
        while (curB != None and curB.next != None):
            if (curA.P == curB.P):
                # find the same point and crossover
                temp = curA.next
                curA.next = curB.next
                curB.next = temp

                newA.update_pool_after_change()
                newB.update_pool_after_change()

                log.logger.info("-----Crossover Successfully!-----")
                log.logger.info("Chromesomes before crossover is :")
                log.logger.info(chromA.chrom.travel())
                log.logger.info("Chromesomes after crossover is :")
                log.logger.info(newA.chrom.travel())
                log.logger.info('--'*4)

                return newA, newB
            else:
                curB = curB.next
        curA = curA.next

def has_same_node(chromA, chromB):
    """Judge chromA and chromB has the same point to crossover
    Args:
        chromA (class AChromsomes): _description_
        chromB (class AChromsomes): _description_

    Returns:
        bool: true or false
    """
    curA = chromA.chrom._head
    while (curA != None and curA.next != None):
        curB = chromB.chrom._head
        while (curB != None and curB.next != None):
            if (curA.P == curB.P):
                if (curA == chromA.chrom._head and curB == chromB.chrom._head):
                    curB = curB.next
                    continue
                elif (curA.next.next == None and curB.next.next == None):
                    curA = curA.next
                    continue
                else:
                    return True
            curB = curB.next
        curA = curA.next

    return False

def perform_crossover(chrom1, chrom2):
    chrom1 = copy.deepcopy(chrom1)  # 确保深拷贝
    chrom2 = copy.deepcopy(chrom2)

    curA = chrom1.chrom._head
    while curA is not None and curA.next is not None:
        curB = chrom2.chrom._head
        while curB is not None and curB.next is not None:
            if curA.P == curB.P:  # 找到交叉点
                temp = curA.next
                curA.next = curB.next
                curB.next = temp
                chrom1.update_pool_after_change()  # 更新染色体
                chrom2.update_pool_after_change()  # 更新染色体

                log.logger.info("-----Crossover Successfully!-----")
                log.logger.info("Chromosomes before crossover is :")
                log.logger.info(chrom1.chrom.travel())  # 打印交叉前染色体路径
                log.logger.info("Chromosomes after crossover is :")
                log.logger.info(chrom1.chrom.travel())  # 打印交叉后染色体路径
                return chrom1, chrom2  # 返回交叉后的染色体

            curB = curB.next
        curA = curA.next

    return chrom1, chrom2  # 如果没有交叉点，返回原始染色体
def perform_crossover(chrom1, chrom2):
    chrom1 = copy.deepcopy(chrom1)  # 确保深拷贝
    chrom2 = copy.deepcopy(chrom2)

    curA = chrom1.chrom._head
    while curA is not None and curA.next is not None:
        curB = chrom2.chrom._head
        while curB is not None and curB.next is not None:
            if curA.P == curB.P:  # 找到交叉点
                temp = curA.next
                curA.next = curB.next
                curB.next = temp
                chrom1.update_pool_after_change()  # 更新染色体
                chrom2.update_pool_after_change()  # 更新染色体

                log.logger.info("-----Crossover Successfully!-----")
                log.logger.info("Chromosomes before crossover is :")
                log.logger.info(chrom1.chrom.travel())  # 打印交叉前染色体路径
                log.logger.info("Chromosomes after crossover is :")
                log.logger.info(chrom1.chrom.travel())  # 打印交叉后染色体路径
                return chrom1, chrom2  # 返回交叉后的染色体

            curB = curB.next
        curA = curA.next

    return chrom1, chrom2  # 如果没有交叉点，返回原始染色体
def perform_crossover(chrom1, chrom2):
    chrom1 = copy.deepcopy(chrom1)  # 确保深拷贝
    chrom2 = copy.deepcopy(chrom2)

    curA = chrom1.chrom._head
    while curA is not None and curA.next is not None:
        curB = chrom2.chrom._head
        while curB is not None and curB.next is not None:
            if curA.P == curB.P:  # 找到交叉点
                temp = curA.next
                curA.next = curB.next
                curB.next = temp
                chrom1.update_pool_after_change()  # 更新染色体
                chrom2.update_pool_after_change()  # 更新染色体

                log.logger.info("-----Crossover Successfully!-----")
                log.logger.info("Chromosomes before crossover is :")
                log.logger.info(chrom1.chrom.travel())  # 打印交叉前染色体路径
                log.logger.info("Chromosomes after crossover is :")
                log.logger.info(chrom1.chrom.travel())  # 打印交叉后染色体路径
                return chrom1, chrom2  # 返回交叉后的染色体

            curB = curB.next
        curA = curA.next

    return chrom1, chrom2  # 如果没有交叉点，返回原始染色体

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

# 精英选择策略
def elite_select(all_pop, all_fitness, all_length):
    best_index = np.argmin(all_fitness)
    elite_individual = all_pop[best_index]
    elite_fitness = all_fitness[best_index]
    elite_length = all_length[best_index]

    total_fitness = sum(all_fitness)
    probabilities = [1 / fit for fit in all_fitness]
    total_probability = sum(probabilities)
    normalized_probabilities = [prob / total_probability for prob in probabilities]

    selected_index = np.random.choice(np.arange(len(all_fitness)), size=pop_size - 1, replace=True, p=normalized_probabilities)

    selected_pop = [all_pop[idx] for idx in selected_index]
    selected_fitness = [all_fitness[idx] for idx in selected_index]
    selected_length = [all_length[idx] for idx in selected_index]

    pop = np.concatenate(([elite_individual], selected_pop))
    fitness = np.concatenate(([elite_fitness], selected_fitness))
    length = np.concatenate(([elite_length], selected_length))

    return pop, fitness, length

# 初始化参数
pop_size = 15
max_iter = 20
count = 0
chrom_list = []

# 开始计时
start_time = time.time()
log.logger.info("-----Sampling-----")

# 初始化染色体列表
chrom_list = []
while count < pop_size:
    chrom = Chromesome(args, cfg, log)
    if chrom.get_a_chrom():
        log.logger.info('--This is the %dth chromsome---' % (count + 1))
        chrom_list.append(chrom)  # 直接使用append来存储Chromesome对象
        count += 1

log.logger.info("-----Sampling Done-----")

# 保存最优结果
with open(result_dir + "best.csv", "w", newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["Generation", "Best Fitness"])

# 初始化适应度
fitness = []
for idx, chrom in enumerate(chrom_list):
    chrom = trim_chrom(chrom)
    fit = evaluate(chrom)
    fitness.append(fit)

# 初始化个人最优和全局最优
personal_best = chrom_list.copy()
personal_best_fit = fitness.copy()
global_best = chrom_list[np.argmin(fitness)]
global_best_fit = min(fitness)

# PSO迭代
for i in range(max_iter):
    log.logger.info(f"--- Iteration {i+1} ---")
    
    # 计算当前代每个染色体的长度
    chrom_length = [chrom.chrom.length() for chrom in chrom_list]

    for idx, chrom in enumerate(chrom_list):
        # 使用 copy 模块复制对象
        original_chrom = copy.copy(chrom)  # 存储原始染色体
        original_fitness = fitness[idx]  # 存储原始适应度

        log.logger.info(f'=====> Mutation of the {i+1}th iteration, {idx+1}th particle <=====') 

        # 变异操作
        cond, mut_chrom = mut(chrom)
        if cond:
            chrom = mut_chrom  # 更新染色体
        else:
            chrom = mut_chrom  # 即使没有变异，也更新

        # 变异后的染色体修剪和适应度评估
        chrom = trim_chrom(chrom)
        fit_mut = evaluate(chrom)

        if fit_mut >= fitness[idx]:
            # 如果变异后的适应度较好，进行局部交叉
            log.logger.info("----- Crossover with local -----")
            chrom1, chrom2 = perform_crossover(chrom, personal_best[idx])  # 调用交叉，返回两个染色体
            chrom1 = trim_chrom(chrom1)  # 修剪交叉后的染色体
            chrom2 = trim_chrom(chrom2)  # 修剪交叉后的染色体
            fit_cross_loc = evaluate(chrom1)  # 评估第一个染色体

            if fit_cross_loc >= fitness[idx]:
                # 如果局部交叉后的适应度较好，进行全局交叉
                log.logger.info("----- Crossover with global -----")
                chrom1, chrom2 = perform_crossover(chrom1, global_best)  # 全局交叉
                chrom1 = trim_chrom(chrom1)  # 修剪交叉后的染色体
                chrom2 = trim_chrom(chrom2)  # 修剪交叉后的染色体
                fit_cross_glo = evaluate(chrom1)  # 评估全局交叉后的染色体

                if fit_cross_glo < original_fitness:
                    fitness[idx] = fit_cross_glo
                    chrom_list[idx] = chrom1
                else:
                    # 否则恢复原始染色体
                    chrom_list[idx] = original_chrom
                    fitness[idx] = original_fitness
            else:
                # 如果局部交叉后的适应度较好，更新
                fitness[idx] = fit_cross_loc
                chrom_list[idx] = chrom1
        else:
            # 如果变异后的适应度较差，恢复原始染色体
            fitness[idx] = fit_mut
            chrom_list[idx] = chrom

        # 更新个人最优
        if fitness[idx] < personal_best_fit[idx]:
            personal_best_fit[idx] = fitness[idx]
            personal_best[idx] = chrom

        # 更新全局最优
        if fitness[idx] < global_best_fit:
            global_best_fit = fitness[idx]
            global_best = chrom

    # 每代记录最优结果
    with open(result_dir + "best.csv", "a", newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow([i + 1, min(fitness)])
    
    # 保存每代最优个体详细信息
    best_chrom_idx = np.argmin(fitness)
    best_chrom = chrom_list[best_chrom_idx]
    save_every_generation(i + 1, chrom_list, fitness, chrom_length)

    # 保存每代最优个体的详细路径信息
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

# 结束计时并输出
end_time = time.time()
log.logger.info('--Execution Time: %f hours---' % ((end_time - start_time) / 3600))
