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


#sys.path.append('.')
# 获取当前工作目录
cwd = os.getcwd()
# 获取上级目录
parent_dir = os.path.abspath(os.path.join(cwd, os.pardir))
sys.path.append(parent_dir)

moo_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
print("moo_dir:",moo_dir)

sys.path.append(moo_dir) #添加project路径
sys.path.append(moo_dir + "/mooseeker/model") #添加project路径下的mooseeker module路径
sys.path.append(moo_dir + "/mooseeker/utils") #添加project路径下的mooseeker module路径

# 这个是我自己加的 不然运行nohup命令的时候运行不了
sys.path.append("..")

def parser():
    
    params = {
        "task_id": "013",
        "algorithm": "single",
        # "pop_size": 10,  # 原本是10
        # "gen": 20,  # 原本是20
        "NrMax": 15,
        "ob_substrate": "C00082",
        "ob_product": "C00755",
        "abundant": ["C00001", "C00002", "C00080"],
        "database": ["KEGG", "BiGG"],
        "host": ["ecoli"],
        # "eva_func": ["yield"]
        # "eva_func": ["length", "gibbs", "yield"]
    }

    args = argparse.Namespace(**params)
    args.save_path = ".."
    args.project = '_'.join((args.task_id, args.ob_substrate, args.ob_product))
    args.save_path = '/'.join((args.save_path, args.project)) + '/'
    print("main_single:args.save_path")
    print(args.save_path)
    if not os.path.exists(args.save_path): os.mkdir(args.save_path)
    return args

args = parser()
cfg = get_config()
unique_paths = set()
# log = Logger(cfg['file_path']['log_dir'] + 'GA.log', fmt='%(message)s')

pop_size=15
max_iter = 20
CROSS_RATE = 0.8
MUTATION_RATE = 0.4 # 我觉得根据这个路径的特殊性 变异的概率可以大一点


# 获取当前时间戳
timestamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
project = '_'.join(('GA', args.ob_product, timestamp))
result_dir = cfg['file_path']['result_dir'] + project + '/'

if not os.path.exists(result_dir):
    os.mkdir(result_dir)

log_file_path = result_dir + 'GA.log'
log = Logger(log_file_path, fmt='%(message)s')

if os.path.exists(result_dir + 'all_ga.txt'):
    os.remove(result_dir + 'all_ga.txt')

def evaluate(chrom):
    # Ensure chrom is the actual chromosome object, not a numpy ndarray
    chrom = chrom[0]  # This accesses the first element of the ndarray

    weight1 = 0.3
    weight2 = 0.2
    weight3 = 0.3
    weight4 = 0.2
    product = 'C00755'

    length = chrom.chrom.length()  # Now 'chrom' should be the Chromesome object
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

    with open(result_dir + 'all_ga.txt', 'a') as f:
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

def mut(Chrome):
    chrom = copy.deepcopy(Chrome)

    while True:
        m = random.randint(0, chrom[0].chrom.length()-2)
        cur = chrom[0].chrom._head
        for i in range(m):
            cur = cur.next

        s = cur.P

        # generate mut_chrom
        mut_chrom = np.full((1, 1), None, dtype=object)
        while (1):
            temp_chrom = Chromesome(args, cfg, log, initial_substrate=s)
            if temp_chrom.get_a_chrom():
                mut_chrom[0] = temp_chrom
                break

        # connect chrom and mut_chrom
        cur = chrom[0].chrom._head
        mut_cur = mut_chrom[0][0].chrom._head

        for i in range(m):
            cur = cur.next

        cur.next = mut_cur

        # trim new chrom
        chrom[0].trim_chrom(chrom[0].chrom)

        # 检查修剪后的路径长度
        if chrom[0].chrom.length() <= 20:
            break
    # update chrom's pool
    chrom[0].update_pool_after_change()
    if chrom[0].chrom.get_rxn_list() != Chrome[0].chrom.get_rxn_list():

        log.logger.info("-----Mutation Successfully!-----")
        log.logger.info("Chromesomes before mutation is :")
        log.logger.info(Chrome[0].chrom.travel())
        log.logger.info("Chromesomes after mutation is :")
        log.logger.info(chrom[0].chrom.travel())
        log.logger.info('--' * 10)

        return [True, chrom]
    else:
        # problem.log.logger.info(chrom.travel())
        log.logger.info("-----Mutation Failed!-----")
        return [False, chrom]

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
    # newA, newB = chromA, chromB
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
                log.logger.info('--' * 5)
                log.logger.info(chromB.chrom.travel())
                log.logger.info("Chromesomes after crossover is :")
                log.logger.info(newA.chrom.travel())
                log.logger.info('--' * 5)
                log.logger.info(newB.chrom.travel())
                log.logger.info('--'*10)

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
    # 遍历到链表的前一个节点 除掉最后一个节点
    while (curA != None and curA.next != None):
        curB = chromB.chrom._head
        while (curB != None and curB.next != None):
            if (curA.P == curB.P):
                # 如果都是头节点 则不满足有相同节点的条件
                if (curA == chromA.chrom._head and curB == chromB.chrom._head):
                    curB = curB.next
                    continue

                # 如果都是倒数第二节点，也不满足有相同节点的条件
                elif (curA.next.next == None and curB.next.next == None):
                    curB = curB.next
                    continue
                else:
                    return True
            else:
                curB = curB.next
        curA = curA.next
    return False

def perform_crossover(chrom1, chrom2):
    chrom1 = np.array(chrom1) # (1,)
    chrom2 = np.array(chrom2)
    stacked = np.vstack((chrom1, chrom2))
    # 垂直方向拼接为一个新的数组 chrom1=stacked[0][0] chrom2=stacked[1][0]
    n_parents_, n_matings = stacked.shape # 2,1
    output = np.full_like(stacked, None, dtype=object)
    # 生成一个形状和stacked一样的array  (2,1)
    # 初始为none array([[None],[None]], dtype=object)
    for i in range(n_matings):
        # # 如果有则交叉
        if has_same_node(stacked[0, i], stacked[1, i]):
            output[0, i], output[1, i] = cross(stacked[0, i], stacked[1, i])
        else:
            # 没有就返回原始的chrom1
            output[0, i], output[1, i] = stacked[0, i], stacked[1, i]
            log.logger.info("-----Has no same node, Crossover failed!-----")

    return output # 返回交叉后的个体 最后返回的是交叉后的两个个体

def trim_chrom(input_chrom):
    chrom = copy.deepcopy(input_chrom)

    head = chrom[0].chrom._head
    # 修剪情况: a->b b->c c->a a->d ===> a->d
    if (chrom[0].chrom != None and chrom[0].chrom.length() > 2):
        dummy = SingleReaction(S=None, P=None, reaction={})
        dummy.next = head
        pre = dummy
        left = pre.next
        while(left != None):
            right = left.next
            while(right != None):
                if left.S == right.S:
                    pre.next = right
                    if (pre.S == None):
                        chrom[0].chrom._head = pre.next
                    left = pre.next
                    # break
                right = right.next
            pre = pre.next
            left = left.next
            pre.next = left  # 添加的一行
        return chrom
    else:
        return input_chrom

## 我觉得师姐这个缺少了 选择 的一步
def select(pop, fitness):
    # 计算适应度总和
    total_fitness = sum(fitness)
    # 计算每个个体被选择的概率，适应度值越小，概率越大
    probabilities = [1 / fit for fit in fitness]
    total_probability = sum(probabilities)
    # 归一化概率，使概率之和为1
    normalized_probabilities = [prob / total_probability for prob in probabilities]
    
    # 使用numpy进行随机选择，确保最优个体一定会被选中
    best_idx = np.argmin(fitness)
    # 保证最优个体进入下一代
    selected_index = np.random.choice(np.arange(len(fitness)), size=pop_size - 1, replace=True, p=normalized_probabilities)
    pop = [pop[best_idx]] + [pop[idx] for idx in selected_index]
    return pop


def save_every_generation(index, pop, fitness, chrom_length):
    with open(result_dir + 'every_generation_childs.txt', 'a') as f:
        f.write(f"===== Generation {index} =====\n")
        for i in range(len(fitness)):
            length = chrom_length[i]
            fit = fitness[i]
            f.write(str([length, fit]) + '\n')
            chrom = pop[i][0]  # 获取 Chromesome 对象
            travel_list = chrom.chrom.travel()  # 访问 Chromesome 对象的 chrom 属性
            for line in travel_list:
                for l in line:
                    f.write(l)
                    f.write('\t')
                f.write('\n')


count = 0
# chrom_list = []
start_time = time.time()  # 开始计时

log.logger.info("-----Sampling-----")
# 生成初始的个体
chrom_list = np.full((pop_size, 1), None, dtype=object)
while count < pop_size:
    chrom = Chromesome(args, cfg, log)
    # 这里在生成chrom时，会调用目标底物与目标产物
    if chrom.get_a_chrom():
        log.logger.info('--This is the %dth chromsome---' % (count+1))
        # chrom_list.append(chrom)
        chrom_list[count, 0] = chrom
        count += 1
log.logger.info("-----Sampling Done-----")

# 写入 CSV 文件，记录每代的最优适应度
with open(result_dir + "best.csv", "w", newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["Generation", "Best Fitness"])

# 计算初始个体的适应度值
fitness = []
for idx, chrom in enumerate(chrom_list):
    chrom = trim_chrom(chrom)
    fit = evaluate(chrom)
    fitness.append(fit)

chrom_next = chrom_list

# 开始迭代
for i in range(max_iter):
    # 获取当前代的最优个体和适应度
    best_chrom_idx = np.argmin(fitness)
    best_chrom = chrom_list[best_chrom_idx]
    best_fitness = fitness[best_chrom_idx]

    log.logger.info(f"-----Generation {i+1}: Best Fitness So Far: {best_fitness}-----")
    
    # 从chrom_list中选择两个不同的索引值(除最优个体）
    remaining_indices = [idx for idx in range(len(chrom_list)) if idx != best_chrom_idx]
    
    # 随机选择进行交叉的两个个体
    if len(remaining_indices) > 1:
        indices = np.random.choice(remaining_indices, 2, replace=False)
        index1 = indices[0]
        index2 = indices[1]

        # 交叉
        if np.random.rand() < CROSS_RATE:
            log.logger.info("----- start Crossover -----")
            # 如果概率在交叉的概率内才交叉
            Y = perform_crossover(chrom_list[index1], chrom_list[index2])
            chrom_list[index1] = Y[0]  # 替换为交叉后的个体
            chrom_list[index2] = Y[1]  # 替换为交叉后的个体
            log.logger.info("-----Crossover Done-----")

    # 随机选择一个个体进行变异
    mutate_idx = np.random.choice(remaining_indices)

    if np.random.rand() < MUTATION_RATE:
        cond, mut_cross = mut(chrom_list[mutate_idx])
        if cond:
            chrom_list[mutate_idx] = mut_cross
        else:
            # 如果变异失败，重新生成一个新的个体替代
            newchrom = Chromesome(args, cfg, log)
            if newchrom.get_a_chrom():
                chrom_list[mutate_idx][0] = newchrom
        log.logger.info('=====> Mutation of the %dth iter %dth particle<=====' % (i+1, mutate_idx+1))

    # 保留最优个体到下一代
    chrom_list[0] = best_chrom  # 将上一代的最优个体直接作为下一代的第一个个体

    # 剪枝并重新计算适应度
    chrom_length = []
    for idx in range(len(chrom_list)):
        chrom_list[idx] = trim_chrom(chrom_list[idx])  # 修剪染色体
        fitness[idx] = evaluate(chrom_list[idx])  # 计算适应度
        chrom_length.append(chrom_list[idx][0].chrom.length())

    # 更新适应者生存选择
    chrom_list = select(chrom_list, fitness)

    # 确保最优个体在选择后不会被替换
    chrom_list[0] = best_chrom  # 始终保留最优个体
    fitness[0] = best_fitness

    # 更新 best.csv 文件，记录当前代的最优适应度
    with open(result_dir + "best.csv", "a", newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow([i+1, min(fitness)])

    # 保存每代信息到文件
    save_every_generation(i+1, chrom_list, fitness, chrom_length)

    # 保存当前代最优个体详细信息
    best_chrom_idx = np.argmin(fitness)
    best_chrom = chrom_list[best_chrom_idx][0]

    with open(result_dir + 'best_childs.txt', 'a') as f:
        f.write(f"Generation {i+1}\n")
        # 路径长度，网络中总反应个数，适应度值
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

# 算法结束后输出所有不重复路径到 final.txt
with open(result_dir + 'final.txt', 'w') as f:
    for path_str in unique_paths:
        f.write(path_str + '\n')


end_time = time.time()

log.logger.info('--Execution Time: %f h---' % ((end_time - start_time) / 3600))
