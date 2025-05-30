import os
import sys
import csv
import json
import numpy as np
import time
import argparse

## 以路径长度作为适应度值 越小越好


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

from mooseeker.model.BiOperators import *
from mooseeker.utils import *
from mooseeker.fit_funcs import get_gibbs, get_yield,get_network_totalReactions

def parser():
    # params = {
    #     "task_id": "009",
    #     "algorithm": "single",
    #     "pop_size": 10,  # 原本是10
    #     "gen": 20,  # 原本是20
    #     "NrMax": 15,
    #     "ob_substrate": "C00267",
    #     "ob_product": "C00022",
    #     "abundant": ["C00001", "C00002", "C00080"],
    #     "database": ["KEGG", "BiGG"],
    #     "host": ["ecoli"],
    #     "eva_func": ["yield"]
    #     # "eva_func": ["length", "gibbs", "yield"]
    # }
    params = {
        "task_id": "014",
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
log = Logger(cfg['file_path']['log_dir'] + 'GA2_0407.log', fmt='%(message)s')

weight1 = 1 # 长度
weight2 = 0 # 吉布斯自由能
weight3 = 0 # 理论产量
weight4 = 0 # 设置权重 weight对应的是网络最少反应数

# pop_size=20
pop_size=20
# max_iter = 50
max_iter = 50
CROSS_RATE = 0.8
MUTATION_RATE = 0.4 # 我觉得根据这个路径的特殊性 变异的概率可以大一点


project = '_'.join(('GA2_0407', args.ob_product))
# project = '_'.join(('My', args.production, str(int(args.weight1 * 100)), str(int(args.weight2 * 100)),
#                     str(int(args.weight3 * 100))))
result_dir = cfg['file_path']['result_dir'] + project + '/'
if not os.path.exists(result_dir): os.mkdir(result_dir)

if os.path.exists(result_dir + 'all_ga.txt'):
    os.remove(result_dir + 'all_ga.txt')

def evaluate(X):
    weight1 = 1
    weight2 = 0
    weight3 = 0
    weight4 = 0

    #log.logger.info('---get pathway length---')
    length = X[0].chrom.length()

    #log.logger.info('---get pathway gibbs---')
    # gibbs = get_gibbs(X[0].chrom.get_rxn_list())

    #log.logger.info('---get pathway yield---')
    # _yield = get_yield(X[0].chrom.get_rxn_list())

    # log.logger.info('---get nework total reactions---')
    # total_reactions = get_network_totalReactions(cfg,X[0],log)

    # fit = float(weight1) * length + float(weight2) * gibbs - float(weight3) * _yield
    # fit = float(weight4) * total_reactions
    fit = float(weight1) * length


    #print('The chromesome is:\n', X.chrom.travel())
    # print('Length:%f \nGibbs:%f \nYield:%f \nWeighted Score:%f \n' % (length, gibbs, _yield, fit))
    log.logger.info('Length:%f \nWeighted Score:%f \n' % (length, fit))

    # 每次都需要删除一下
    csv_file_path = result_dir + 'evaluation_data.csv'
    if not os.path.exists(csv_file_path):
        with open(csv_file_path, 'w', newline='') as file:
            writer = csv.writer(file)
            # writer.writerow(["fit", "gibbs", "yield", "length"])
            writer.writerow(["fit", "length"])

    with open(csv_file_path, 'a', newline='') as file:
        writer = csv.writer(file)
        # writer.writerow([fit, gibbs, _yield, length])
        writer.writerow([fit, length])

    # 在评估的地方保存个体和适应度值 注意这个文件每次运行都需要删除一下 不然每次都是追加
    with open(result_dir + 'all_ga_0407.txt', 'a') as f:
        # f.write(str([length, gibbs, _yield, fit]))
        f.write(str([length, fit]))
        f.write('\n')
        f.write(str(X[0].chrom.get_rxn_list()))
        f.write('\n')
        # log.logger.info('---The chromesome is---')
        travel_list = X[0].chrom.travel()
        for line in travel_list:
            for l in line:
                f.write(l)
                f.write('\t')
            f.write('\n')
        f.write('==' * 10)
        f.write('\n')

    X[0].set_fit(fit)
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
def select(pop,fitness):
    # 计算适应度总和
    total_fitness = sum(fitness)
    # 计算每个个体被选择的概率，适应度值越小，概率越大
    probabilities = [1 / fit for fit in fitness]
    total_probability = sum(probabilities)
    # 归一化概率，使概率之和为1
    normalized_probabilities = [prob / total_probability for prob in probabilities]

    # 使用numpy进行随机选择
    selected_index = np.random.choice(np.arange(len(fitness)), size=pop_size, replace=True, p=normalized_probabilities)
    pop = [pop[idx] for idx in selected_index]
    return pop

count = 0
#chrom_list = []
start_time = time.time()  # 开始计时

log.logger.info("-----Sampling-----")
# 生成初始的个体
chrom_list = np.full((pop_size, 1), None, dtype=object)
while (count < pop_size):
    chrom = Chromesome(args, cfg, log)
    if chrom.get_a_chrom():
        log.logger.info('--This is the %dth chromsome---' % (count+1))
        #chrom_list.append(chrom)
        chrom_list[count,0] = chrom
        count += 1
log.logger.info("-----Sampling Done-----")

with open(result_dir+"best.csv", "w", newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    csvwriter.writerow(["Generation", "Best Fitness"])

# 计算初始个体的适应度值
fitness=[]
for idx, chrom in enumerate(chrom_list):
    chrom = trim_chrom(chrom)
    fit = evaluate(chrom)
    fitness.append(fit)

chrom_next = chrom_list

# 迭代过程
for i in range(max_iter):
    best_chrom_idx = np.argmin(fitness)
    best_chrom = chrom_list[best_chrom_idx]

    log.logger.info("-----Crossover-----")
    # 从chrom_list中选择两个不同的索引值(除最优个体）
    remaining_indices = [idx for idx in range(len(chrom_list)) if idx != best_chrom_idx]
    indices = np.random.choice(remaining_indices, 2, replace=False)

    index1 = indices[0]
    index2 = indices[1]

    # 交叉
    if np.random.rand() < CROSS_RATE:
        log.logger.info("----- start Crossover -----")
        # 如果概率在交叉的概率内才交叉
        Y = perform_crossover(chrom_list[index1], chrom_list[index2])
        chrom_list[index1] = Y[0] # 替换为交叉后的个体
        chrom_list[index2] = Y[1] # 替换为交叉后的个体

        log.logger.info("-----Crossover Done-----")

    #idx = np.random.randint(len(chrom_list))
    mutate_idx = np.random.choice(remaining_indices)

    # 变异
    if np.random.rand() < MUTATION_RATE:
        cond, mut_cross = mut(chrom_list[mutate_idx])
        if cond:
            chrom_list[mutate_idx] = mut_cross
        else:
            newchrom = Chromesome(args, cfg, log)
            if newchrom.get_a_chrom():
                chrom_list[mutate_idx][0] = newchrom
        log.logger.info('=====> Mututation of the %dth iter %dth particle<====='%(i+1,idx+1))

    # 剪枝，计算适应值
    for idx in range(len(chrom_list)):
        chrom_list[idx] = trim_chrom(chrom_list[idx])
        fitness[idx] = evaluate(chrom_list[idx])

    # 适者生存 适应度值越小越好
    chrom_list = select(chrom_list,fitness)

    #
    with open(result_dir+"best.csv", "a", newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow([i+1, min(fitness)])

    # 记录当前代 最好的子代 子代travel，长度，总反应个数，适应度值，酶组合
    best_chrom_idx = np.argmin(fitness)
    best_chrom = chrom_list[best_chrom_idx][0]

    with open(result_dir + 'best_childs.txt', 'a') as f:
        f.write(str(i+1))
        f.write('\n')
        # 路径长度，网络中总反应个数，适应度值（其实现在就是网络中总反应个数）
        length = best_chrom.chrom.length()
        # total_reactions = best_chrom.total_reactions
        fit = best_chrom.fit

        f.write(str([length, fit]))
        f.write('\n')
        # 写入路径中的反应及化合物
        travel_list = best_chrom.chrom.travel()
        for line in travel_list:
            for l in line:
                f.write(l)
                f.write('\t')
            f.write('\n')
        # # 写入最好的酶组合
        # f.write(str(best_chrom.best_enzyme_pool))
        # f.write('\n')
        # # 写入当前酶组合产生的总反应
        # f.write(str(best_chrom.best_total_rxn_list))
        # f.write('\n')
        f.write('==' * 10)
        f.write('\n')

end_time = time.time()

log.logger.info('--Execution Time: %f h---' % ((end_time-start_time)/3600))