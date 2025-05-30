import os
import sys
import csv
import json
import numpy as np
import time
import argparse
import random
import shutil # 删除文件夹

## 以吉布斯自由能作为适应度值 使用烟花算法

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
from mooseeker.fit_funcs import get_gibbs_mlp, get_yield,get_network_totalReactions

def parser(target_compound):
    if(target_compound=="C00022"):
        params = {
            # "task_id": "009",
            # "algorithm": "single",
            # "pop_size": 10,  # 原本是10
            # "gen": 20,  # 原本是20
            "NrMax": 15,
            "ob_substrate": "C00267",
            "ob_product": "C00022",
            "abundant": ["C00001", "C00002", "C00080"],
            "database": ["KEGG", "BiGG"],
            "host": ["ecoli"],
            "eva_func": ["yield"]
            # "eva_func": ["length", "gibbs", "yield"]
        }
    elif(target_compound=="C00755"):
        params = {
            # "task_id": "010",
            # "algorithm": "single",
            # "pop_size": 10,  # 原本是10
            # "gen": 20,  # 原本是20
            "NrMax": 15,
            # "ob_substrate": "C00082", # 酪氨酸
            "ob_substrate": "C01494", # 阿魏酸
            "ob_product": "C00755",
            "abundant": ["C00001", "C00002", "C00080"],
            "database": ["KEGG", "BiGG"],
            "host": ["ecoli"],
            # "eva_func": ["yield"]
            # "eva_func": ["length", "gibbs", "yield"]
        }

    args = argparse.Namespace(**params)
    # args.save_path = ".."
    # args.project = '_'.join((args.task_id, args.ob_substrate, args.ob_product))
    # args.save_path = '/'.join((args.save_path, args.project)) + '/'
    # print("main_single:args.save_path")
    # print(args.save_path)
    # if not os.path.exists(args.save_path): os.mkdir(args.save_path)
    return args

weight1 = 0 # 长度
weight2 = 1 # 吉布斯自由能
weight3 = 0 # 理论产量
weight4 = 0 # 设置权重 weight对应的是网络最少反应数

pop_size=20 # 每个子代的烟花数
# pop_size= 3 # 每个子代的烟花数
max_iter = 50
# max_iter = 5
# gaussianNum = 2 # 高斯变异烟花数
# A = 20 # 调整爆炸半径的常数
# A = 10 # 调整爆炸半径的常数
# M = 20 # 调整爆炸烟花个数的常数
# M = 3 # 调整爆炸烟花个数的常数
# min_radius = 1
CROSS_RATE = 0.8
MUTATION_RATE = 0.4

target_compound = "C00755"
# target_compound = "C00022"

args = parser(target_compound)
cfg = get_config()
log = Logger(cfg['file_path']['log_dir'] + 'FWA2.log', fmt='%(message)s')

project = '_'.join(('ga_gibbs_0423_4', args.ob_product))
# project = '_'.join(('My', args.production, str(int(args.weight1 * 100)), str(int(args.weight2 * 100)),
#                     str(int(args.weight3 * 100))))
result_dir = cfg['file_path']['result_dir'] + project + '/'

if os.path.exists(result_dir):
    # 先删除之前的全部文件
    shutil.rmtree(result_dir)

if not os.path.exists(result_dir):
    # 再创建空的文件夹
    os.mkdir(result_dir)
## 写入运行参数
with open(result_dir + 'parameters.txt', 'w') as f:
    f.write("权重分配 \n")
    f.write("[路径长度，吉布斯自由能，理论产量，网络总反应数] \n")
    f.write(str([weight1,weight2,weight3,weight4]))
    f.write('\n')
    f.write('\n')

    f.write("pop_size = %d \n" % pop_size)
    f.write("max_iter = %d \n" % max_iter)
    f.write("交叉率 CROSS_RATE = %f \n" % CROSS_RATE)
    f.write("变异率 MUTATION_RATE = %f \n" % MUTATION_RATE)
    f.write("ob_substrate ----> %s \n" % args.ob_substrate)
    f.write("ob_product ----> %s \n" % args.ob_product)

def evaluate(X):
    weight1 = 0  # 长度
    weight2 = 1  # 吉布斯自由能
    weight3 = 0  # 理论产量
    weight4 = 0  # 设置权重 weight对应的是网络最少反应数

    #log.logger.info('---get pathway length---')
    length = X[0].chrom.length()

    #log.logger.info('---get pathway gibbs---')
    gibbs = get_gibbs_mlp(cfg,X[0].chrom.get_rxn_list())

    #log.logger.info('---get pathway yield---')
    # _yield = get_yield(cfg,X[0].chrom.get_rxn_list())

    # log.logger.info('---get nework total reactions---')
    # total_reactions = get_network_totalReactions(cfg,X[0],log)

    # fit = float(weight1) * length + float(weight2) * gibbs - float(weight3) * _yield
    fit = float(weight2) * gibbs
    # fit = float(weight4) * total_reactions
    # fit = float(weight1) * length

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
    # 每次评估 就会输出的个体
    with open(result_dir + 'all_ga.txt', 'a') as f:
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

## 对某个个体实行变异 # 也是随机选择的变异点 只不过终点都是目标化合物
def mut(Chrome):
    #
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

        return [True, chrom[0]]
    else:
        # problem.log.logger.info(chrom.travel())
        log.logger.info("-----Mutation Failed!-----")
        return [False, chrom[0]]

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

## 根据适应度值 执行轮盘赌选择策略
## 后续可以调整为 精英选择策略
def select(all_pop,all_fitness,all_length):
    # 计算适应度总和
    total_fitness = sum(all_fitness)
    # 计算每个个体被选择的概率，适应度值越小，概率越大
    probabilities = [1 / fit for fit in all_fitness]
    total_probability = sum(probabilities)
    # 归一化概率，使概率之和为1
    normalized_probabilities = [prob / total_probability for prob in probabilities]

    # 使用numpy进行随机选择 size 应该就是pop_size
    selected_index = np.random.choice(np.arange(len(all_fitness)), size=pop_size, replace=True, p=normalized_probabilities)
    pop = [all_pop[idx] for idx in selected_index]
    fitness = [all_fitness[idx] for idx in selected_index]
    length = [all_length[idx] for idx in selected_index]
    return pop,fitness,length

# 因为吉布斯自由能值不能为负。所以说如果是正的，那么就将其概率设为0.
def elite_select(all_pop,all_fitness):
    # 精英选择出的个体
    best_index = np.argmin(all_fitness)
    elite_individual = all_pop[best_index]
    elite_fitness = all_fitness[best_index]


    # 计算每个个体被选择的概率，适应度值越小，概率越大
    # 计算每个个体被选择的概率
    probabilities = []
    # 为正时 舍弃掉这个个体
    for fit in all_fitness:
        if fit >= 0:
            probabilities.append(0)
        else:
            # 负值 取绝对值。也就是绝对值越大，概率越大
            probabilities.append(abs(fit))

    total_probability = sum(probabilities)
    # 归一化概率，使概率之和为1
    normalized_probabilities = [prob / total_probability for prob in probabilities]

    # 使用numpy进行随机选择 size 应该就是pop_size-1
    # 剩下的N-1个烟花的选择使用轮盘赌的方法在候选者集合中进行选择。不需要剔除精英个体
    selected_index = np.random.choice(np.arange(len(all_fitness)), size=pop_size-1, replace=True,
                                      p=normalized_probabilities)

    # 轮盘赌选择的个体
    selected_pop = [all_pop[idx] for idx in selected_index]
    selected_fitness = [all_fitness[idx] for idx in selected_index]

    # 返回下一代个体、其对应的适应度值和个体长度
    pop = np.concatenate(([elite_individual], selected_pop))
    fitness = np.concatenate(([elite_fitness], selected_fitness))

    return pop, fitness

def save_every_generation(index,pop,fitness):
    with open(result_dir + 'every_generation_childs.txt', 'a') as f:
        f.write('==' * 5)
        # f.write(' ')
        f.write(' '+str(index)+' ')
        # f.write(' ')
        f.write('==' * 5)
        f.write('\n')
        for i in range(len(fitness)):

            # 路径长度，网络中总反应个数，适应度值（其实现在就是网络中总反应个数）
            # total_reactions = best_chrom.total_reactions
            fit = fitness[i]

            f.write(str([fit]))
            f.write('\n')
            # 写入路径中的反应及化合物
            chrom = pop[i][0]
            travel_list = chrom.chrom.travel()
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
            f.write('\n')

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
    # len是一个函数名 所以这里不能直接用len
    fit = evaluate(chrom)
    # 这个必须有 解决写入文件内best_enzyme_pool为空的问题
    chrom_list[idx] = chrom
    fitness.append(fit)

chrom_next = chrom_list

# 存储每一代的爆炸烟花
explode_fireworks = []
# 存储每一代的变异烟花
mutate_fireworks = []


# 写入初代烟花 最好的个体
with open(result_dir+"best.csv", "a", newline='') as csvfile:
    csvwriter = csv.writer(csvfile)
    # fiteness
    csvwriter.writerow([0, min(fitness)])

# 记录当前代 最好的子代 子代travel，长度，总反应个数，适应度值，酶组合
best_chrom_idx = np.argmin(fitness)
best_chrom = chrom_list[best_chrom_idx][0]

with open(result_dir + 'best_childs.txt', 'a') as f:
    f.write(str(0))
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

# 写入初代烟花的 所有个体
save_every_generation(0,chrom_list,fitness)

# 迭代过程
for i in range(max_iter):
    ## 上一代中最好的适应度值
    # 用一个路径长度数组存储路径长度
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
        chrom_list[index1] = Y[0]  # 替换为交叉后的个体
        chrom_list[index2] = Y[1]  # 替换为交叉后的个体

        log.logger.info("-----Crossover Done-----")

    # idx = np.random.randint(len(chrom_list))
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
        log.logger.info('=====> Mututation of the %dth iter %dth particle<=====' % (i + 1, idx + 1))

    # 剪枝，计算适应值
    for idx in range(len(chrom_list)):
        chrom_list[idx] = trim_chrom(chrom_list[idx])
        fitness[idx] = evaluate(chrom_list[idx])
    # print("fitness:")
    # print(fitness)
    # print(np.shape(fitness))

    # 适者生存 适应度值越小越好
    chrom_list,fitness = elite_select(chrom_list, fitness)


    #
    with open(result_dir+"best.csv", "a", newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        # fiteness
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

    # 记录当前代select后的个体
    save_every_generation(i+1, chrom_list, fitness)

end_time = time.time()

log.logger.info('--Execution Time: %f h---' % ((end_time-start_time)/3600))