import os
import sys
import csv
import json
import numpy as np
import time
import argparse
import random
import shutil # 删除文件夹

## 以路径长度作为适应度值 烟花算法

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
    #     # "task_id": "009",
    #     # "algorithm": "single",
    #     # "pop_size": 10,  # 原本是10
    #     # "gen": 20,  # 原本是20
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
        # "task_id": "010",
        # "algorithm": "single",
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
    # args.save_path = ".."
    # args.project = '_'.join((args.task_id, args.ob_substrate, args.ob_product))
    # args.save_path = '/'.join((args.save_path, args.project)) + '/'
    # print("main_single:args.save_path")
    # print(args.save_path)
    # if not os.path.exists(args.save_path): os.mkdir(args.save_path)
    return args

args = parser()
cfg = get_config()
log = Logger(cfg['file_path']['log_dir'] + 'FWA2.log', fmt='%(message)s')

weight1 = 1 # 长度
weight2 = 0 # 吉布斯自由能
weight3 = 0 # 理论产量
weight4 = 0 # 设置权重 weight对应的是网络最少反应数

pop_size=20 # 每个子代的烟花数
# pop_size= 3 # 每个子代的烟花数
# max_iter = 50  # 观察结果30代之后就收敛了 基本上没变化了
max_iter = 30
# max_iter = 5
gaussianNum = 2 # 高斯变异烟花数
A = 20 # 调整爆炸半径的常数
# A = 10 # 调整爆炸半径的常数
M = 20 # 调整爆炸烟花个数的常数
# M = 3 # 调整爆炸烟花个数的常数
min_radius = 2

project = '_'.join(('FWA2_0415', args.ob_product))
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
    f.write("变异烟花数 gaussianNum = %d \n" % gaussianNum)
    f.write("调整爆炸半径的常数 A = %d \n" % A)
    f.write("爆炸最小半径 min_radius = %d \n" % min_radius)
    f.write("调整爆炸烟花个数的常数 M = %d \n" % M)

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
    # 每次评估 就会输出的个体
    with open(result_dir + 'all_fwa.txt', 'a') as f:
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
    return fit,length

# 爆炸算子
# 当某个烟花爆炸了100次，还没有爆炸成功， 则选取原烟花作为爆炸烟花。
def generate_explode_fireworks(chrom_list,explode_radius,explode_nums,log):
    # chrom_list, 上一代所有的个体
    # explode_radius, 上一代每个烟花对应的爆炸半径
    # explode_nums 上一代每个烟花对应的爆炸后生成的烟花的个数
    # explode_fireworks = []
    explode_fireworks = np.full((np.sum(explode_nums), 1), None, dtype=object)
    n = 0

    ## 对于每一个烟花
    for k in range(len(chrom_list)):
        radius = explode_radius[k]
        num = explode_nums[k]
        # 需要先深拷贝一下之前的节点 不然可能会出错
        Chrom = copy.deepcopy(chrom_list[k][0])
        log.logger.info("第 %d 个烟花发生爆炸，爆炸半径是 %d,爆炸火花数量为 %d，反应物列表长度为 %d" % (k+1, radius, num,len(Chrom.chrom.get_cpd_list())))
        log.logger.info("爆炸前的烟花是：")
        log.logger.info(Chrom.chrom.travel())
        # 遍历生成多个爆炸火花
        for j in range(num):
            log.logger.info("第 %d 个烟花发生爆炸，产生第 %d 个烟花" % (k+1, j+1))

            # 获取化合物 对于第n个反应，其起始化合物是n，终止化合物是n+1
            # 对于爆炸半径 radius 起始化合物范围 应该是 len(cpd_list)
            Chrom = copy.deepcopy(chrom_list[k][0]) # 这里一定要再声明一次 不然就会被上一次爆炸后的烟花覆盖掉
            cpd_list = Chrom.chrom.get_cpd_list()

            # 确定起始化合物 索引
            # 'a' cannot be empty unless no samples are taken 报错。因为len(cpd_list)-radius=0
            # 但是怎么会是0呢？化合物的长度，一定是大于反应的长度的。radius最大是反应长度。
            ## 设置一下 如果计算为0 则=1
            random_num = len(cpd_list)-radius
            if (random_num == 0):
                random_num = 1
                log.logger.info("len(cpd_list)-radius 计算为0")
            # cpd_start_idx = np.random.choice(np.arange(len(cpd_list)-radius),size=1)[0]
            cpd_start_idx = np.random.choice(np.arange(random_num),size=1)[0]
            # 确定终止化合物 索引 可能会出现爆炸半径比路径长的情况
            cpd_end_idx = cpd_start_idx + radius

            # 确定起始化合物和终止化合物
            cpd_start = cpd_list[cpd_start_idx]
            cpd_end = cpd_list[cpd_end_idx]

            log.logger.info("爆炸的起始化合物索引是：%d，起始化合物是 %s" % (cpd_start_idx, cpd_start))
            log.logger.info("爆炸的终止化合物索引是：%d, 终止化合物是 %s" % (cpd_end_idx, cpd_end))

            # 要拼接的节点
            # 原链表的第 cpd_start_idx 个节点的指针要指向 explode_middle_chrom 火花的第1个节点
            # explode_middle_chrom 火花的最后一个节点要指向 原链表的第 cpd_end_idx+1 个节点

            # 生成一个Chrom
            # explode_middle_chrom = Chromesome(args,cfg,log,initial_substrate=cpd_start,initial_product=cpd_end)
            explode_firework_n = 0
            while (1):
                explode_middle_chrom = Chromesome(args, cfg, log, initial_substrate=cpd_start,initial_product=cpd_end)
                # 会在 get_a_chrom 函数持续生成多次路径，直到生成合适的路径。
                # 但是不应该重复多次。所有这里传递进去一个数值，指定生成的次数。
                if (explode_firework_n == 5):
                    break
                explode_firework_n += 1
                if explode_middle_chrom.get_a_chrom(times=100):
                    # 如果爆炸了100次 还没有爆炸成功


                    # 已经找到了一条 中间两个化合物的路径
                    # 先剪枝
                    # explode_middle_chrom.trim_chrom(explode_middle_chrom.chrom)
                    # 都先指向头节点
                    cur = Chrom.chrom._head
                    explode_middle_chrom_cur = explode_middle_chrom.chrom._head
                    # 声明哑节点 哑节点指向cur
                    dummy = SingleReaction(S=None, P=None, reaction={})
                    dummy.next = cur
                    cur = dummy

                    # 遍历cur 指向第 cpd_start_idx 个节点。
                    # 但因为是从第一个节点开始往下遍历的，所以需要-1
                    for i in range(cpd_start_idx):
                        cur = cur.next

                    ## 必须先存储好Chrom的 第cpd_start_idx 个节点 与 第 cpd_end+1 个节点
                    cur_middle_start = cur

                    cur = Chrom.chrom._head

                    # 遍历cur 指向第 cpd_end+1 个节点
                    # 但因为是从第一个节点开始往下遍历的，所以需要-1
                    for i in range(cpd_end_idx):
                        # 因为有报错 所以这里添加了if 但是按理说 不应该有报错才对
                        # cur = cur.next
                        if cur!=None:
                            cur = cur.next
                    cur_middle_end = cur

                    # 将 cur 的next节点指向 explode_middle_chrom 的头节点 explode_middle_chrom_cur
                    cur_middle_start.next = explode_middle_chrom_cur # 报错 但是怎么会为null呢
                    # 遍历 explode_middle_chrom_cur 指向尾节点
                    while(explode_middle_chrom_cur.next != None):
                        explode_middle_chrom_cur = explode_middle_chrom_cur.next
                    # 将 explode_middle_chrom_cur 尾节点的 next指针指向 Chrom的 第 cpd_end+1 个节点
                    explode_middle_chrom_cur.next = cur_middle_end

                    ## 最后再指定一下头节点
                    Chrom.chrom._head = dummy.next
                    # # 拼接后 需要再做一次剪枝
                    # Chrom.trim_chrom(Chrom.chrom)
                    # 拼接后 需要再做一次剪枝 忘记赋值了
                    Chrom.chrom = Chrom.trim_chrom(Chrom.chrom)[0]

                    ## 应该再添加一个判断 如果爆炸后的路径长度大于20 重新爆炸
                    if(Chrom.chrom.length()>20):
                        ## !!! 长度大于20 后。chrom已经被改变了，应该恢复为原始的chrom。
                        Chrom = copy.deepcopy(chrom_list[k][0])
                        continue
                    else:
                        break

            # explode_fireworks.append(Chrom)
            explode_fireworks[n][0] = Chrom
            n += 1

            log.logger.info("爆炸后的烟花是：")
            log.logger.info(Chrom.chrom.travel())

    return explode_fireworks

# 产生变异的烟花
def generate_mutation_fireworks(chrom_list,gaussianNum):
    mutate_fireworks = np.full((gaussianNum, 1), None, dtype=object)
    n = 0

    # mutate_fireworks = []
    # 随机选择需要变异gaussianNum个烟花的序号
    mut_idx = np.random.choice(np.arange(len(chrom_list)), size=gaussianNum, replace=False)

    for idx in mut_idx:
        log.logger.info("第 %d 个烟花发生变异" % (idx))
        # 随机选择需要变异的两个烟花
        flag, mutation_firework = mut(chrom_list[idx])
        # mutate_fireworks.append(mutation_firework)
        mutate_fireworks[n][0] = mutation_firework
        n+=1

    return mutate_fireworks

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

def elite_select(all_pop,all_fitness,all_length):
    # 精英选择出的个体
    best_index = np.argmin(all_fitness)
    elite_individual = all_pop[best_index]
    elite_fitness = all_fitness[best_index]
    elite_length = all_length[best_index]

    # 计算适应度总和
    total_fitness = sum(all_fitness)
    # 计算每个个体被选择的概率，适应度值越小，概率越大
    probabilities = [1 / fit for fit in all_fitness]
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
    selected_length = [all_length[idx] for idx in selected_index]

    # 返回下一代个体、其对应的适应度值和个体长度
    pop = np.concatenate(([elite_individual], selected_pop))
    fitness = np.concatenate(([elite_fitness], selected_fitness))
    length = np.concatenate(([elite_length], selected_length))

    return pop, fitness, length

def save_every_generation(index,pop,fitness,chrom_length):
    with open(result_dir + 'every_generation_childs.txt', 'a') as f:
        f.write('==' * 5)
        # f.write(' ')
        f.write(' '+str(index)+' ')
        # f.write(' ')
        f.write('==' * 5)
        f.write('\n')
        for i in range(len(fitness)):

            # 路径长度，网络中总反应个数，适应度值（其实现在就是网络中总反应个数）
            length = chrom_length[i]
            # total_reactions = best_chrom.total_reactions
            fit = fitness[i]

            f.write(str([length, fit]))
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
chrom_length = []
for idx, chrom in enumerate(chrom_list):
    chrom = trim_chrom(chrom)
    # len是一个函数名 所以这里不能直接用len
    fit,len0 = evaluate(chrom)
    # 这个必须有 解决写入文件内best_enzyme_pool为空的问题
    chrom_list[idx] = chrom
    fitness.append(fit)
    chrom_length.append(len0)

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
save_every_generation(0,chrom_list,fitness,chrom_length)

# 迭代过程
for i in range(max_iter):
    log.logger.info("----- 第 %d 代烟花-----" %(i+1))
    explode_fireworks = [] # 爆炸烟花个体
    mutate_fireworks = [] # 变异烟花个体
    explode_fireworks_fitness = [] # 爆炸烟花适应度值
    mutate_fireworks_fitness = [] # 变异烟花适应度值

    explode_fireworks_length = []  # 爆炸烟花适应度值
    mutate_fireworks_length = []  # 变异烟花适应度值

    ## 上一代中最好的适应度值
    # 用一个路径长度数组存储路径长度
    best_chrom_idx = np.argmin(fitness)
    best_chrom = chrom_list[best_chrom_idx]

    # 当前种群中适应度值最小值
    fit_min = np.min(fitness)
    # 当前种群中适应度值最大值
    fit_max = np.max(fitness)

    ### 计算每个烟花的爆炸半径
    explode_radius = np.ceil(A * (fitness-fit_min+10e-6) / (np.sum(fitness-fit_min) + 10e-6))
    explode_radius = np.where(explode_radius <= 1, min_radius, explode_radius)  # 爆炸半径要>=2
    explode_radius = np.where(explode_radius >= chrom_length, chrom_length, explode_radius)  # 爆炸半径不能超过原路径长度 若超过则将其换位路径长度

    explode_radius = explode_radius.astype(int) # 转为整数
    ### 计算每个烟花可产生的爆炸烟花的个数
    explode_nums = np.ceil(M * (fit_max - fitness + 10e-6) / (np.sum(fit_max - fitness) + 10e-6))
    explode_nums = explode_nums.astype(int)  # 转为整数

    # 生成 爆炸烟花
    log.logger.info("----- explode -----")
    explode_fireworks = generate_explode_fireworks(chrom_list,explode_radius,explode_nums,log)

    # 生成变异烟花 gaussianNum 需要生成的变异烟花的个数
    log.logger.info("----- mutate -----")
    mutate_fireworks = generate_mutation_fireworks(chrom_list,gaussianNum)

    # 爆炸烟花剪枝并计算适应值
    log.logger.info("----- evaluate explode_fireworks -----")
    log.logger.info("共 %d 个爆炸火花" %(len(explode_fireworks)))
    for idx in range(len(explode_fireworks)):

        explode_fireworks[idx] = trim_chrom(explode_fireworks[idx])
        # explode_fireworks_fitness[idx] = evaluate(explode_fireworks[idx])
        fit,len0 = evaluate(explode_fireworks[idx])
        explode_fireworks_fitness.append(fit)
        explode_fireworks_length.append(len0)
    # 变异烟花剪枝并计算适应值
    log.logger.info("----- evaluate mutate_fireworks -----")
    for idx in range(len(mutate_fireworks)):
        # 报错  'numpy.ndarray' object has no attribute 'chrom'
        mutate_fireworks[idx] = trim_chrom(mutate_fireworks[idx])
        # mutate_fireworks_fitness.append(evaluate(mutate_fireworks[idx]))
        fit, len0 = evaluate(mutate_fireworks[idx])
        mutate_fireworks_fitness.append(fit)
        mutate_fireworks_length.append(len0)

    # 汇总所有的烟花 上一代烟花+爆炸火花+变异火花
    all_fireworks = np.concatenate([chrom_list, explode_fireworks, mutate_fireworks])
    # 汇总所有的适应度值
    all_fiteness = np.concatenate([fitness,explode_fireworks_fitness,mutate_fireworks_fitness])
    # ValueError: zero-dimensional arrays cannot be concatenated 报错
    all_length = np.concatenate([chrom_length,explode_fireworks_length,mutate_fireworks_length])

    # 选择下一代
    # 轮盘赌
    # chrom_list,fitness,chrom_length = select(all_fireworks,all_fiteness,all_length)
    # 精英保留策略
    chrom_list,fitness,chrom_length = elite_select(all_fireworks,all_fiteness,all_length)

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
    save_every_generation(i+1, chrom_list, fitness, chrom_length)

end_time = time.time()

log.logger.info('--Execution Time: %f h---' % ((end_time-start_time)/3600))