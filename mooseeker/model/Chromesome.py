# -*- coding: utf-8 -*-
# @Author: Caoyh
# @Date:   2022-05-28 20:45:14
# @Last Modified time: 2022-06-15 12:02:42

import sys
import os

script_dir = os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
sys.path.append(script_dir)

import numpy as np
import copy

from mooseeker.utils import get_config
from mooseeker.utils import SingleLinkList
from mooseeker.utils import SingleReaction
from mooseeker.log import Logger


class Chromesome(object):
    def __init__(self, args, cfg, log, initial_substrate="",initial_product = "",chrom=None):

        self.args = args
        self.cfg = cfg
        self.log = log

        # 目标底物。这里需要考虑当变异的时候，是从中间随便一个点开始变异的，也就是起始化合物需要是中间的一个点。
        # 而不是最出定好的底物。
        if (initial_substrate == ""):
            self.ob_substrate = args.ob_substrate
        else:
            self.ob_substrate = initial_substrate

        # 在生成爆炸烟花的时候 需要生成两个中间化合物之间的路径
        if (initial_product == ""):
            self.ob_product = args.ob_product
        else:
            self.ob_product = initial_product

        # self.ob_substrate = args.ob_substrate # 目标底物
        # self.ob_product = args.ob_product # 目标产物
        self.abundant = args.abundant # 底物池
        self.NrMax = args.NrMax # 路径最大长度

        # 加载池子 MYPOOL_20221120
        self.all_pool = np.load(cfg['file_path']['mypool'], allow_pickle=True).item()
        
        self.ava_pool = dict() # 空字典
        if chrom == None:
            self.chrom = SingleLinkList()
        else:
            self.chrom = chrom

        self.count = 0
        # 存储当前chrom 最好的酶组合情况
        self.best_enzyme_pool = []
        self.total_reactions = 0 #网络中的总反应个数
        self.best_total_rxn_list = []
        self.fit = 0 # 记录当前chrom适应度值
        self.total_compounds = []
        self.best_total_rxn_dict = {} 

    def set_best_network(self,best_enzyme_pool,best_total_rxn_list,total_reactions,total_compounds):
        self.best_enzyme_pool = best_enzyme_pool
        self.total_reactions = total_reactions
        self.best_total_rxn_list = best_total_rxn_list
        self.total_compounds = total_compounds

    def set_fit(self,fit):
        self.fit = fit

    ## 循环了10000次还没找到的话，就是没找到
    def get_a_chrom(self,times=10000):
        count = 0
        while 1:
            # 这里改小一点 防止一直递归
            if count == times:
                self.count = count
                return False
            elif self.get_a_chrom_roulette():
                # with open(self.ages.save_path+"chrom_count.txt", 'a') as f:
                #     f.write(count)
                #     f.write('\n')
                self.count = count
                return True
            else:
                count += 1

    ## 加载底物池
    def init_ava_pool(self):
        """
        重新初始化 
        :return:
        """
        all_pool_cpds = list(self.all_pool.keys()) # 长度6790 这个键是化合物，对应的值是一个SingleReaction对象 所有的反应对
        for c in self.abundant + [self.ob_substrate]: # 遍历底物池中的所有 化合物
            if c not in all_pool_cpds: # 如果c不是有效的化合物（这里其实是排除了辅因子和error）
                continue
            else:
                temp = list(self.all_pool[c]) # 加载反应物c 对应的所有反应 多个SingleReaction 对象
                self.ava_pool[c] = temp # ava_pool 是一个字典。键，化合物c。值，化合物c对应的所有反应。
                # 因为是遍历 所以ava_pool中就存储了，当前底物池中，所有化合物对应的所有反应。

    def update_ava_pool(self, Node):
        """
        根据给定的化合物 更新 pool
        :param compound:
        :return:
        """
        # substrates, products = split_equation(Node.reaction['equation'])
        # Node 是我们选出来的一个反应 SingleReaction
        reactants, products = Node.reaction['reactants'], Node.reaction['products']
        all_pool_cpds = list(self.all_pool.keys()) # 长度6790 这个键是化合物，对应的值是一个SingleReaction对象
        for c in reactants + products:
            if c not in all_pool_cpds:
                continue
            temp = list(self.all_pool[c]) # 存储的是 选定反应物所对应的 所有 反应物——底物池子
            if c in self.ava_pool: # 如果 化合物c 原本就在底物池中 则可以跳过
                continue
            else: # 如果化合物c原本不在底物池中
                self.ava_pool[c] = temp # 这里是说 我选出来了一条路径后 我就视为，参与反应的所有反应物，不管之前有没有在底物池中，我都视为它在底物池中了。

    def update_pool_after_change(self):
        """ 
        After trim, crossover and mutation
        Update pool of chrom
        return self.ava_pool, self.compounds
        """
        self.init_ava_pool()
        cur = self.chrom._head 
        while (cur != None):
            self.update_ava_pool(cur)
            cur = cur.next
        # self.get_all_compounds()

    def roulette(self, NodeList):
        """轮盘赌算法
        :param NodeList: NodeList 是一个链表
        :return: 返回索引值
        """
        sum_val = 0
        for i in range(len(NodeList)):
            ## Tanimoto 有个工具可以直接计算相似度 直接根据MACCS分子指纹来计算的 如果相似性高，则被选择的概率就大。
            # 相似度比较的是，底物s1，在底物-产物对池子中，目标产物与生成物之间的相似性
            if NodeList[i].Tanimoto != 0 and NodeList[i].Tanimoto != None:
                sum_val = sum_val + NodeList[i].Tanimoto

        random_val = np.random.random()
        probability = 0  #累计概率
        for i in range(len(NodeList)):
            if NodeList[i].Tanimoto != 0 and NodeList[i].Tanimoto != None:
                probability = probability + NodeList[i].Tanimoto / sum_val
            if probability >= random_val:
                return i
            else:
                continue

    def get_a_chrom_roulette(self):
        # 按照底物相似性 根据轮盘赌算法（多样性） 生成个体

        self.init_ava_pool() # 初始化

        # 创建一个单链表格式的变量 存储 第一个是值，第二个存储的是指向下一个节点的地址
        mychrom = SingleLinkList()
        # roulette函数 轮盘赌法 选择出来一个 index  也就是一个反应。
        index = self.roulette(list(self.ava_pool[self.ob_substrate])) # ob_substrate 目标底物
        # 可以看出来 Node就是一个 SingleReaction 对象
        Node = copy.deepcopy(list(self.ava_pool[self.ob_substrate])[index]) # ava_pool底物池 键是所有的化合物，值是化合物对应的所有反应
        self.update_ava_pool(Node) # 使用选出来的这条反应 更新池子。
        mychrom.append(Node) # 添加节点

        ## 设置的最大路径是100 路径长度小于100时
        while (Node.P != self.ob_product and mychrom.length() < 100):
            index = self.roulette(list(self.ava_pool[Node.P])) #  将 产物 所涉及的所有对子，进行轮盘赌，选出一条反应index
            Node = copy.deepcopy(list(self.ava_pool[Node.P])[index]) # 选出一条反应 SingleReaction
            self.update_ava_pool(Node) # 用这个反应更新池子
            mychrom.append(Node) # 链表添加节点

        # 找到目标产物 判断中间产物 是不是目标化合物
        if Node.P == self.ob_product:
            # Get a chromesome sucessfully and trim the chrom
            ## 剪枝 就是师姐说的 会有冗余
            chrom,flag = self.trim_chrom(mychrom)
            # 没有重复化合物的路径 是否符合最初设定的长度 NrMax设定为20  如果剪枝后>20 那么就不要它
            # check_chrom_available 有重复反应
            # 如果经过剪枝了,输出剪枝后的路径
            if flag:
                self.log.logger.info('=====> after trim_chrom:  长度是%s---%s<=====' %(chrom.length(),chrom.travel()))
            if (chrom.length() < self.NrMax) and self.check_chrom_available(chrom):
                self.chrom = chrom # 得到一条路径 也就是一个 SingleLinkList 对象。
                return True
        # 清空池子
        self.ava_pool.clear()
        return False

    # 如果剪枝后有相同的反应 那就说明发生了可逆反应 应该是单向的才对 所以要删除这种情况下产生的子代。
    def check_chrom_available(self, chrom):
        """检查是否存在同一个reaction出现两次的情况
        """
        head = chrom._head
        left = head
        if (chrom != None and chrom.length() > 2):
            while(left.next != None):
                if left.reaction['kegg_id'] == left.next.reaction['kegg_id']:
                    return False
                else:
                    left = left.next
        return True

    # 做了一些修改
    def trim_chrom(self, input_chrom):
        chrom = copy.deepcopy(input_chrom)

        head = chrom._head
        # 修剪情况: a->b b->c c->a a->d ===> a->d 
        if (chrom != None and chrom.length() > 2):
            dummy = SingleReaction(S=None, P=None, reaction={})
            dummy.next = head
            pre = dummy
            left = pre.next
            while(left != None):
                right = left.next
                while(right != None):
                    if left.S == right.S: 
                        pre.next = right
                        # 自己修改的部分 原本的无法修改开头冗余的路径
                        if (pre.S == None):
                            chrom._head = pre.next
                        # pre = dummy
                        left = pre.next
                        # break  # 去掉break 因为有可能重复两次 a->b->a->b->a->c->d
                    right = right.next
                pre = pre.next
                left = left.next
                pre.next = left  # 添加的一行
            # true表示被剪枝。false表示没有被剪枝
            return chrom,True
        else:
            return input_chrom,False

if __name__ == '__main__':

    params = {
        'task_id': '001',
        'pop_size': 20,  # 种群大小
        'gen': 30,  # 迭代次数
        'NrMax': 10,  # 路径最大长度
        'ob_substrate': 'C00082',
        'ob_product': 'C00755',  # 目标化合物
        'abundant': ['C00001', 'C00002', 'C00080'],  # 可选初始化合物
        'database': ['KEGG', 'BiGG'],  # 数据库
        'host': 'ecoli',  # 宿主细胞
        'eva_func': ['length', 'gibbs', 'yield', ],  # 评估函数
    }

    from mooseeker.utils import parser_args
    args = parser_args(params)

    cfg = get_config()
    log = Logger('../cache/log/chromesome.log')

    cla = Chromesome(args, cfg, log)
    cla.get_a_chrom()

    rxn_list = cla.chrom.travel()
    print(rxn_list)
