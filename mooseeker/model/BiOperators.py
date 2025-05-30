import time
import copy
import random
import numpy as np
from pymoo.core.sampling import Sampling
from pymoo.core.crossover import Crossover
from pymoo.core.mutation import Mutation
from pymoo.core.callback import Callback
from pymoo.core.duplicate import ElementwiseDuplicateElimination

from mooseeker.utils import get_config
from mooseeker.model import Chromesome
from mooseeker.utils import SingleLinkList, SingleReaction
from mooseeker.log import Logger


class BioSampling(Sampling):
    def _do(self, problem, n_samples, **kwargs):

        problem.log.logger.info("-----Sampling-----")
        X = np.full((n_samples, 1), None, dtype=object) # n_samples行 1列 初始值为None 类型为object
        count = 0
        start_time = time.time()
        while count < n_samples:
            chrom = Chromesome(problem.args, problem.cfg, problem.log)
            if chrom.get_a_chrom():
                X[count, 0] = chrom 
                problem.log.logger.info('This is the %dth chromsome'%(count))
                problem.log.logger.info(chrom.chrom.travel())
                count = count + 1
            with open(problem.args.save_path+"biosampling_chrom.txt", 'a+') as f:
                f.write(f"This is the {count}th chrom, it has tried {chrom.count} times. Chrom is {chrom.chrom.travel()}")
                f.write('\n')

        time_cost = time.time() - start_time
        with open(problem.args.save_path+"biosampling_chrom.txt", 'a+') as f:
                f.write(f"The computational cost of Sampling is {time_cost}s = {time_cost/60}min = {time_cost/3600}h.\n")

        problem.log.logger.info("-----Sampling Done-----")
        return X

# 交叉
class BioCrossover(Crossover):
    def __init__(self):
        # define the crossover: number of parents and number of offsprings
        super().__init__(2, 2)

    def _do(self, problem, X, **kwargs):
        problem.log.logger.info("-----Crossover-----")
        # The input of has the following shape (n_parents, n_matings, n_var)
        n_parents_, n_matings, n_var = X.shape

        # The output owith the shape (n_offsprings, n_matings, n_var)
        # Because there the number of parents and offsprings are equal it keeps the shape of X
        Y = np.full_like(X, None, dtype=object)

        for idx in range(n_matings):
            # whether there is the same point between parents
            if self.has_same_node(X[0, idx, 0], X[1, idx, 0]):
                Y[0, idx, 0], Y[1, idx, 0] = self.cross(X[0, idx, 0], X[1, idx, 0], problem)
            else:
                Y[0, idx, 0], Y[1, idx, 0] = X[0, idx, 0], X[1, idx, 0]
        problem.log.logger.info("-----Crossover Done-----")
        return Y

    def cross(self, chromA, chromB, problem):
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

                    problem.log.logger.info("-----Crossover Successfully!-----")
                    problem.log.logger.info("Chromesomes before crossover is :")
                    problem.log.logger.info(chromA.chrom.travel())
                    problem.log.logger.info('--'*5)
                    problem.log.logger.info(chromB.chrom.travel())
                    problem.log.logger.info("Chromesomes after crossover is :")
                    problem.log.logger.info(newA.chrom.travel())
                    problem.log.logger.info('--'*5)
                    problem.log.logger.info(newB.chrom.travel())
                    problem.log.logger.info('--'*10)

                    # 这里的newA因该是是Chromsome对象
                    # 对newA 和 newB进行剪枝
                    newA.chrom,flagA = newA.trim_chrom(newA.chrom)
                    newB.chrom, flagB = newA.trim_chrom(newB.chrom)
                    if(flagA):
                        problem.log.logger.info('=====> crossover after trim_chrom A:  %s<=====' % (newA.chrom.travel()))
                    if(flagB):
                        problem.log.logger.info('=====> crossover after trim_chrom B:  %s<=====' % (newB.chrom.travel()))
                    return newA, newB
                else:
                    curB = curB.next
            curA = curA.next

    def has_same_node(self, chromA, chromB):
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

# 变异
class BioMutation(Mutation):
    def __init__(self):
        super().__init__()

    def _do(self, problem, X, **kwargs):
        problem.log.logger.info("-----Mutation-----")
        try:
            Mut_X = self.Muta(X, problem)
            return Mut_X
        except:
            mut_x = Chromesome(problem.args, problem.cfg, problem.log)
            mut_x.get_a_chrom()
            idx = np.random.randint(low=0, high=len(X))
            X[idx] = mut_x
            return X

    def Muta(self, chromesomes, problem):
        count = 0
        while (count < 100):
            idx = np.random.randint(len(chromesomes))
            cond, mut_cross = self.mut(chromesomes[idx], problem)
            if cond:
                chromesomes[idx] = mut_cross
                return chromesomes

    def mut(self, Chrome, problem):
        chrom = copy.deepcopy(Chrome)
        # 确定了底物和产物
        m, n = random.sample(range(0, chrom[0].chrom.length()), 2)
        cur = chrom[0].chrom._head
        for i in range(0, chrom[0].chrom.length()):
            if i == min(m, n):
                s = cur.P
                cur = cur.next
            elif i == max(m, n):
                p = cur.P
                break
            else:
                cur = cur.next

        # generate mut_chrom
        mut_chrom = np.full((1, 1), None, dtype=object)
        while (1):
            temp_chrom = Chromesome(problem.args, problem.cfg, problem.log)
            if temp_chrom.get_a_chrom():
                mut_chrom[0] = temp_chrom
                break

        # connect chrom and mut_chrom
        cur = chrom[0].chrom._head
        mut_cur = mut_chrom[0][0].chrom._head

        for i in range(1, max(m, n) + 1):
            if i == min(m, n):
                temp = cur.next
                cur.next = mut_cur
                while (mut_cur != None):
                    if mut_cur.next == None:
                        mut_cur.next = temp
                        # trim new chrom
                        chrom[0].trim_chrom()
                        # update chrom's pool
                        chrom[0].update_pool_after_change()
                        if chrom[0].chrom != Chrome[0].chrom:
                            
                            problem.log.logger.info("-----Mutation Successfully!-----")
                            problem.log.logger.info("Chromesomes before mutation is :")
                            problem.log.logger.info(Chrome.chrom.travel())
                            problem.log.logger.info("Chromesomes after mutation is :")
                            problem.log.logger.info(chrom.chrom.travel())
                            problem.log.logger.info('--'*10)

                            return [True, chrom]
                        else:
                            # problem.log.logger.info(chrom.travel())
                            problem.log.logger.info("-----Mutation Failed! Initialize a new chromesome!-----")
                            return [False, chrom]
                    else:
                        mut_cur = mut_cur.next
            else:
                cur = cur.next 

class BioDuplicateElimination(ElementwiseDuplicateElimination):
    def is_equal(self, a, b):
        return set(a.X[0].chrom.get_rxn_list()) == set(b.X[0].chrom.get_rxn_list())

class BioCallback(Callback):
    def __init__(self) -> None:
        super().__init__()
        self.data["best"] = []

    def notify(self, algorithm):
        self.data["best"].append(algorithm.pop.get("F").min())




def cross(chromA, chromB):
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
                print("-----Crossover Successfully!-----")
                return newA, newB
            else:
                curB = curB.next
        curA = curA.next
    print("-----Crossover Failed!-----")

def has_same_node(chromA, chromB):
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

def biomut(Chrome, args, cfg, log):
    chrom = copy.deepcopy(Chrome)
    # 确定了底物和产物
    m, n = random.sample(range(0, chrom[0].chrom.length()), 2)
    cur = chrom[0].chrom._head
    for i in range(0, chrom[0].chrom.length()):
        if i == min(m, n):
            s = cur.P
            cur = cur.next
        elif i == max(m, n):
            p = cur.P
            break
        else:
            cur = cur.next

    # generate mut_chrom
    mut_chrom = np.full((1, 1), None, dtype=object)
    while (1):
        temp_chrom = Chromesome(args, cfg, log)
        if temp_chrom.get_a_chrom():
            mut_chrom[0] = temp_chrom
            break

    # connect chrom and mut_chrom
    cur = chrom[0].chrom._head
    mut_cur = mut_chrom[0][0].chrom._head

    for i in range(1, max(m, n) + 1):
        if i == min(m, n):
            temp = cur.next
            cur.next = mut_cur
            while (mut_cur != None):
                if mut_cur.next == None:
                    mut_cur.next = temp
                    # trim new chrom
                    chrom[0].trim_chrom()
                    # update chrom's pool
                    chrom[0].update_pool_after_change()
                    if chrom[0].chrom != Chrome[0].chrom:
                        print("-----Mutation Sucessfully!-----")
                        return [True, chrom]
                    else:
                        # problem.log.logger.info(chrom.travel())
                        print("-----Mutation Failed! Initialize a new chromesome!-----")
                        return [False, chrom]
                else:
                    mut_cur = mut_cur.next
        else:
            cur = cur.next 




