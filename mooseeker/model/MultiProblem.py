"""
Autor: caoyh
Date: 2022-11-17 14:43:06
LastEditTime: 2022-11-20 16:44:54
"""
import numpy as np
from pymoo.core.problem import ElementwiseProblem
from mooseeker.fit_funcs import get_gibbs
# 扩展成网络后的理论产量计算
from mooseeker.fit_funcs import get_yield


class MultiProblem(ElementwiseProblem):
    def __init__(self, args, cfg, log):
        super().__init__(n_var=1, n_obj=3, n_constr=0)
        self.args = args
        self.cfg = cfg
        self.log = log

    def _evaluate(self, x, out, *args, **kwargs):
        F1 = self.func_1(x)
        F2 = self.func_2(x)
        F3 = self.func_3(x) # 这里写的理论产量计算是 相反数

        # 在评估的地方保存个体和适应度值
        with open(self.args.save_path+'all_x_f.txt', 'a') as f:
            # f.write(str([F1, F2, F3]))
            f.write(str([F1, F2, -F3])) # f3 是理论产量的负值
            f.write('\n')
            f.write(str(x[0].chrom.get_rxn_list()))
            f.write('\n')

            travel_list = x[0].chrom.travel()
            for line in travel_list:
                for l in line:
                    f.write(l)
                    f.write('\t')
                f.write('\n')
            f.write('==' * 10)
            f.write('\n')

        out["F"] = np.array([F1, F2, F3], dtype=float)

    def func_1(self, Chrom):
        # length of the pathway
        # less is better
        self.log.logger.info('---get pathway lenght---')
        return Chrom[0].chrom.length() #  chrom是一个SingleLinkList对象，包含length函数，可以计算反应个数

    def func_2(self, Chrom):
        # Gibbs of the pathway MDF
        # less is better
        self.log.logger.info('---get pathway gibbs---')
        return get_gibbs(self.cfg, Chrom[0].chrom.get_rxn_list()) # get_rxn_list 是一个只有Rxxxxx的反应列表

    def func_3(self, Chrom):
        # yelid of pathway
        # less is better
        rxn_list = Chrom[0].chrom.get_rxn_list()
        self.log.logger.info('---get pathway yield---')
        self.log.logger.info(Chrom[0].chrom.travel())
        self.log.logger.info(rxn_list)

        chrom = Chrom[0].chrom

        # return -get_yield(self.cfg, rxn_list)
        return -get_yield(self.cfg,chrom,self.log,self.args)

        # return -np.random.rand()



