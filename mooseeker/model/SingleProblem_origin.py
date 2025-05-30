"""
Autor: caoyh
Date: 2022-11-17 16:37:08
LastEditTime: 2022-11-20 22:01:19
"""
import numpy as np
from pymoo.core.problem import ElementwiseProblem
from mooseeker.utils import SingleLinkList, SingleReaction
from mooseeker.fit_funcs import get_gibbs
from mooseeker.fit_funcs import get_yield

class SingleProblem(ElementwiseProblem):
    def __init__(self, args, cfg, log):
        super().__init__(n_var=1, n_obj=1, n_constr=0)
        self.args = args
        self.cfg = cfg
        self.log = log

    def _evaluate(self, x, out, *args, **kwargs):
        F1 = self.func_1(x)
        F2 = self.func_2(x)
        F3 = self.func_3(x)

        # 在评估的地方保存个体和适应度值
        with open(self.args.save_path+'all_x_f.txt', 'a') as f:
            f.write(str([F1, F2, -F3]))
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

        ## 适应度值
        fit = float(self.args.weight3) * F1 + float(self.args.weight3) * F2 + float(self.args.weight3) * F3
        out["F"] = np.array([fit], dtype=float)

    def func_1(self, Chrom):
        # length of the pathway
        # less is better

        self.log.logger.info('---get pathway lenght---')
        return Chrom[0].chrom.length()

    def func_2(self, Chrom):
        # Gibbs of the pathway MDF
        # less is better

        self.log.logger.info('---get pathway gibbs---')
        return get_gibbs(self.cfg, Chrom[0].chrom.get_rxn_list())

    def func_3(self, Chrom):
        # yelid of pathway
        # less is better
        rxn_list = Chrom[0].chrom.get_rxn_list()
        self.log.logger.info('---get pathway yield---')
        self.log.logger.info(Chrom[0].chrom.travel())
        self.log.logger.info(rxn_list)
        return -get_yield(self.cfg, rxn_list)
        # return -np.random.rand()


