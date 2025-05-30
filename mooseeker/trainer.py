# -*- coding:utf-8 -*-
# @Author     : caoyh
# @Time       : 2023/2/24 11:28
# Description :
import sys
import os

script_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))

sys.path.append(script_dir)

import time
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.optimize import minimize
from mooseeker.utils import get_config
from mooseeker.model import MultiProblem
from mooseeker.model import SingleProblem
from mooseeker.model.BiOperators import *
from mooseeker.log import Logger
import numpy as np
import dill

# 根据request_dict构建的args，在utils中
def trainer(args,save_checkpoint=True,save_history=True):

    cfg = get_config()

    log = Logger(args.save_path + time.strftime("%Y%m%d") + '.log')

    log.logger.info('=====> Train for %s<=====' % args.project)
    start_time = time.time()
    # if args.algorithm == 'multi':
    # algorithm = NSGA2(pop_size=args.pop_size, # 2
    #                 sampling=BioSampling(),
    #                 crossover=BioCrossover(),
    #                 mutation=BioMutation(),
    #                 eliminate_duplicates=BioDuplicateElimination()) # 这四个Bio函数都是在BIOperators中的

    ## 这里用到的pymoo包是一个优化算法包，里面包含了遗传算法GA，多目标优化算法NSGA2
    # res = minimize(problem=MultiProblem(args, cfg, log),
    #             algorithm=algorithm,
    #             termination=('n_gen', args.gen), # 终止标准
    #             seed=1,
    #             callback=BioCallback(), # 在BIOperators中
    #             verbose=True,
    #             save_history=True)

    # if args.algorithm == 'single':
    algorithm = GA(pop_size=args.pop_size, #cfg['init']['pop_size'],
                   sampling=BioSampling(),
                   crossover=BioCrossover(),
                   mutation=BioMutation(),
                   eliminate_duplicates=BioDuplicateElimination())

    res = minimize(problem=SingleProblem(args, cfg, log),
                   algorithm=algorithm,
                   # termination=('n_gen', cfg['init']['gen']),
                   termination=('n_gen', args.gen),
                   seed=1,
                   callback=BioCallback(),
                   verbose=True,
                   save_history=True)

    # save checkpoint
    if save_checkpoint:
        checkpoint_file = args.save_path + time.strftime("%Y%m%d")
        log.logger.info('=====> save checkpoint to file:  %s<=====' % (checkpoint_file))
        with open(checkpoint_file, "wb") as f:
            dill.dump(algorithm, f) # dill可以认为是增强版的pickle dill.dump(obj, file)  obj 要保存的变量对象   file 要保存至的文件

    # ------------------ save result ------------------

    # save history
    if save_history:
        hist = res.history
        n_evals = []  # corresponding number of function evaluations\
        hist_X = []  # the individuals in each generation
        hist_F = []  # the objective space values in each generation
        hist_cv = []  # constraint violation in each generation
        hist_cv_avg = []  # average constraint violation in the whole population

        for algo in hist:
            # store the number of function evaluations
            n_evals.append(algo.evaluator.n_eval)

            # retrieve the optimum from the algorithm
            opt = algo.opt

            # store the least contraint violation and the average in each population
            hist_cv.append(opt.get("CV").min())
            hist_cv_avg.append(algo.pop.get("CV").mean())

            # get all the individuals
            hist_X.append(opt.get("X").tolist())

            # filter out only the feasible and append and objective space values
            feas = np.where(opt.get("feasible"))[0]
            hist_F.append(opt.get("F")[feas].tolist())

            # save to file
        log.logger.info('=====> Save result and history ! <=====')
        np.savez(args.save_path + 'result.npz', X=res.X, F=res.F, CV=res.CV, G=res.G, dtype=object)
        np.savez(args.save_path + 'history.npz', n_evals=n_evals, X=hist_X, F=hist_F, cv=hist_cv, cv_avg=hist_cv_avg,
                 dtype=object)
    time_cost = time.time() - start_time

    with open(args.save_path + 'result.txt', 'w') as f:
        f.write(f"This is for {args.project}.\n")
        f.write(f"task_id={args.task_id}\n")
        # f.write(f"algorithm type is {args.algorithm}\n")
        f.write(f"population size={args.pop_size}\n")
        f.write(f"gen={args.gen}\n")
        f.write(f"The computational cost is {time_cost}s = {time_cost/60}min = {time_cost/3600}h.\n")

    log.logger.info('=====Done!=====')

    # algorithm: multi 多目标算法会有多个结果，所以需要使用两次引用才能得到最后的结果
    # res_dict = {
    #     'task_id': args.task_id,
    #     'result': {
    #         'best_path_rxn_list': res.X[0][0].chrom.get_rxn_list(),
    #         'best_path_cpd_list': res.X[0][0].chrom.get_cpd_list(),
    #         'best_path_eva': {
    #             'length': res.F[0][0],
    #             'gibbs': res.F[0][1],
    #             'yield': res.F[0][2],
    #         }
    #     }
    # }


    # 单目标算法 只使用理论产量作为目标值
    # 对于 理论产量 应该输出其相反数 才是真实值。
    print("----------- res.X[0].chrom.get_rxn_list -----------")
    print(res.X[0].chrom.get_rxn_list)

    print("----------- res.F -----------")
    print(np.shape(res.F))
    print(res.F[0])


    res_dict = {
        'task_id': args.task_id,
        'result': {
            'best_path_rxn_list': res.X[0].chrom.get_rxn_list(),
            'best_path_cpd_list': res.X[0].chrom.get_cpd_list(),
            'best_path_eva': {
                # 'length': res.F[0][0],
                # 'gibbs': res.F[0][1],
                # 'yield': -res.F[0][2],
                'yield': -res.F[0],
            }
        }
    }

    # single：单目标算法最后的结果只有一个，所以最后只需要引用一次即可
    # res_dict = {
    #     'task_id': args.task_id,
    #     'result': {
    #         'best_path_rxn_list': res.X[0][0].chrom.get_rxn_list(),
    #         'best_path_cpd_list': res.X[0][0].chrom.get_cpd_list(),
    #         'best_path_eva': {
    #             'length': res.F[0][0],
    #             'gibbs': res.F[0][1],
    #             'yield': res.F[0][2],
    #         }
    #     }
    # }

    return res_dict 