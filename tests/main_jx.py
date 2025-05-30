'''
Autor: caoyh
Date: 2022-11-18 09:37:01
LastEditTime: 2022-11-20 22:06:09
'''
import os
import sys
import dill
import time
import numpy as np
import argparse
from tabnanny import verbose
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.optimize import minimize
from argparse import ArgumentParser

# sys.path.append('/home/jiax/anaconda3/chile/pathway-search')

moo_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
print("moo_dir:",moo_dir)


sys.path.append(moo_dir) #添加project路径
sys.path.append(moo_dir + "/mooseeker") #添加project路径下的mooseeker module路径
# sys.path.append(moo_dir + "/utils")
# sys.path.append(moo_dir + "/log")
# 这个是我自己加的 不然运行nohup命令的时候运行不了
sys.path.append("..")

print("sys.path:",sys.path)


from mooseeker.utils import get_config
# from model import *
from mooseeker.model import MultiProblem
from mooseeker.model import SingleProblem
from mooseeker.log import Logger
from mooseeker.model.BiOperators import *

# def parser():
#     parser = ArgumentParser('MooSeeker')
#     parser.add_argument('-alg', '--algorithm', default='single',
#                         help='Multi or Single objective algorithm to use')  ## multi OR single
#     parser.add_argument('-p', '--production', default='vanillin',
#                         help='Glycolysis or Vanillin to produce')  # glycolysis or vanillin
#     parser.add_argument('-w1', '--weight1', default=0, type=float,
#                         help='The weight of length for single objective algorithm')
#     parser.add_argument('-w2', '--weight2', default=0.3, type=float,
#                         help='The weight of gibbs for single objective algorithm')
#     parser.add_argument('-w3', '--weight3', default=0.7, type=float,
#                         help='The weight of yield for single objective algorithm')
#     args = parser.parse_args(args=[])
#
#     if args.algorithm == 'multi':
#         args.project = '_'.join((args.algorithm, args.production))
#     elif args.algorithm == 'single':
#         args.project = '_'.join(('1027_GA', args.production, str(int(args.weight1 * 100)), str(int(args.weight2 * 100)),
#                                  str(int(args.weight3 * 100))))
#     return args

# parser 是jx师姐论文中添加的函数 师姐论文中主要适用getconfig 也就是config.yml来更改参数
def parser():
    params = {
        "task_id": "010",
        "algorithm": "single",
        "pop_size": 10,  # 原本是10
        "gen": 20,  # 原本是20
        "NrMax": 15,
        "ob_substrate": "C00267",
        "ob_product": "C00022",
        # "ob_substrate": "C00082",
        # "ob_product": "C00755", # 香兰素
        "abundant": ["C00001", "C00002", "C00080"],
        "database": ["KEGG", "BiGG"],
        "host": ["ecoli"],
        "eva_func": ["yield"]
        # "eva_func": ["length", "gibbs", "yield"]
    }

    args = argparse.Namespace(**params)
    # args.w1 = 0
    # args.w2 = 0.3
    # args.w3 = 0.7

    # 因为我在singleproblem中把F1 和 F2注释掉了，相当于只计算理论产量。所以这里不传入w1和w2 没有影响
    # args.weight1 = 0
    # args.weight2 = 0.3
    args.weight3 = 0.7
    args.save_path = ".."

    args.project = '_'.join((args.task_id, args.ob_substrate, args.ob_product))
    args.save_path = '/'.join((args.save_path, args.project)) + '/'
    print("main_jx:args.save_path")
    print(args.save_path)
    if not os.path.exists(args.save_path): os.mkdir(args.save_path)

    # parser = ArgumentParser('MooSeeker')
    # parser.add_argument('-alg', '--algorithm', default='single',
    #                     help='Multi or Single objective algorithm to use')  ## multi OR single
    # parser.add_argument('-p', '--production', default='vanillin',
    #                     help='Glycolysis or Vanillin to produce')  # glycolysis or vanillin
    # parser.add_argument('-w1', '--weight1', default=0, type=float,
    #                     help='The weight of length for single objective algorithm')
    # parser.add_argument('-w2', '--weight2', default=0.3, type=float,
    #                     help='The weight of gibbs for single objective algorithm')
    # parser.add_argument('-w3', '--weight3', default=0.7, type=float,
    #                     help='The weight of yield for single objective algorithm')
    # args = parser.parse_args(args=[])
    #
    # if args.algorithm == 'multi':
    #     args.project = '_'.join((args.algorithm, args.production))
    # elif args.algorithm == 'single':
    #     args.project = '_'.join(('1027_GA', args.production, str(int(args.weight1 * 100)), str(int(args.weight2 * 100)),
    #                              str(int(args.weight3 * 100))))
    return args

def train(args, cfg, log):
    log.logger.info('=====> Train for %s<=====' % (args.project))
    result_dir = cfg['file_path']['result_dir'] + args.project + '/'
    if not os.path.exists(result_dir): os.mkdir(result_dir)
    if args.algorithm == 'multi':
        algorithm = NSGA2(pop_size=cfg['init']['pop_size'],
                          sampling=BioSampling(),
                          crossover=BioCrossover(),
                          mutation=BioMutation(),
                          eliminate_duplicates=BioDuplicateElimination())

        res = minimize(problem=MultiProblem(args, cfg, log),
                       algorithm=algorithm,
                       termination=('n_gen', cfg['init']['gen']),
                       seed=1,
                       callback=BioCallback(),
                       verbose=True,
                       save_history=True)

    # 师姐用的是单目标算法
    if args.algorithm == 'single':
        algorithm = GA(pop_size=cfg['init']['pop_size'],
                       sampling=BioSampling(),

                       crossover=BioCrossover(),
                       mutation=BioMutation(),
                       eliminate_duplicates=BioDuplicateElimination())

        res = minimize(problem=SingleProblem(args, cfg, log),
                       algorithm=algorithm,
                       termination=('n_gen', cfg['init']['gen']),
                       seed=1,
                       callback=BioCallback(),
                       verbose=True,
                       save_history=True)

    # save checkpoint
    checkpoint_file = cfg['file_path']['checkpoint_dir'] + args.project + '_' + time.strftime("%Y%m%d")
    log.logger.info('=====> save checkpoint to file:  %s<=====' % (checkpoint_file))
    with open(checkpoint_file, "wb") as f:
        dill.dump(algorithm, f)

    # print("Best solution found: \nF = %s\nX = %s" % (res.F, res.X[0][0].chrom.travel()))
    # print(res.X[0][0].chrom.get_rxn_list())

    # ------------------ save result ------------------

    # save history
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
    np.savez(result_dir + 'result.npz', X=res.X, F=res.F, CV=res.CV, G=res.G, dtype=object)
    np.savez(result_dir + 'history.npz', n_evals=n_evals, X=hist_X, F=hist_F, cv=hist_cv, cv_avg=hist_cv_avg,
             dtype=object)
    log.logger.info('=====Done!=====')


if __name__ == '__main__':
    old_time = time.time()
    print("运行时间" + str(old_time))
    args = parser()
    # print(args.project)
    cfg = get_config()
    log = Logger(cfg['file_path']['log_dir'] + 'real_1027_GA.log', fmt='%(message)s')
    train(args, cfg, log)
    current_time = time.time()
    print("运行时间为" + str(current_time - old_time) + "s")
    result_dir = cfg['file_path']['result_dir'] + args.project + '/'
    with open(result_dir + 'time_single_van_ga.txt', 'a') as file:
        file.write("运行时间为" + str(current_time - old_time) + "s")
        file.write('\n')