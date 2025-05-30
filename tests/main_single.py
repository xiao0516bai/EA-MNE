import sys
import os
import argparse
import json
moo_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
print("moo_dir:",moo_dir)


sys.path.append(moo_dir) #添加project路径
sys.path.append(moo_dir + "/mooseeker") #添加project路径下的mooseeker module路径

# 这个是我自己加的 不然运行nohup命令的时候运行不了
sys.path.append("..")

print("sys.path:",sys.path)
from mooseeker.trainer import trainer
from argparse import ArgumentParser

# jx师姐添加的一个函数 在pso.py中有调用
def parser():
    params = {
        "task_id": "009",
        "algorithm": "single",
        "pop_size": 10,  # 原本是10
        "gen": 20,  # 原本是20
        "NrMax": 15,
        "ob_substrate": "C00267",
        "ob_product": "C00022",
        "abundant": ["C00001", "C00002", "C00080"],
        "database": ["KEGG", "BiGG"],
        "host": ["ecoli"],
        "eva_func": ["yield"]
        # "eva_func": ["length", "gibbs", "yield"]
    }

    args = argparse.Namespace(**params)
    args.w1 = 0
    args.w2 = 0.3
    args.w3 = 0.7
    # args.weight1 = 0
    # args.weight2 = 0
    # args.weight3 = 1
    args.save_path = ".."

    args.project = '_'.join((args.task_id, args.ob_substrate, args.ob_product))
    args.save_path = '/'.join((args.save_path, args.project)) + '/'
    print("main_single:args.save_path")
    print(args.save_path)

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

# from trainer import trainer
if __name__ == "__main__":
    print(moo_dir)
    ##########################################
    # 单独运行算法的方法
    # args = parser_args()
    # print(args.__dict__)
    # res = trainer(args)
    # print(res)
    ##########################################

    #########################################
    # 运行项目的方法
    # params = {
    #     "task_id": "005",
    #     # "algorithm": "single",
    #     "pop_size": 10,
    #     "gen": 20,
    #     "NrMax": 15,
    #     "ob_substrate": "C00082",
    #     "ob_product": "C00755",
    #     "abundant": ["C00001", "C00002", "C00080"],
    #     "database": ["KEGG", "BiGG"],
    #     "host": ["ecoli"],
    #     "eva_func": ["length", "gibbs", "yield"]
    # }

    # 糖酵解路径
    params = {
        "task_id": "009",
        "algorithm": "single",
        "pop_size": 10, # 原本是10
        "gen": 20,  # 原本是20
        "NrMax": 15,
        "ob_substrate": "C00267",
        "ob_product": "C00022",
        "abundant": ["C00001", "C00002", "C00080"],
        "database": ["KEGG", "BiGG"],
        "host": ["ecoli"],
        "eva_func": ["yield"]
        # "eva_func": ["length", "gibbs", "yield"]
    }

    args = argparse.Namespace(**params)
    # args.w1 = 0.33
    # args.w2 = 0.33
    # args.w3 = 0.33
    # args.weight1 = 0
    # args.weight2 = 0
    args.weight3 = 1
    args.save_path = ".."

    args.project = '_'.join((args.task_id, args.ob_substrate, args.ob_product))
    args.save_path = '/'.join((args.save_path, args.project)) + '/'
    print("main_single:args.save_path")
    print(args.save_path)
    # args.save_path = '/'.join((base_dir, 'static/result', args.project)) + '/'
    if not os.path.exists(args.save_path): os.mkdir(args.save_path)
    res = trainer(args)
    print(res)

    # 我想把这个res追加到result.txt里面
    with open(args.save_path + 'result.txt', 'a') as f:
        f.write(f"This is result")
        # f.write(res)
        f.write(json.dumps(res))

