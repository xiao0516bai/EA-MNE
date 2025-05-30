# -*- coding:utf-8 -*-
# @Author     : caoyh
# @Time       : 2023/2/24 11:45
# Description :
import argparse
import os


def parser_args(params={}):
    # args = argparse.Namespace(**params)

    parser = argparse.ArgumentParser()
    parser.add_argument('--task_id', type=str, default='001', help='The unique id for task.')
    parser.add_argument('--algorithm', type=str, default='single', help='The type of the optimizer.')
    parser.add_argument('--pop_size', default=20, type=int, help='The size of population.')
    parser.add_argument('--gen', default=20, type=int,
                        help='The generation number to evolve in evaluationary algorithm.')
    parser.add_argument('--NrMax', default=10, type=int, help='The max number of the pathway length.')
    parser.add_argument('--ob_substrate', default='C00082', help='The objective substrate for production.')
    parser.add_argument('--ob_product', default='C00755', help='The objective product for production.')
    parser.add_argument('--abundant',
                        default=["C00001", "C00002", "C00080"],
                        nargs='+', help='The available compounds can be chosen at the initial step.')
    parser.add_argument('--database', default=['KEGG', 'BiGG'], nargs='+', help='The available database.')
    parser.add_argument('--host', default='ecoli', help='The host cell aims to calculate yield for the pathway.')
    parser.add_argument('--eva_func', default=['total_reactions', 'net_tox', 'network_yield', 'gibbs'], nargs='+',
                        help='The evaluation function to evaluate pathway.')
    parser.add_argument('-w1', '--weight1', default=0.4, type=float,
                        help='The weight of length for single objective algorithm')
    parser.add_argument('-w2', '--weight2', default=0.2, type=float,
                        help='The weight of gibbs for single objective algorithm')
    parser.add_argument('-w3', '--weight3', default=0.4, type=float,
                        help='The weight of yield for single objective algorithm')
    parser.add_argument('-w4', '--weight4', default=0.2, type=float,
                        help='The weight of yield for single objective algorithm')
    parser.add_argument('--save_path', default='./tests/result0508', help='The dir to save the result.')

    parser.set_defaults(**params) #如果传入默认参数，则设置，否者可以使用命令行输入
    args = parser.parse_args()

    args.project = '_'.join((args.task_id, args.algorithm, args.ob_product, str(args.pop_size), str(args.gen), str(int(args.weight1*100))))
    args.save_path = '/'.join((args.save_path, args.project)) + '/'
    if not os.path.exists(args.save_path):
        os.makedirs(args.save_path)
        # os.mkdir(args.save_path)

    return args
