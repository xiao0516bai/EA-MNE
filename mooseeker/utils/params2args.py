# -*- coding:utf-8 -*-
# @Author     : caoyh
# @Time       : 2023/3/3 09:57
"""
django cannot process argparse
so let this file to get the params to args
"""
import argparse
import os

def params2args(params):
    args = argparse.Namespace(**params)
    args.w1 = 0.4 
    args.w2 = 0.2
    args.w3 = 0.4
    args.w3 = 0.2
    args.save_path = "."
    args.project = '_'.join((args.task_id, args.ob_substrate, args.ob_product))
    args.save_path = '/'.join((args.save_path, args.project)) + '/'
    if not os.path.exists(args.save_path): os.mkdir(args.save_path)
    return args


