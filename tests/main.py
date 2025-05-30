import sys
import os
import argparse
moo_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))
print("moo_dir:",moo_dir)


sys.path.append(moo_dir) #添加project路径
sys.path.append(moo_dir + "/mooseeker") #添加project路径下的mooseeker module路径

# 这个是我自己加的 不然运行nohup命令的时候运行不了
sys.path.append("..")

print("sys.path:",sys.path)
from mooseeker.trainer import trainer
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
        "task_id": "005",
        # "algorithm": "single",
        "pop_size": 10,
        "gen": 20,
        "NrMax": 15,
        "ob_substrate": "C00267",
        "ob_product": "C00022",
        "abundant": ["C00001", "C00002", "C00080"],
        "database": ["KEGG", "BiGG"],
        "host": ["ecoli"],
        "eva_func": ["length", "gibbs", "yield"]
    }

    args = argparse.Namespace(**params)
    # args.w1 = 0.33
    # args.w2 = 0.33
    # args.w3 = 0.33
    args.weight1 = 0.33
    args.weight2 = 0.33
    args.weight3 = 0.33
    args.save_path = ".."

    args.project = '_'.join((args.task_id, args.ob_substrate, args.ob_product))
    args.save_path = '/'.join((args.save_path, args.project)) + '/'
    # args.save_path = '/'.join((base_dir, 'static/result', args.project)) + '/'
    if not os.path.exists(args.save_path): os.mkdir(args.save_path)
    res = trainer(args)
    print(res)
