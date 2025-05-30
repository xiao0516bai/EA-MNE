'''
Description: 
Autor: caoyh
Date: 2022-11-16 12:59:17
LastEditTime: 2022-11-20 16:06:27
'''
import os
import json
import pdb
from rdkit import DataStructs
from rdkit import Chem
import re, time
import numpy as np
import pandas as pd
import sys
# sys.path.append('.')

script_dir = os.path.abspath(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))
sys.path.append(script_dir)
# print(script_dir)

from mooseeker.utils import get_config
from mooseeker.utils import read_json
from mooseeker.utils import SingleReaction
from mooseeker.kegg_helper.pykegg import split_equation
   
## 第一步 基于KEGG数据库 初始化pair池子
class MyPool(object):
    def __init__(self, cfg) -> None:
        self.cfg = cfg
        self.rxn_dict = read_json(cfg['file_path']['all_rxn_dict'])
        self.cpd_dict = read_json(cfg['file_path']['all_cpd_dict'])

        self.save_path = self.cfg['file_path']['data_dir']
        if not os.path.exists(self.save_path): os.mkdir(self.save_path) # /cache/data
        
        self.pool = dict()
        self._get_pool()

    def get_Tanimoto(self, s,p):
        def check_cpd(cpd_name):
            # 1. check cpd name
            pat = re.compile('^C\d{5}')
            res = pat.findall(cpd_name)
            if res != []:
                return True
            else:
                return False
        assert (check_cpd(s) and check_cpd(p))
        
        # 2. check smile
        if self.cpd_dict[s]['smile']!=None and self.cpd_dict[p]['smile']!=None:
            try:
                # 检查smile格式是否存在 存在的话 根据smile格式获得分子的mol文件
                s_smile = self.cpd_dict[s]['smile'] # 这里就是，从反应数据库拿到了kegg化合物id，然后拿着id，去化合物数据库中 取smile
                s_mol = Chem.MolFromSmiles(s_smile)
                p_smile = self.cpd_dict[p]['smile']
                p_mol = Chem.MolFromSmiles(p_smile)
                mols = [s_mol, p_mol] # 这里存的就是两个化合物 底物——产物 对
                fps = [Chem.RDKFingerprint(x) for x in mols]
                t = DataStructs.FingerprintSimilarity(fps[0], fps[1]) # 计算相似性
            except Exception as e:
                print(e)
                print("s is %s and p is %s" % (s, p))
                t = 0
        else:
            t = 0
        print('Tanimoto of %s and %s is %s.' %(s, p, str(t)))
        return t

    # 无效的化合物列表？
    def _get_invalid_cpd_list(self):
        # 读取化合物csv
        def get_cpd_csv(file_path):
            cpds = pd.read_csv(file_path, header=None)
            cpds = cpds.drop_duplicates() # 删除重复数据
            cpds_list = cpds[0].values.tolist()
            return cpds_list
        cofactors_list = get_cpd_csv(self.save_path+'cofactors.csv')
        error_cpds_list = get_cpd_csv(self.save_path+'error_cpds.csv')
        return cofactors_list+error_cpds_list 

    # 从数据库得到pair池子 保存至MYPOOL_20221120.npy里面

    def _get_pool(self):
        invalid_cpds = self._get_invalid_cpd_list() # 返回cache/data下，cofactors.csv和error_cpds.csv合并后的文件

        KEYS = list(self.rxn_dict.keys()) # 得到kegg反应id
        len_pool = 0
        for R in KEYS:
            s_list, p_list = split_equation(self.rxn_dict[R]["equation"]) # 获取每个反应的方程式 然后再将其划分为反应物和生成物列表

            ## kegg中的反应都是可逆的 所以要分别以反应物和产物为起始，进行循环。 得到 反应物——底物对 不包含辅因子所在的对。
            # 这里的s和p就是单独的 一个化合物了
            for s in s_list:
                for p in p_list:
                    if (s in invalid_cpds) or (p in invalid_cpds): continue
                    srn_front = SingleReaction(s, p, self.rxn_dict[R], Tanimoto=self.get_Tanimoto(s,p)) # Tanimoto=self.get_Tanimoto(s,p) 计算反应物和生成物之间的相似性。
                    # 如果键不存在于字典中，将会添加该键并将default的值设为该键的默认值，如果键存在于字典中，将读出该键原来对应的值，default的值不会覆盖原来已经存在的键的值。
                    # 其实就是往字典中添加键值对。
                    # 如果pool中有没有键s则，生成默认值set()空集合，并将此srn_front对象添加到集合中
                    # 若pool中有s键，则直接将srn_front对象添加到集合中
                    # 也就是说，每个反应物，对应多个SingleReaction对象，也就是对应多个 反应物——生成物 对
                    self.pool.setdefault(s, set()).add(srn_front)
                    len_pool += 1
            for p in p_list:
                for s in s_list:
                    if (s in invalid_cpds) or (p in invalid_cpds): continue
                    srn_front = SingleReaction(p, s, self.rxn_dict[R], Tanimoto=self.get_Tanimoto(p,s))
                    self.pool.setdefault(p, set()).add(srn_front)
                    len_pool += 1
            print("Reaction %s is done!" %(R))

        pool_all_file = self.save_path + 'MYPOOL_' + time.strftime("%Y%m%d") + '.npy'   
        np.save(pool_all_file, self.pool) 

        print('===== Done! ======')
        print("There are %s pairs in the pool" %(len_pool))

def main():
    cfg = get_config() 
    MyPool(cfg)
    

if __name__=='__main__':
    main()
    

