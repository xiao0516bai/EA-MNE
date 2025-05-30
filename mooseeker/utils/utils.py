'''
Description: 
Autor: caoyh
Date: 2022-11-16 11:30:14
LastEditTime: 2022-11-20 16:28:54
'''
# -*- coding: utf-8 -*-
# @Author: Caoyh
# @Date:   2022-06-15 09:56:12
# @Last Modified time: 2022-07-17 15:57:44
import sys
import os
import json
import csv
import pandas as pd
import numpy as np
import yaml

def write_json(data, file):
    with open(file, 'w') as f:
        json.dump(data, f)

def read_json(file):
    with open(file, 'r') as f:
        data = json.load(f)
    return data

## 获取参数
def get_config():
    # project_dir = os.path.abspath(os.path.dirname(__file__)) 获得当前目录
    project_dir = os.path.abspath(os.path.dirname(os.path.dirname(__file__)))  # 上级目录
    yaml_file = os.path.join(project_dir, 'data/config.yaml')
    with open(yaml_file, 'r', encoding='utf-8') as f:
        data_dict = yaml.load(f, Loader=yaml.FullLoader)

    for key in data_dict['file_path']:
        data_dict['file_path'][key] = os.path.join(project_dir, data_dict['file_path'][key])

    ## 存储酶列表——反应列表——网络总反应长度的文件 每次运行都先删除之前的备份
    if os.path.exists(data_dict['file_path']['enzyme_rxn_count_dir']):  # 检查文件是否存在
        os.remove(data_dict['file_path']['enzyme_rxn_count_dir'])  # 如果存在，删除文件

    return data_dict

class SingleReaction(object):
    def __init__(self, S, P, reaction, Tanimoto=0.0, next=None):
        # S 是底物 P是产物
        self.S = S
        self.P = P
        self.reaction = reaction
        self.Tanimoto = Tanimoto
        self.next = next

class SingleLinkList(object):
    def __init__(self):
        self._head = None

    def is_empty(self):
        return self._head == None

    def clear(self):
        self._head = None

    ## 得到反应字典
    def get_rxn_dict(self):
        rxn_dict = {}
        cur = self._head 
        while cur != None:
            # 用字典存储这条链表上的反应
            # 键：kegg反应id
            # 值：方程式
            rxn_dict[cur.reaction['kegg_id']] = cur.reaction['equation']
            cur = cur.next
        return rxn_dict

    # 得到反应列表
    def get_rxn_list(self):
        rxn_list = []
        cur = self._head 
        while cur != None:
            rxn_list.append(cur.reaction['kegg_id'])
            cur = cur.next
        return rxn_list

    # 得到化合物列表
    def get_cpd_list(self):
        cpd_list = []
        cur = self._head
        while cur != None:
            if cur.S not in cpd_list:
                cpd_list.append(cur.S)
            if cur.P not in cpd_list:
                cpd_list.append(cur.P)
            cur = cur.next
        return cpd_list

    # [kegg反应id，反应物，生成物]
    def travel(self):
        """遍历整个链表"""
        travel_list = []
        cur = self._head
        while cur != None:
            travel_list.append([cur.reaction['kegg_id'], cur.S, cur.P])
            # print(cur.reaction['kegg_id'], ":", cur.S, " --> ", cur.P, ",")
            cur = cur.next
        return travel_list

    # 计算路径长度
    def length(self):
        count = 0
        cur = self._head
        while cur != None:
            count += 1
            cur = cur.next
        return count

    def add(self, node):
        """在头部添加节点，节点是已知，更新头节点"""
        # node = SingleNode()
        # node = Node(item)
        # 将node节点的next指向头节点。然后将头节点再设置为node即可。
        node.next = self._head
        self._head = node

    def get_first_node(self):
        return self._head

    # 先把指针指向第一个节点 然后不断next直到下一个next为null时停止。 并且将next指向需要添加的节点node。
    def append(self, node):
        """在尾部添加节点"""
        if self.is_empty():
            self._head = node
        else:
            cur = self._head
            # node = Node(item)
            while cur.next != None:
                cur = cur.next
            cur.next = node

    def insert(self, pos, node):
        """在选定的位置添加节点"""
        cur = self._head
        # node = Node(item)
        count = 0
        if pos <= 0:
            self.add(node)
        elif pos > (self.length() - 1):
            self.append(node)
        else:
            while count < (pos - 1):
                count += 1
                cur = cur.next
            node.next = cur.next
            cur.next = node

    def trim_chrom(self, input_chrom):
        chrom = copy.deepcopy(input_chrom)

        head = chrom._head
        # 修剪情况: a->b b->c c->a a->d ===> a->d
        if (chrom != None and chrom.length() > 2):
            dummy = SingleReaction(S=None, P=None, reaction={})
            dummy.next = head
            pre = dummy
            left = pre.next
            while(left != None):
                right = left.next
                while(right != None):
                    if left.S == right.S:
                        pre.next = right
                        # 自己修改的部分 原本的无法修改开头冗余的路径
                        if (pre.S == None):
                            chrom._head = pre.next
                        # pre = dummy
                        left = pre.next
                        # break  # 去掉break 因为有可能重复两次 a->b->a->b->a->c->d
                    right = right.next
                pre = pre.next
                left = left.next
                pre.next = left  # 添加的一行
            # true表示被剪枝。false表示没有被剪枝
            return chrom,True
        else:
            return input_chrom,False

