#!/bin/bash 
# multi-
nohup python tests/main.py --task_id 202305081947 --gen 10 > tests/05081747.log 2>&1 & 
nohup python tests/main.py --task_id 202305081948 --gen 15 > tests/05081948.log 2>&1 & 
nohup python tests/main.py --task_id 202305081949 --gen 20 > tests/05081949.log 2>&1 & 

# # single-
nohup python tests/main.py --task_id 20230508190110 --gen 10 --algorithm single -w1 0.1 -w2 0.4 -w3 0.5 > tests/0508650110.log 2>&1 & 
nohup python tests/main.py --task_id 20230508190210 --gen 10 --algorithm single -w1 0.2 -w2 0.3 -w3 0.5 > tests/0508650210.log 2>&1 & 
nohup python tests/main.py --task_id 20230508190310 --gen 10 --algorithm single -w1 0.33 -w2 0.33 -w3 0.33 > tests/05086560310.log 2>&1 & 

# # single-
nohup python tests/main.py --task_id 20230508190115 --gen 15 --algorithm single -w1 0.1 -w2 0.4 -w3 0.5 > tests/0508650115.log 2>&1 & 
nohup python tests/main.py --task_id 20230508190215 --gen 15 --algorithm single -w1 0.2 -w2 0.3 -w3 0.5 > tests/0508650215.log 2>&1 & 
nohup python tests/main.py --task_id 20230508190315 --gen 15 --algorithm single -w1 0.33 -w2 0.33 -w3 0.33 > tests/05086560315.log 2>&1 & 

# # single-
nohup python tests/main.py --task_id 20230508190120 --gen 20 --algorithm single -w1 0.1 -w2 0.4 -w3 0.5 > tests/0508650120.log 2>&1 & 
nohup python tests/main.py --task_id 20230508190220 --gen 20 --algorithm single -w1 0.2 -w2 0.3 -w3 0.5 > tests/0508650220.log 2>&1 & 
nohup python tests/main.py --task_id 20230508190320 --gen 20 --algorithm single -w1 0.33 -w2 0.33 -w3 0.33 > tests/05086560320.log 2>&1 & 