<!--
 * @Author: xiao0516bai
 * @Date: 2025-5-30 09:27:27
 * @LastEditTime: 2022-5-30 10:19:53
-->
# EA_MNE
## Installion

```bash
conda create -n ea python=3.7 -y
conda activate ea
conda conda install -c rdkit rdkit -y
pip install -r requirement.txt
```

## Introduction
```
Within the EA_MNE method, we propose two new evaluation criteria: 

    (1) the number of effective branching reactions, which assesses the extent of branching impacts, 
    (2) the network theoretical yield, which precisely quantifies yield losses caused by branching reactions.

Furthermore, we develop a comprehensive metabolic pathway design method based on evolutionary algorithms integrating four key criteria: the number of effective branching reactions, network
theoretical yield, network toxicity, and Gibbs free energy. This integrated approach provides a systematic solution for addressing branching reaction challenges, significantly improving both the accuracy of pathway evaluation and the synthetic efficiency of microbial systems.
```

## KEGG Database

> The KEGG database is an open access database that can download through [KEGG API](https://www.kegg.jp/kegg/rest/keggapi.html).
We provide `mooseeker/kegg_helper/kegg.py` to get necessary information from KEGG.

## Usage

```bash
python sfla.py
```
Here, a set of paramters can be provided as `dict` or through `argparse`.

The list of parameters is shown as follows:
```bash
'--task_id', 'The unique id for task.'
'--NrMax', 'The max number of the pathway length.'    
'--ob_substrate', 'The objective substrate for production.'
'--ob_product', 'The objective product for production.'   
'--abundant', 'The available compounds can be chosen at the initial step.'
'--database', 'The available database.'
'--host', 'The host cell aims to calculate yield for the pathway.'
```