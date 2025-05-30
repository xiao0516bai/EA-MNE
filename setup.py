# -*- coding:utf-8 -*-
# @Author     : caoyh
# @Time       : 2023/2/25 17:05
# Description :
import setuptools
from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="mooseeker",
    # version="0.0.1",
    version="0.0.2", # for project
    author="caoyh",
    author_email="caoyahui@tju.edu.cn",
    description="MooSeeker package",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://gitee.com/bellacaoyh/newmooseeker.git",
    packages=setuptools.find_packages(),
    license='MIT',
    install_requires=['pymoo','cobra','bs4','requests','dill','streamlit','pyyaml','pandas','ipykernel'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,  # 打包包含静态文件标识
)
