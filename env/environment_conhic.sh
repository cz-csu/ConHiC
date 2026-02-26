#!/bin/bash

# ConHiC 环境安装脚本
# 用于安装 ConHiC 所需的依赖环境

echo "开始安装 ConHiC 环境..."

# 添加 conda-forge 通道
conda config --add channels conda-forge

# 创建 ConHiC 环境
conda create -n conhic python=3.10.14 -y

# 激活环境
conda activate conhic

# 安装基础依赖
conda install mkl -y
conda install sparse_dot_mkl "numpy<2.0.0" -y

# 安装 Python 包
pip3 install scipy matplotlib
pip3 install scikit-learn networkx "pysam==0.20.0"
pip3 install portion

# 导出环境配置
conda env export > environment_conhic.yml

echo "ConHiC 环境安装完成！"
echo "环境配置文件已保存为: environment_conhic.yml"