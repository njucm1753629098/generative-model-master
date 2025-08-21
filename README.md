## Introduction
A novel AI algorithm of TransGenGRU is proposed to perform de novo natural products design that employed as novel NLRP3 inhibitors. The code in this repository is inspired on [QBMG](https://github.com/SYSU-RCDD/QBMG)
## Usage
### 1.perform pre-training process
`python train.py`
### 2.perform transfer learning process
`python transfer_learning.py`
### 3.generate molecules with the pretrained model
`python sample.py data/50_epoch.ckpt`
### 4.generate molecules with the transfer learning model
`python sample.py data/200_epochs_transfer_augument3_one.ckpt`
### 5.calculate the properties of the pre-training dataset, transfer learning dataset, and generated molecules
`python calculate_property.py`
### 6.compare the properties between the transfer learning dataset and generated molecules
`python pca_plot.py`
`python KDE.py`


