#!/usr/bin/env python

import torch
from torch.utils.data import DataLoader
from rdkit import Chem
from rdkit import rdBase
import data_struct as ds
from data_struct import MolData, Vocabulary
from model import RNN
import sys
import random
import numpy as np


def Sample(filename):
    
    voc = Vocabulary(init_from_file="./Voc")  
    Prior = RNN(voc)
    # Can restore from a saved RNN
    Prior.rnn.load_state_dict(torch.load(filename))
    totalsmiles = []
     
    molecules_total = 0
    
    for epoch in range(1, 11):
        seqs, likelihood, _ = Prior.sample(100)
        valid = 0
        for i, seq in enumerate(seqs.cpu().numpy()):
            smile = voc.decode(seq)
            if Chem.MolFromSmiles(smile):
                valid += 1
                totalsmiles.append(smile)
                       
        molecules_total = len(totalsmiles)
        print(("\n{:>4.1f}% valid SMILES".format(100 * valid /len(seqs))))
        print(valid, molecules_total, epoch)
    print("\n{:>4.1f}% valid SMILES".format(molecules_total * 100/1000))
    return totalsmiles
    
def set_seed(seed=10):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


if __name__ == "__main__":
    set_seed(4)
    filename = sys.argv[1]
    
    print(filename)
    totalsmiles=Sample(filename)
    f = open('data/sample_augument3_standard.smi', 'w')
      
    for smile in totalsmiles:
        f.write(smile + "\n")
    f.close()
    print('Sampling completed')
    