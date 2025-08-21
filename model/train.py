#!/usr/bin/env python

import torch
from torch.utils.data import DataLoader
import pickle
from rdkit import Chem
from tqdm import tqdm
import data_struct as ds
from data_struct import MolData, Vocabulary
from data_struct import Variable, decrease_learning_rate
from model import RNN
# import sys

def pretrain(restore_from = None):
    """Trains the Prior RNN"""

    voc = Vocabulary(init_from_file="./Voc")
    moldata = MolData("data/zinc250k_filtered.txt", voc)
    data = DataLoader(moldata, batch_size=128, shuffle=True, drop_last=False,
                      collate_fn=MolData.collate_fn)
    Prior = RNN(voc)
    
    if restore_from:
        Prior.rnn.load_state_dict(torch.load(restore_from))

    optimizer = torch.optim.Adam(Prior.rnn.parameters(), lr = 0.001)
    for epoch in range(1, 51):
       
        for step, batch in tqdm(enumerate(data), total=len(data)):

            seqs = batch.long()

            log_p, _ = Prior.likelihood(seqs)
            loss = - log_p.mean()

            optimizer.zero_grad() 
            loss.backward() 
            optimizer.step() 

            # Every 500 steps we decrease learning rate and print some information
            if step % 500 == 0 and step != 0:
                #decrease_learning_rate(optimizer, decrease_by=0.003)
                decrease_learning_rate(optimizer, decrease_by=0.001)
            '''
            if step % 30 == 0 and step != 0:
                tqdm.write("*" * 50)
                tqdm.write("Epoch {:3d}   step {:3d}    loss: {:5.2f}\n".format(epoch, step, loss.item()))
                seqs, likelihood, _ = Prior.sample(100)
                valid = 0
                f = open('test_output.smi', 'a')
                for i, seq in enumerate(seqs.cpu().numpy()):
                    smile = voc.decode(seq)
                    if Chem.MolFromSmiles(smile):
                        valid += 1
                        f.write(smile + "\n")
                    if i < 10:
                        tqdm.write(smile)
                f.close()
                tqdm.write("\n{:>4.1f}% valid SMILES".format(100 * valid / len(seqs)))
                tqdm.write("*" * 50 + "\n")
            '''               
        # Save the Prior
        torch.save(Prior.rnn.state_dict(), "data/50_epoch.ckpt")

if __name__ == "__main__":
  
    smiles_file = 'data/zinc250k.txt'
    print("Reading smiles...")
    smiles_list = ds.canonicalize_smiles_from_file(smiles_file)
    print("Constructing vocabulary...")
    voc_chars = ds.construct_vocabulary(smiles_list, "./Voc")
    ds.write_smiles_to_file(smiles_list, "data/zinc250k_filtered.txt")
    pretrain()
