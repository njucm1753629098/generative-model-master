import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, QED, rdMolDescriptors, RDConfig
import os
import sys
sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

def calculate_properties(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None, None, None, None, None, None, None, None
        logp = Descriptors.MolLogP(mol)
        mw = Descriptors.MolWt(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rotbonds = Descriptors.NumRotatableBonds(mol)
        qed = QED.qed(mol)
        sas = sascorer.calculateScore(mol)
        return logp, mw, tpsa, hbd, hba, rotbonds, qed, sas
if __name__ == '__main__':
    file_path = 'data/sample_augument3_standard_unique_novel.txt'
    #file_path = 'data/inhibitor2_filtered.txt'
    #file_path = 'data/zinc250k_filtered.txt'
    with open(file_path, 'r') as file:
        smiles_list = file.readlines()
    smiles_list = [smiles.strip() for smiles in smiles_list]
    property_data = []
    for smiles in smiles_list:
        logp, mw, tpsa, hbd, hba, rotbonds, qed, sas = calculate_properties(smiles)
        property_data.append([smiles, logp, mw, tpsa, hbd, hba, rotbonds, qed, sas])

    property_df = pd.DataFrame(property_data, columns=['smiles', 'logp', 'mw', 'tpsa', 'hbd', 'hba', 'rotbonds', 'qed', 'sas'])
    property_df.to_csv('data/sample_augument3_standard_unique_novel_with_property.csv', index=False)
    #property_df.to_csv('data/inhibitor2_property.csv', index=False)
    #property_df.to_csv('data/zinc250k_filtered_property.csv', index=False)