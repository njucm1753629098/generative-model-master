
if __name__ == "__main__":

    with open('data/sample_pretrain_generated.smi', 'r') as file:
        smiles_lines = file.readlines()

    unique_smiles = set()
    unique_smiles_ordered = []
    for smiles in smiles_lines:
        smiles = smiles.strip()   
        if smiles not in unique_smiles:
            unique_smiles.add(smiles)
            unique_smiles_ordered.append(smiles)   

    with open('data/sample_pretrain_generated_unique.txt', 'w') as file:
        for smiles in unique_smiles_ordered:
            file.write(smiles + '\n')