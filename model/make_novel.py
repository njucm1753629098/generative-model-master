if __name__ == "__main__":
    '''
    with open('data/tl_filtered.smi', 'r') as file:
        rl_filtered_smiles = [line.strip() for line in file]
    '''
    with open('data/inhibitor2_filtered.txt', 'r') as file:
        rl_filtered_smiles = [line.strip() for line in file]
    
    
    
    with open('data/sample_pretrain_generated_unique.txt', 'r') as file:
        tran_output_smiles = [line.strip() for line in file]
    
    
    new_output_smiles = []
    
   
    for smiles in tran_output_smiles:
        if smiles not in rl_filtered_smiles:
            new_output_smiles.append(smiles)
    
    
    with open('data/sample_pretrain_generated_unique_novel.txt', 'w') as file:
        for smiles in new_output_smiles:
            file.write(smiles + '\n')