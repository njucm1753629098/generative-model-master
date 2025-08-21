from rdkit import Chem
from rdkit.Chem import AllChem


def count_chiral(smiles_list):
    """
    Count how many SMILES strings in the input list contain at least one chiral center.

    Parameters
    ----------
    smiles_list : list[str]
        List of SMILES strings.

    Returns
    -------
    tuple[int, int, int]
        total_molecules : total number of molecules processed
        chiral_count    : number of molecules with at least one chiral center
        failed_count    : number of molecules that could not be parsed
    """
    total_molecules = len(smiles_list)
    chiral_count = 0
    failed_count = 0

    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            failed_count += 1
            continue

        mol = Chem.AddHs(mol)
        try:
            AllChem.EmbedMolecule(mol)
            chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
            if chiral_centers:
                chiral_count += 1
        except Exception:
            continue

    return total_molecules, chiral_count, failed_count


if __name__ == "__main__":
    with open('data/generate_329.txt') as f:
        lines = [l.strip() for l in f if l.strip()]

    total, chiral_cnt, failed = count_chiral(lines)
    valid = total - failed

    print(f'Total: {total} | Failed: {failed}')
    print(f'Contains chiral center: {chiral_cnt} ({chiral_cnt / valid * 100:.2f}% of valid)')