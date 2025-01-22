from rdkit import Chem
from rdkit.Chem import AllChem
from typing import Dict, Optional
import json
import threading
import sqlite3
import os
import time


class PeptideStructure:
    def __init__(self, database_path: str):
        self._database_path = database_path
        self._cache_file = 'peptide_structures_cache.json'
        self._structure_cache = {}
        self._cache_lock = threading.Lock()
        self._aa_smiles = {
            'A': 'CC(N)C(=O)O',  # Alanine
            'C': 'NC(CS)C(=O)O',  # Cysteine
            'D': 'NC(CC(=O)O)C(=O)O',  # Aspartic acid
            'E': 'NC(CCC(=O)O)C(=O)O',  # Glutamic acid
            'F': 'NC(Cc1ccccc1)C(=O)O',  # Phenylalanine
            'G': 'NCC(=O)O',  # Glycine
            'H': 'NC(Cc1cnc[nH]1)C(=O)O',  # Histidine
            'I': 'NC(C(CC)C)C(=O)O',  # Isoleucine
            'K': 'NC(CCCCN)C(=O)O',  # Lysine
            'L': 'NC(CC(C)C)C(=O)O',  # Leucine
            'M': 'NC(CCSC)C(=O)O',  # Methionine
            'N': 'NC(CC(=O)N)C(=O)O',  # Asparagine
            'P': 'O=C(O)C1CCCN1',  # Proline
            'Q': 'NC(CCC(=O)N)C(=O)O',  # Glutamine
            'R': 'NC(CCCNC(=N)N)C(=O)O',  # Arginine
            'S': 'NC(CO)C(=O)O',  # Serine
            'T': 'NC(C(O)C)C(=O)O',  # Threonine
            'V': 'NC(C(C)C)C(=O)O',  # Valine
            'W': 'NC(Cc1c[nH]c2ccccc12)C(=O)O',  # Tryptophan
            'Y': 'NC(Cc1ccc(O)cc1)C(=O)O'  # Tyrosine
        }
        self._valid_sequences = self._load_valid_sequences()
        self._load_or_generate_structures()

    def _load_valid_sequences(self) -> set:
        """从数据库加载所有有效的肽序列"""
        valid_sequences = set()
        try:
            conn = sqlite3.connect(self._database_path)
            cur = conn.cursor()
            cur.execute("SELECT DISTINCT Sequence FROM peptides")
            valid_sequences = {row[0] for row in cur.fetchall()}
        except Exception as e:
            print(f"Error loading sequences: {e}")
        finally:
            if 'conn' in locals():
                conn.close()
        return valid_sequences

    def _load_cache_from_file(self) -> Dict:
        """从文件加载缓存的结构数据"""
        if os.path.exists(self._cache_file):
            try:
                with open(self._cache_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                print(f"Error loading cache file: {e}")
        return {}

    def _save_cache_to_file(self):
        """将结构数据保存到缓存文件"""
        try:
            with open(self._cache_file, 'w') as f:
                json.dump(self._structure_cache, f)
        except Exception as e:
            print(f"Error saving cache file: {e}")

    def _create_peptide_mol(self, sequence: str) -> Optional[Chem.Mol]:
        """根据序列创建肽分子"""
        try:
            if not sequence or not all(aa in self._aa_smiles for aa in sequence):
                return None

            mol = Chem.MolFromSmiles(self._aa_smiles[sequence[0]])
            if len(sequence) == 1:
                return mol

            for aa in sequence[1:]:
                next_aa = Chem.MolFromSmiles(self._aa_smiles[aa])
                mol = Chem.CombineMols(mol, next_aa)

            return mol
        except Exception as e:
            print(f"Error in create_peptide_mol: {e}")
            return None

    def _generate_structure(self, sequence: str) -> Dict:
        """为单个序列生成结构数据"""
        try:
            mol = self._create_peptide_mol(sequence)
            if mol is None:
                return {"error": "Invalid sequence"}

            mol = Chem.AddHs(mol)
            AllChem.EmbedMolecule(mol, randomSeed=42)
            AllChem.MMFFOptimizeMolecule(mol)

            conf = mol.GetConformer()
            structure_data = {
                "atoms": [],
                "bonds": [],
                "secondary_structure": self._predict_secondary_structure(mol)
            }

            for i in range(mol.GetNumAtoms()):
                atom = mol.GetAtomWithIdx(i)
                pos = conf.GetAtomPosition(i)
                structure_data["atoms"].append({
                    "id": i,
                    "element": atom.GetSymbol(),
                    "coordinates": {
                        "x": float(pos.x),
                        "y": float(pos.y),
                        "z": float(pos.z)
                    }
                })

            for bond in mol.GetBonds():
                structure_data["bonds"].append({
                    "start": bond.GetBeginAtomIdx(),
                    "end": bond.GetEndAtomIdx(),
                    "order": int(bond.GetBondTypeAsDouble())
                })

            return structure_data

        except Exception as e:
            print(f"Error generating structure: {e}")
            return {"error": str(e)}

    def _predict_secondary_structure(self, mol) -> Dict:
        """预测二级结构倾向"""
        try:
            conf = mol.GetConformer()
            num_atoms = mol.GetNumAtoms()

            dihedrals = []
            for i in range(num_atoms - 3):
                angle = AllChem.GetDihedralDeg(conf, i, i + 1, i + 2, i + 3)
                dihedrals.append(angle)

            helix_tendency = sum(1 for angle in dihedrals if -60 <= angle <= -45)
            sheet_tendency = sum(1 for angle in dihedrals if -140 <= angle <= -110)

            total = len(dihedrals) if dihedrals else 1
            return {
                "alpha_helix_probability": helix_tendency / total,
                "beta_sheet_probability": sheet_tendency / total,
                "disordered_probability": 1 - (helix_tendency + sheet_tendency) / total
            }

        except Exception as e:
            print(f"Error in secondary structure prediction: {e}")
            return {
                "alpha_helix_probability": 0,
                "beta_sheet_probability": 0,
                "disordered_probability": 1
            }

    def _load_or_generate_structures(self):
        """加载或生成结构数据"""
        print("Loading structures from cache...")
        self._structure_cache = self._load_cache_from_file()

        # 检查是否需要生成新的结构
        missing_sequences = self._valid_sequences - set(self._structure_cache.keys())
        if missing_sequences:
            print(f"Generating structures for {len(missing_sequences)} new sequences...")
            count = 0
            total = len(missing_sequences)
            start_time = time.time()

            for seq in missing_sequences:
                self._structure_cache[seq] = self._generate_structure(seq)
                count += 1
                if count % 10 == 0:
                    elapsed = time.time() - start_time
                    eta = (elapsed / count) * (total - count)
                    print(f"Progress: {count}/{total} ({(count / total) * 100:.1f}%) - ETA: {eta:.0f}s")

            print("Saving new structures to cache...")
            self._save_cache_to_file()
        else:
            print("All structures loaded from cache.")

    def get_structure(self, sequence: str) -> Optional[Dict]:
        """获取指定序列的结构（仅当序列在数据库中存在时）"""
        if sequence not in self._valid_sequences:
            return None
        return self._structure_cache.get(sequence)


# 测试代码
if __name__ == "__main__":
    DATABASE_PATH = 'D:/HYPLPDB/data/peptides.db'  # 替换为实际的数据库路径
    generator = PeptideStructure(DATABASE_PATH)

    # 测试获取结构
    test_seq = "ACDEG"  # 替换为数据库中实际存在的序列
    result = generator.get_structure(test_seq)
    if result:
        print(json.dumps(result, indent=2))
    else:
        print("Sequence not found in database")
