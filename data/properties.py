from Bio.SeqUtils.ProtParam import ProteinAnalysis
import pandas as pd
import sqlite3
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem
import re


def clean_sequence(sequence):
    """清理序列,保留标准氨基酸并记录修改"""
    standard_aas = set('ACDEFGHIKLMNPQRSTVWY')
    original = sequence

    # 移除特殊字符和空格
    sequence = re.sub(r'[^A-Z]', '', sequence.upper())

    # 检查是否包含非标准氨基酸
    non_standard = set(sequence) - standard_aas
    if non_standard:
        print(f"序列 {original} 包含非标准氨基酸 {non_standard}, 这些将被忽略")
        sequence = ''.join(aa for aa in sequence if aa in standard_aas)

    return sequence if sequence else None


def generate_structure(sequence, generator):
    """生成分子结构并返回SMILES和InChI"""
    try:
        # 使用rdkit生成分子
        mol = generator._create_peptide_mol(sequence)
        if mol is None:
            return None, None

        # 生成SMILES和InChI
        smiles = Chem.MolToSmiles(mol)
        inchi = Chem.MolToInchi(mol)

        return smiles, inchi
    except Exception as e:
        print(f"生成结构时出错: {str(e)}")
        return None, None


def calculate_peptide_properties(sequence, structure_generator=None):
    """计算肽序列的所有理化性质"""
    try:
        # 清理序列
        cleaned_sequence = clean_sequence(sequence)
        if not cleaned_sequence:
            print(f"序列 {sequence} 清理后没有有效的氨基酸")
            return None

        # 创建ProteinAnalysis对象
        peptide = ProteinAnalysis(cleaned_sequence)

        # 基本性质
        properties = {
            'original_sequence': sequence,
            'cleaned_sequence': cleaned_sequence,
            'sequence_length': len(cleaned_sequence),
            'molecular_weight': round(peptide.molecular_weight(), 2),
            'isoelectric_point': round(peptide.isoelectric_point(), 2),
            'charge_at_pH7': round(peptide.charge_at_pH(7.0), 2),
            'gravy': round(peptide.gravy(), 3),
            'instability_index': round(peptide.instability_index(), 2),
            'aromaticity': round(peptide.aromaticity(), 3)
        }

        # 如果提供了结构生成器,生成SMILES和InChI
        if structure_generator:
            smiles, inchi = generate_structure(cleaned_sequence, structure_generator)
            properties['smiles'] = smiles
            properties['inchi'] = inchi

        return properties

    except Exception as e:
        print(f"处理序列时出错: {sequence}")
        print(f"错误信息: {str(e)}")
        return None


def update_database(db_path, table_name, sequence_column, data):
    """将计算的性质更新到数据库中"""
    if not data:
        print("没有数据需要更新")
        return

    conn = sqlite3.connect(db_path)
    cursor = conn.cursor()

    try:
        # 获取现有的列名
        cursor.execute(f"PRAGMA table_info({table_name})")
        existing_columns = [row[1] for row in cursor.fetchall()]

        # 为新的性质添加列
        sample_data = next(iter(data.values()))
        new_columns = []
        for prop in sample_data.keys():
            if prop not in existing_columns:
                col_type = 'TEXT' if isinstance(sample_data[prop], str) else 'REAL'
                try:
                    cursor.execute(f"ALTER TABLE {table_name} ADD COLUMN {prop} {col_type}")
                    new_columns.append(prop)
                except sqlite3.OperationalError as e:
                    print(f"添加列 {prop} 时出错: {str(e)}")
                    continue

        if new_columns:
            print(f"成功添加了以下新列: {', '.join(new_columns)}")

        # 更新数据
        success_count = 0
        error_count = 0
        for seq, props in data.items():
            try:
                # 构建更新语句
                update_pairs = []
                values = []
                for k, v in props.items():
                    if k in existing_columns or k in new_columns:
                        update_pairs.append(f"{k}=?")
                        values.append(v)

                if update_pairs:
                    update_query = f"""
                        UPDATE {table_name} 
                        SET {', '.join(update_pairs)}
                        WHERE {sequence_column}=?
                    """
                    values.append(seq)
                    cursor.execute(update_query, values)
                    success_count += 1
            except Exception as e:
                print(f"更新序列 {seq} 数据时出错: {str(e)}")
                error_count += 1
                continue

        conn.commit()
        print(f"成功更新了 {success_count} 条记录")
        if error_count > 0:
            print(f"更新失败 {error_count} 条记录")

    except Exception as e:
        print(f"更新数据库时出错: {str(e)}")
        conn.rollback()

    finally:
        conn.close()


def main():
    # 配置参数
    db_path = r"D:/HYPLPDB/data/peptides.db"
    table_name = "peptides"
    sequence_column = "Sequence"

    try:
        # 创建结构生成器
        from structure_search import PeptideStructure
        structure_generator = PeptideStructure(db_path)

        # 从数据库读取序列
        conn = sqlite3.connect(db_path)
        sequences = pd.read_sql_query(f"SELECT {sequence_column} FROM {table_name}", conn)
        conn.close()

        print(f"成功读取到 {len(sequences)} 条序列记录")

        # 计算所有序列的性质
        results = {}
        total = len(sequences)
        for idx, row in sequences.iterrows():
            sequence = row[sequence_column]
            print(f"处理序列 {idx + 1}/{total}")

            props = calculate_peptide_properties(sequence, structure_generator)
            if props:
                # 使用原始序列作为键
                results[sequence] = props

        print(f"成功处理了 {len(results)} 条序列")

        # 更新数据库
        update_database(db_path, table_name, sequence_column, results)
        print("处理完成！")

    except Exception as e:
        print(f"运行出错: {str(e)}")


if __name__ == "__main__":
    main()