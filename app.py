from rdkit import DataStructs
from torch.utils.data import functional_datapipe
from werkzeug.utils import secure_filename
import tempfile
import csv
import json
from datetime import datetime
from flask import Flask, request, jsonify, send_file
import os
import pandas as pd
from collections import Counter
from flask_cors import CORS
from config import Config
import logging
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem import AllChem
import base64
import io
from Bio import Align
import sqlite3
from typing import List, Dict, Union

# 配置日志
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
app = Flask(__name__)
app.config.from_object(Config)

# 修改CORS配置，只允许特定origin
CORS(app, supports_credentials=True, resources={
    r"/*": {
        "origins": ["http://localhost:5174", "http://localhost:5173", "http://127.0.0.1:5174", "http://127.0.0.1:5173",
                    "http://192.168.0.13:5173/", ],
        "methods": ["GET", "POST", "PUT", "DELETE", "OPTIONS"],
        "allow_headers": ["Content-Type", "Authorization"],
        "expose_headers": ["Content-Range", "X-Content-Range"]
    }
})
# 配置请求限制
limiter = Limiter(
    app=app,
    key_func=get_remote_address,
    default_limits=["200 per day", "50 per hour"],
    storage_uri="memory://"
)


class DatabasePool:
    _instance = None

    def __new__(cls):
        if cls._instance is None:
            cls._instance = super().__new__(cls)
            cls._instance.connections = {}
        return cls._instance

    def get_connection(self):
        thread_id = os.getpid()
        if thread_id not in self.connections or not self._is_connection_valid(self.connections[thread_id]):
            conn = sqlite3.connect(Config.DATABASE)
            conn.row_factory = sqlite3.Row
            self.connections[thread_id] = conn
        return self.connections[thread_id]

    def _is_connection_valid(self, conn):
        try:
            conn.execute("SELECT 1")
            return True
        except (sqlite3.Error, sqlite3.ProgrammingError):
            return False

    def close_all(self):
        for conn in self.connections.values():
            conn.close()
        self.connections.clear()


db_pool = DatabasePool()


def get_db():
    return db_pool.get_connection()


@app.teardown_appcontext
def close_db(error):
    db_pool.close_all()


# 添加自定义错误处理
class APIError(Exception):
    def __init__(self, message, status_code=400):
        self.message = message
        self.status_code = status_code

    def to_dict(self):
        return {'error': self.message}


@app.errorhandler(APIError)
def handle_api_error(error):
    response = jsonify(error.to_dict())
    response.status_code = error.status_code
    return response


class StructureSearch:
    def __init__(self, database_path):
        self.database_path = database_path

    def _create_mol_from_input(self, input_type, input_data):
        """从不同输入格式创建RDKit分子对象"""
        mol = None
        try:
            if not input_data:
                raise ValueError(f"Empty {input_type} input")

            # 检测输入类型是否正确
            if input_type == 'smiles' and input_data.startswith('InChI='):
                input_type = 'inchi'
            elif input_type == 'inchi' and not input_data.startswith('InChI='):
                input_type = 'smiles'

            # 清理输入数据
            cleaned_data = input_data.strip()

            if input_type == 'smiles':
                mol = Chem.MolFromSmiles(cleaned_data)
                if mol is None:
                    # 尝试作为 InChI 解析
                    if cleaned_data.startswith('InChI='):
                        mol = Chem.MolFromInchi(cleaned_data)
                        if mol is None:
                            raise ValueError(f"Invalid InChI string: {cleaned_data}")
                    else:
                        raise ValueError(f"Invalid SMILES string: {cleaned_data}")

            elif input_type == 'inchi':
                if not cleaned_data.startswith('InChI='):
                    cleaned_data = f"InChI={cleaned_data}"
                mol = Chem.MolFromInchi(cleaned_data)
                if mol is None:
                    raise ValueError(f"Invalid InChI string: {cleaned_data}")

            elif input_type == 'pdb':
                mol = Chem.MolFromPDBBlock(cleaned_data)
                if mol is None:
                    raise ValueError(f"Invalid PDB format: {cleaned_data[:100]}...")

            if mol is None:
                raise ValueError(f"Failed to create molecule from {input_type}")

            # 标准化处理
            try:
                # 移除氢原子以标准化结构
                mol = Chem.RemoveHs(mol)

                # 添加氢原子进行3D构象生成
                mol = Chem.AddHs(mol)

                # 生成3D构象
                AllChem.EmbedMolecule(mol, randomSeed=42)

                # 优化几何构型
                AllChem.MMFFOptimizeMolecule(mol)

                # 确保分子有属性名称
                if not mol.HasProp('_Name'):
                    mol.SetProp('_Name', '')

                logger.info("Successfully created and optimized molecule")
                return mol

            except Exception as e:
                logger.warning(f"Could not optimize molecule: {str(e)}")
                # 即使优化失败也返回基本分子
                return mol

        except Exception as e:
            logger.error(f"Error processing {input_type}: {str(e)}")
            raise ValueError(f"Error processing {input_type}: {str(e)}")

    def calculate_fingerprint(self, mol):
        """计算分子指纹，使用多种指纹类型的组合"""
        try:
            if mol is None:
                raise ValueError("Invalid molecule")

            from rdkit.Chem import AllChem, MACCSkeys

            # 计算多个类型的指纹
            morgan_fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, 2048)
            maccs_fp = MACCSkeys.GenMACCSKeys(mol)
            topological_fp = AllChem.RDKFingerprint(mol)

            logger.info("Generated multiple fingerprint types")
            return {
                'morgan': morgan_fp,
                'maccs': maccs_fp,
                'topological': topological_fp
            }

        except Exception as e:
            logger.error(f"Error calculating fingerprint: {str(e)}")
            raise ValueError(f"Error calculating fingerprint: {str(e)}")

    def calculate_similarity(self, query_fps, target_fps):
        """计算综合相似度分数"""
        try:
            # 计算每种指纹类型的相似度
            morgan_sim = DataStructs.TanimotoSimilarity(query_fps['morgan'], target_fps['morgan'])
            maccs_sim = DataStructs.TanimotoSimilarity(query_fps['maccs'], target_fps['maccs'])
            topo_sim = DataStructs.TanimotoSimilarity(query_fps['topological'], target_fps['topological'])

            # 使用加权平均计算最终相似度
            weights = {'morgan': 0.4, 'maccs': 0.3, 'topological': 0.3}
            final_similarity = (
                    morgan_sim * weights['morgan'] +
                    maccs_sim * weights['maccs'] +
                    topo_sim * weights['topological']
            )

            return final_similarity

        except Exception as e:
            logger.error(f"Error calculating similarity: {str(e)}")
            return 0.0

    def search_similar_structures(self, query_mol, threshold=0.3):  # 降低默认阈值
        conn = None
        try:
            query_fps = self.calculate_fingerprint(query_mol)
            logger.info("Query fingerprints generated")

            conn = sqlite3.connect(self.database_path)
            cur = conn.cursor()

            cur.execute("""
                SELECT ID, Species, Sequence, SMILES, InChI, Molecular_Weight, 
                       Function_Description, Protein_Source
                FROM peptides 
                WHERE (SMILES IS NOT NULL AND SMILES != '') 
                   OR (InChI IS NOT NULL AND InChI != '')
            """)
            rows = cur.fetchall()
            logger.info(f"Found {len(rows)} peptides in database")

            results = []
            similarity_distribution = {
                "0.0-0.2": 0,
                "0.2-0.4": 0,
                "0.4-0.6": 0,
                "0.6-0.8": 0,
                "0.8-1.0": 0
            }

            for row in rows:
                peptide_id, species, sequence, smiles, inchi, mw, function, source = row
                try:
                    target_mol = None
                    if smiles:
                        target_mol = Chem.MolFromSmiles(smiles)
                    if target_mol is None and inchi:
                        target_mol = Chem.MolFromInchi(inchi)

                    if target_mol is None:
                        continue

                    target_fps = self.calculate_fingerprint(target_mol)
                    similarity = self.calculate_similarity(query_fps, target_fps)

                    # 记录更细粒度的相似度分布
                    if similarity < 0.2:
                        similarity_distribution["0.0-0.2"] += 1
                    elif similarity < 0.4:
                        similarity_distribution["0.2-0.4"] += 1
                    elif similarity < 0.6:
                        similarity_distribution["0.4-0.6"] += 1
                    elif similarity < 0.8:
                        similarity_distribution["0.6-0.8"] += 1
                    else:
                        similarity_distribution["0.8-1.0"] += 1

                    if similarity >= threshold:
                        results.append({
                            'id': peptide_id,
                            'species': species,
                            'sequence': sequence,
                            'smiles': smiles,
                            'inchi': inchi,
                            'molecular_weight': mw,
                            'function': function,
                            'protein_source': source,
                            'similarity': round(similarity, 3)
                        })

                except Exception as e:
                    logger.error(f"Error processing peptide {peptide_id}: {str(e)}")
                    continue

            logger.info(f"Similarity distribution: {similarity_distribution}")
            logger.info(f"Found {len(results)} matches above threshold {threshold}")

            results.sort(key=lambda x: x['similarity'], reverse=True)
            return results

        except Exception as e:
            logger.error(f"Structure search failed: {str(e)}")
            raise Exception(f"Structure search failed: {str(e)}")
        finally:
            if conn:
                conn.close()


class BlastComparison:
    def __init__(self, database_path: str):
        self.database_path = database_path
        self.aligner = Align.PairwiseAligner()
        self.aligner.mode = 'global'
        self.aligner.open_gap_score = -10
        self.aligner.extend_gap_score = -0.5
        self.aligner.match_score = 2
        self.aligner.mismatch_score = -1
    def validate_sequence(self, sequence: str) -> bool:
        """Enhanced sequence validation function"""
        valid_aa = set('ACDEFGHIKLMNPQRSTVWY')
        # Clean sequence (remove whitespace)
        cleaned_sequence = ''.join(sequence.strip().split())
        return all(aa in valid_aa for aa in cleaned_sequence.upper())

    def parse_fasta(self, fasta_content: str) -> List[Dict[str, str]]:
        """Parse FASTA format content into a list of sequences"""
        sequences = []
        current_header = None
        current_sequence = []

        for line in fasta_content.strip().split('\n'):
            if line.startswith('>'):
                if current_header:
                    sequences.append({
                        'header': current_header,
                        'sequence': ''.join(current_sequence)
                    })
                current_header = line[1:].strip()
                current_sequence = []
            else:
                current_sequence.append(line.strip())

        if current_header:  # Add the last sequence
            sequences.append({
                'header': current_header,
                'sequence': ''.join(current_sequence)
            })

        return sequences

    def validate_fasta_sequences(self, sequences: List[Dict[str, str]]) -> None:
        """Validate all sequences in a FASTA file"""
        if not sequences:
            raise ValueError("No sequences found in input")

        for seq_entry in sequences:
            sequence = seq_entry['sequence'].upper()
            if not sequence:
                raise ValueError(f"Empty sequence found for header: {seq_entry['header']}")
            if not self.validate_sequence(sequence):
                raise ValueError(
                    f"Invalid amino acid sequence for {seq_entry['header']}. "
                    "Please use standard amino acids (ACDEFGHIKLMNPQRSTVWY)"
                )

    def compare_fasta_sequences(self, fasta_content: str) -> List[Dict]:
        """Compare multiple sequences from FASTA format"""
        sequences = self.parse_fasta(fasta_content)
        self.validate_fasta_sequences(sequences)  # Validate all sequences first

        all_results = []
        for seq_entry in sequences:
            results = self.compare_sequences(seq_entry['sequence'])
            all_results.append({
                'query_header': seq_entry['header'],
                'query_sequence': seq_entry['sequence'],
                'matches': results
            })

        return all_results

    def compare_sequences(self, query_sequence: str) -> List[Dict]:
        try:
            sequence = query_sequence.strip().upper()
            conn = sqlite3.connect(self.database_path)
            conn.row_factory = sqlite3.Row  # 添加这行
            try:
                cur = conn.cursor()
                cur.execute("""
                    SELECT ID, Species, Sequence, Function_Description, Protein_Source
                    FROM peptides 
                    WHERE Sequence IS NOT NULL AND Sequence != ''
                """)

                results = []
                for row in cur.fetchall():
                    if row['Sequence']:  # 使用字典方式访问
                        score, identity = self.calculate_alignment_score(sequence, row['Sequence'])

                        if identity >= 30:
                            results.append({
                                'id': row['ID'],
                                'species': row['Species'],
                                'sequence': row['Sequence'],
                                'function': row['Function_Description'],
                                'source': row['Protein_Source'],  # 确保这个字段存在
                                'score': round(score, 2),
                                'identity': round(identity / 100, 3),
                                'sequenceLength': len(row['Sequence'])
                            })
                            # 添加日志
                            logger.info(f"BLAST result with source: {results[-1]}")

                results.sort(key=lambda x: x['score'], reverse=True)
                return results[:50]  # Return top 50 matches

            finally:
                conn.close()

        except Exception as e:
            logger.error(f"Error comparing sequences: {str(e)}")
            raise ValueError(f"Error comparing sequences: {str(e)}")

    def calculate_alignment_score(self, seq1: str, seq2: str) -> tuple:
        """Calculate alignment score and identity using PairwiseAligner."""
        alignments = self.aligner.align(seq1, seq2)
        if not alignments:
            return 0.0, 0.0

        best_alignment = alignments[0]
        score = best_alignment.score

        # Calculate identity
        aligned_seq1, aligned_seq2 = str(best_alignment[0]), str(best_alignment[1])
        matches = sum(1 for a, b in zip(aligned_seq1, aligned_seq2) if a == b)
        identity = (matches / max(len(seq1), len(seq2))) * 100

        return float(score), float(identity)

    def extract_sequence_from_fasta(self, fasta_content: str) -> str:
        """从FASTA内容中提取序列"""
        lines = fasta_content.strip().split('\n')
        if not lines:
            raise ValueError("Empty input")

        if not lines[0].startswith('>'):
            raise ValueError("Invalid FASTA format: Missing header")

        sequence_lines = [line.strip() for line in lines[1:]]
        sequence = ''.join(sequence_lines)

        if not sequence:
            raise ValueError("Empty sequence")

        return sequence

blast_comparison = BlastComparison(Config.DATABASE)


def handle_structure_file(file):
    """处理上传的结构文件"""
    if not file or file.filename == '':
        raise ValueError("No file selected")

    filename = secure_filename(file.filename)
    file_ext = os.path.splitext(filename)[1].lower()
    logger.info(f"Processing file {filename} with extension {file_ext}")

    if file_ext not in ['.mol2', '.sdf', '.pdb']:
        raise ValueError(f"Unsupported file format: {file_ext}")

    try:
        with tempfile.NamedTemporaryFile(suffix=file_ext, delete=False) as temp_file:
            file.save(temp_file.name)
            logger.info(f"File saved to temporary location: {temp_file.name}")

            # 读取文件内容进行验证
            with open(temp_file.name, 'r') as f:
                content = f.read()
                logger.info(f"File content length: {len(content)}")
                logger.debug(f"First 100 characters: {content[:100]}")

            mol = None
            if file_ext == '.mol2':
                logger.info("Attempting to parse MOL2 file...")
                mol = Chem.MolFromMol2File(temp_file.name)

                if mol is not None:
                    logger.info(f"MOL2 file parsed successfully. Atom count: {mol.GetNumAtoms()}")
                    # 添加一些分子信息的日志
                    logger.debug(f"Formula: {Chem.rdMolDescriptors.CalcMolFormula(mol)}")

            elif file_ext == '.sdf':
                logger.info("Attempting to parse SDF file...")
                supplier = Chem.SDMolSupplier(temp_file.name)
                if len(supplier) > 0:
                    mol = supplier[0]
            elif file_ext == '.pdb':
                logger.info("Attempting to parse PDB file...")
                mol = Chem.MolFromPDBFile(temp_file.name)

            if mol is None:
                logger.error(f"Failed to parse molecule from {file_ext} file")
                raise ValueError(f"Failed to parse molecule from {file_ext} file")

            return mol

    except Exception as e:
        logger.error(f"Error processing file: {str(e)}")
        raise ValueError(f"Error processing file: {str(e)}")
    finally:
        try:
            os.unlink(temp_file.name)
            logger.info("Temporary file deleted")
        except:
            logger.warning("Failed to delete temporary file")


def handle_file_upload(file, allowed_extensions):
    if not file or file.filename == '':
        raise APIError("No file selected")

    if not allowed_file(file.filename):
        raise APIError(f"File types not allowed. Supported types: {allowed_extensions}")

    try:
        filename = secure_filename(file.filename)
        temp_file = tempfile.NamedTemporaryFile(delete=False)
        file.save(temp_file.name)
        return temp_file.name
    except Exception as e:
        logger.error(f"File upload error: {str(e)}")
        raise APIError(f"Failed to save file: {str(e)}")


def allowed_file(filename):
    return '.' in filename and filename.rsplit('.', 1)[1].lower() in Config.ALLOWED_EXTENSIONS


def validate_peptide_data(data):
    """验证提交的肽数据"""
    required_fields = ['sequence', 'link', 'submittedBy']

    # 检查必需字段
    for field in required_fields:
        if field not in data:
            raise ValueError(f"Missing required field: {field}")
        if not data[field]:
            raise ValueError(f"Field cannot be empty: {field}")

    # 验证序列格式
    sequence = data['sequence'].upper()
    valid_amino_acids = set('ACDEFGHIKLMNPQRSTVWY')
    if not all(aa in valid_amino_acids for aa in sequence):
        raise ValueError("Invalid amino acid sequence")

    return True


@app.route('/')
def home():
    response = {
        'status': 'running',
        'message': 'HYPLPDB API Server',
        'version': '1.0',
        'timestamp': datetime.now().isoformat()
    }
    return jsonify(response)


@app.route('/api/peptides/<id>', methods=['GET'])
def get_peptide(id):
    try:
        logger.info(f"Attempting to fetch peptide with ID: {id}")
        conn = get_db()
        cur = conn.cursor()

        if id.isdigit():
            id = f"HYPLPDB{int(id):04d}"
        elif id.startswith('HYPLPDB'):
            id = f"HYPLPDB{id[7:]}"  # 替换前缀
        logger.info(f"Final query ID: {id}")

        # 使用实际的列名查询所有需要的字段
        query = '''
            SELECT 
                ID, Species, Sequence, Protein_Source, Method___Extraction,
                Function1, Function2, Function3, Function4, Function5,
                Function_Description, IC50, Validation,
                Length, Molecular_Weight, Isoelectric_Point,
                Charge, Gravy, Instability_Index, Aromaticity,
                SMILES, InChI,
                PMIDs_DOI1, PMIDs_DOI1_link, PMIDs_DOI2, PMIDs_DOI2_link,
                UNIPROT, UNIPROT_link,
                BIOPEPUWM, BIOPEPUWM_link,
                FermFooDb, FermFooDb_link,
                BIOPEPDB, BIOPEPDB_link,
                PlantPepDB, PlantPepDB_link,
                CyclicPepedia, CyclicPepedia_link,
                JPO, JPO_link
            FROM peptides 
            WHERE ID = ?
        '''
        cur.execute(query, (id,))
        row = cur.fetchone()

        if row is None:
            logger.info(f"No peptide found with ID: {id}")
            return jsonify({'error': 'Peptide not found'}), 404

        cur.execute(query, (id,))
        row = cur.fetchone()

        if row is None:
            return jsonify({'error': 'Peptide not found'}), 404

        structure_data = None
        if row['SMILES'] or row['InChI']:
            mol = None
            if row['SMILES']:
                mol = Chem.MolFromSmiles(row['SMILES'])
            elif row['InChI']:
                mol = Chem.MolFromInchi(row['InChI'])

            if mol:
                # 生成2D结构SVG
                img = Draw.MolToImage(mol)
                # 将图像转换为base64
                import io
                import base64
                img_io = io.BytesIO()
                img.save(img_io, 'PNG')
                img_io.seek(0)
                img_base64 = base64.b64encode(img_io.read()).decode()

                # 生成3D结构
                mol_3d = Chem.AddHs(mol)
                AllChem.EmbedMolecule(mol_3d, randomSeed=42)
                AllChem.MMFFOptimizeMolecule(mol_3d)
                pdb = Chem.MolToPDBBlock(mol_3d)

                structure_data = {
                    '2d_image': f"data:image/png;base64,{img_base64}",
                    '3d_data': pdb
                }

        # 构建返回数据
        peptide = {
            'basicInfo': {
                'id': row['ID'],
                'species': row['Species'],
                'sequence': row['Sequence'],
                'length': row['Length'],
                'source': row['Protein_Source'],
                'extractionMethod': row['Method___Extraction'],
                'activity': {
                    'functions': [
                        row[f'Function{i}'] for i in range(1, 6)
                        if row[f'Function{i}'] is not None and row[f'Function{i}'].strip() != ''
                    ],
                    'functionDescription': row['Function_Description'],
                    'ic50': row['IC50'],
                    'validationModel': row['Validation']
                }
            },
            'physicalProperties': {
                'molecularWeight': row['Molecular_Weight'],
                'isoelectricPoint': row['Isoelectric_Point'],
                'chargeAtPH7': row['Charge'],
                'hydrophobicity': row['Gravy'],
                'instabilityIndex': row['Instability_Index'],
                'aromaticity': row['Aromaticity']
            },
            'structure': {
                'smiles': row['SMILES'],
                'inchi': row['InChI'],
                'visualization': structure_data
            },
            'externalReferences': {
                'literature': [
                    {
                        'id': row['PMIDs_DOI1'],
                        'link': f"https://pubmed.ncbi.nlm.nih.gov/{row['PMIDs_DOI1_link']}/" if row[
                            'PMIDs_DOI1'] else None,
                        'display': row['PMIDs_DOI1']
                    },
                    {
                        'id': row['PMIDs_DOI2'],
                        'link': f"https://pubmed.ncbi.nlm.nih.gov/{row['PMIDs_DOI2_link']}/" if row[
                            'PMIDs_DOI2'] else None,
                        'display': row['PMIDs_DOI2']
                    }
                ],
                'databases': {
                    'uniprot': {
                        'id': row['UNIPROT'],
                        'link': f"https://www.uniprot.org/uniprotkb/{row['UNIPROT']}/entry" if row['UNIPROT'] else None,
                        'display': row['UNIPROT']
                    },
                    'biopepuwm': {
                        'id': row['BIOPEPUWM'],
                        'link': f"https://biochemia.uwm.edu.pl/biopep/peptide_data_page1.php?zm_ID={row['BIOPEPUWM']}" if
                        row['BIOPEPUWM'] else None,
                        'display': row['BIOPEPUWM']
                    },
                    'fermfoodb': {
                        'id': row['FermFooDb'],
                        'link': f"https://webs.iiitd.edu.in/raghava/fermfoodb/display_sub.php?details={row['FermFooDb']}" if
                        row['FermFooDb'] else None,
                        'display': row['FermFooDb']
                    },
                    'biopepdb': {
                        'id': row['BIOPEPDB'],
                        'link': f"https://bis.zju.edu.cn/biopepdbr/index.php?p=detail&ID={row['BIOPEPDB']}" if row[
                            'BIOPEPDB'] else None,
                        'display': row['BIOPEPDB']
                    },
                    'plantpepdb': {
                        'id': row['PlantPepDB'],
                        'link': f"http://14.139.61.8/PlantPepDB/pages/information.php?id={row['PlantPepDB']}" if row[
                            'PlantPepDB'] else None,
                        'display': row['PlantPepDB']
                    },
                    'cyclicpepedia': {
                        'id': row['CyclicPepedia'],
                        'link': f"https://www.biosino.org/iMAC/cyclicpepedia/detail?id={row['CyclicPepedia']}" if row[
                            'CyclicPepedia'] else None,
                        'display': row['CyclicPepedia']
                    },
                    'jpo': {
                        'id': row['JPO'],
                        'link': f"https://www.ebi.ac.uk/Tools/dbfetch/dbfetch?db=patent_prt;id={row['JPO'],}" if row[
                            'JPO'] else None,
                        'display': row['JPO']
                    }
                }
            }
        }

        return jsonify(peptide)

    except Exception as e:
        logger.error(f"Error fetching peptide: {str(e)}", exc_info=True)
        return jsonify({'error': 'Failed to fetch peptide'}), 500

    # 添加辅助函数来格式化链接
    def _format_doi_link(link):
        if not link:
            return None
        if link.startswith('http'):
            return link
        return f"https://doi.org/{link}"

    def _format_external_link(link):
        if not link:
            return None
        if link.startswith('http'):
            return link
        if link.startswith('www.'):
            return f"https://{link}"
        return link



@app.route('/api/peptides/search', methods=['GET'])
def search_peptides():
    try:
        query = request.args.get('q', '')
        page = int(request.args.get('page', 1))
        per_page = int(request.args.get('per_page', 10))
        exact_match = request.args.get('exact_match', 'false').lower() == 'true'

        id_query = request.args.get('id', '')
        species_query = request.args.get('species', '')
        sequence_query = request.args.get('sequence', '')
        functionDescription_query = request.args.get('functionDescription', '')

        mw_min = request.args.get('mw_min', type=float)
        mw_max = request.args.get('mw_max', type=float)
        pi_min = request.args.get('pi_min', type=float)
        pi_max = request.args.get('pi_max', type=float)
        charge_min = request.args.get('charge_min', type=float)
        charge_max = request.args.get('charge_max', type=float)

        ic50 = request.args.get('ic50')
        pmid = request.args.get('pmid')
        uniprot = request.args.get('uniprot')
        biopep = request.args.get('biopep')

        conn = get_db()
        cur = conn.cursor()

        query_parts = []
        params = []

        if query:
            if exact_match:
                query_parts.append("(ID = ? OR Species = ? OR Sequence = ? OR Function_Description = ?)")
                normalized_id = query
                if query.isdigit():
                    normalized_id = f"HYPLPDB{int(query):04d}"
                elif query.upper().startswith('HYPLPDB'):
                    normalized_id = f"HYPLPDB{query[7:]}"
                # 修复此处参数数量问题
                params.extend([normalized_id, query, query, query])  # 修复重复扩展
            else:
                query_parts.append("(ID LIKE ? OR Species LIKE ? OR Sequence LIKE ? OR Function_Description LIKE ?)")
                search_param = f'%{query}%'
                params.extend([search_param] * 4)  # 修复多出参数

        if id_query:
            normalized_id = id_query
            if id_query.isdigit():
                normalized_id = f"HYPLPDB{int(id_query):04d}"
            elif id_query.upper().startswith('HYPLPDB'):
                normalized_id = f"HYPLPDB{id_query[7:]}"

            if exact_match:
                query_parts.append("ID = ?")
                params.append(normalized_id)
            else:
                query_parts.append("ID LIKE ?")
                params.append(f"%{normalized_id}%")

        if species_query:
            query_parts.append("Species = ?" if exact_match else "Species LIKE ?")
            params.append(species_query if exact_match else f'%{species_query}%')

        if sequence_query:
            query_parts.append("Sequence = ?" if exact_match else "Sequence LIKE ?")
            params.append(sequence_query if exact_match else f'%{sequence_query}%')

        if functionDescription_query:
            query_parts.append("Function_Description = ?" if exact_match else "Function_Description LIKE ?")
            params.append(functionDescription_query if exact_match else f'%{functionDescription_query}%')

        if mw_min is not None:
            query_parts.append("Molecular_Weight >= ?")
            params.append(mw_min)
        if mw_max is not None:
            query_parts.append("Molecular_Weight <= ?")
            params.append(mw_max)

        if pi_min is not None:
            query_parts.append("Isoelectric_Point >= ?")
            params.append(pi_min)
        if pi_max is not None:
            query_parts.append("Isoelectric_Point <= ?")
            params.append(pi_max)

        if charge_min is not None:
            query_parts.append("Charge >= ?")
            params.append(charge_min)
        if charge_max is not None:
            query_parts.append("Charge <= ?")
            params.append(charge_max)

        if ic50:
            query_parts.append("IC50 LIKE ?")
            params.append(f'%{ic50}%')

        if pmid:
            query_parts.append("(PMIDs_DOI1 = ? OR PMIDs_DOI2 = ?)")
            params.extend([pmid, pmid])  # 确保只扩展两次

        if uniprot:
            query_parts.append("UNIPROT = ?")
            params.append(uniprot)

        if biopep:
            query_parts.append("(BIOPEPUWM = ? OR BIOPEPDB = ?)")
            params.extend([biopep, biopep])

        base_query = '''
            SELECT 
                ID, Species, Sequence, Protein_Source,
                Length, Molecular_Weight, 
                Isoelectric_Point, Charge,
                IC50, Function_Description
            FROM peptides
        '''

        where_clause = " WHERE " + " AND ".join(f"({part})" for part in query_parts) if query_parts else ""
        order_clause = " ORDER BY SPECIES"

        # `LIMIT` 和 `OFFSET` 参数修正
        limit_clause = " LIMIT ? OFFSET ?"
        params.extend([per_page, (page - 1) * per_page])  # 确保最后的参数严格对应 LIMIT 和 OFFSET

        full_query = base_query + where_clause + order_clause + limit_clause

        logger.debug(f"Search query: {full_query}")
        logger.debug(f"Search parameters: {params}")

        cur.execute(full_query, params)

        items = []
        for row in cur.fetchall():
            items.append({
                'id': row['ID'],
                'species': row['Species'],
                'sequence': row['Sequence'],
                'source': row['Protein_Source'],
                'sequenceLength': row['Length'],
                'molecularWeight': row['Molecular_Weight'],
                'isoelectricPoint': row['Isoelectric_Point'],
                'chargeAtPH7': row['Charge'],
                'ic50': row['IC50'],
                'functionDescription': row['Function_Description']
            })

        count_query = 'SELECT COUNT(*) FROM peptides' + where_clause
        cur.execute(count_query, params[:-2])  # 确保移除 LIMIT 和 OFFSET 参数
        total = cur.fetchone()[0]

        return jsonify({
            'items': items,
            'pagination': {
                'total': total,
                'page': page,
                'per_page': per_page,
                'total_pages': (total + per_page - 1) // per_page
            },
            'search_info': {
                'query': query,
                'exact_match': exact_match,
                'filters_applied': len(query_parts)
            }
        })

    except Exception as e:
        logger.error(f"Search error: {str(e)}", exc_info=True)
        return jsonify({'error': 'Search failed', 'message': str(e)}), 500
    finally:
        if conn:
            conn.close()



def get_paginated_results(query, params, page, per_page):
    """优化的分页查询"""
    offset = (page - 1) * per_page

    # 使用 WITH 子句优化分页
    paginated_query = f"""
    WITH RankedResults AS (
        SELECT *, ROW_NUMBER() OVER (ORDER BY SPECIES) as row_num
        FROM ({query})
    )
    SELECT * FROM RankedResults 
    WHERE row_num > ? AND row_num <= ?
    """

    params.extend([offset, offset + per_page])
    return paginated_query


@app.route('/api/peptides/structure-search', methods=['POST'])
@limiter.limit("20 per minute")
def structure_search():
    try:
        # Get pagination parameters
        page = int(request.args.get('page', 1))
        per_page = int(request.args.get('per_page', 10))

        # Initialize structure searcher
        structure_searcher = StructureSearch(Config.DATABASE)
        threshold = float(request.form.get('threshold', 0.3))

        # Handle either file upload or SMILES/InChI input
        mol = None
        if 'file' in request.files:
            file = request.files['file']
            if file and file.filename != '':
                mol = handle_structure_file(file)
        elif 'smiles' in request.form:
            smiles = request.form['smiles']
            if smiles:
                mol = structure_searcher._create_mol_from_input('smiles', smiles)
        elif 'inchi' in request.form:
            inchi = request.form['inchi']
            if inchi:
                mol = structure_searcher._create_mol_from_input('inchi', inchi)

        if mol is None:
            return jsonify({
                'error': 'No valid molecular structure provided. Please provide either a file or SMILES/InChI input.'
            }), 400

        logger.info("Starting similarity search...")
        all_results = structure_searcher.search_similar_structures(mol, threshold=threshold)
        logger.info(f"Search completed. Found {len(all_results)} matches")

        # Calculate pagination
        total = len(all_results)
        total_pages = (total + per_page - 1) // per_page
        start_idx = (page - 1) * per_page
        end_idx = start_idx + per_page

        # Paginate results
        paginated_results = all_results[start_idx:end_idx]

        return jsonify({
            'results': paginated_results,
            'pagination': {
                'total': total,
                'page': page,
                'per_page': per_page,
                'total_pages': total_pages
            }
        })

    except ValueError as ve:
        logger.error(f"Value error in structure search: {str(ve)}")
        return jsonify({'error': str(ve)}), 400
    except Exception as e:
        logger.error(f"Structure search error: {str(e)}")
        return jsonify({'error': 'Search failed', 'message': str(e)}), 500


@app.route('/api/peptides/blast-compare', methods=['POST'])
def blast_comparison_endpoint():
    try:
        app.logger.info('Received BLAST request')

        if 'file' in request.files:
            file = request.files['file']
            if file.filename == '':
                return jsonify({'error': 'No file selected'}), 400

            try:
                file_content = file.read().decode('utf-8')
                if not file_content.strip():
                    return jsonify({'error': 'Empty file'}), 400

                if not file_content.strip().startswith('>'):
                    return jsonify({'error': 'Invalid FASTA format. File must start with ">"'}), 400

                results = blast_comparison.compare_fasta_sequences(file_content)

            except UnicodeDecodeError:
                return jsonify({'error': 'Invalid file format. Please upload a text file.'}), 400

        elif 'sequence' in request.form:
            sequence = request.form['sequence'].strip()
            if not sequence:
                return jsonify({'error': 'Empty sequence'}), 400

            # Format single sequence as FASTA if needed
            if not sequence.startswith('>'):
                sequence = f'>Query\n{sequence}'

            results = blast_comparison.compare_fasta_sequences(sequence)

        else:
            return jsonify({'error': 'No sequence or file provided'}), 400

        # Process results
        all_matches = []
        for result in results:
            all_matches.extend(result['matches'])

        # Sort all matches by score
        all_matches.sort(key=lambda x: x['score'], reverse=True)

        response_data = {
            'results': all_matches[:50],  # Return top 50 matches overall
            'pagination': {
                'total': len(all_matches),
                'page': 1,
                'per_page': 50,
                'total_pages': 1
            }
        }

        return jsonify(response_data)

    except ValueError as ve:
        app.logger.error(f'Validation error: {str(ve)}')
        return jsonify({'error': str(ve)}), 400
    except Exception as e:
        app.logger.error(f'Unexpected error: {str(e)}', exc_info=True)
        return jsonify({'error': 'Server error', 'message': str(e)}), 500

@app.route('/api/sequence-analysis', methods=['GET'])
def analyze_sequences():
    # 读取数据
    try:
        file_path = os.getenv('ANALYSIS_FILE_PATH', 'backend/peptides.csv')
        df = pd.read_csv(file_path, encoding='latin1')

        # 1. 计算序列长度分布
        df['Sequence_Length'] = df['Sequence'].apply(len)
        sequence_length_freq = df['Sequence_Length'].value_counts().sort_index()

        # 2. 计算氨基酸组成
        amino_acid_frequencies = df['Sequence'].apply(lambda seq: Counter(seq))
        total_amino_acids = Counter()
        for amino_acid_count in amino_acid_frequencies:
            total_amino_acids.update(amino_acid_count)

        # 3. 统计功能描述
        functionDescription_freq = df['Function_Description'].value_counts()

        # 4. 统计物种
        species_freq = df['Species'].value_counts()

        return jsonify({
            'sequence_length': {
                'labels': sequence_length_freq.index.tolist(),
                'values': sequence_length_freq.values.tolist()
            },
            'amino_acids': {
                'labels': list(total_amino_acids.keys()),
                'values': list(total_amino_acids.values())
            },
            'functionDescription': {
                'labels': functionDescription_freq.index.tolist(),
                'values': functionDescription_freq.values.tolist()
            },
            'species': {
                'labels': species_freq.index.tolist(),
                'values': species_freq.values.tolist()
            }
        })
    except Exception as e:
        logger.error(f"Analysis error: {str(e)}")
        raise APIError("Analysis failed", 500)


@app.route('/api/peptides/download', methods=['GET'])
def download_peptides():
    filepath = None
    try:
        format_type = request.args.get('format', 'csv')  # 默认为CSV格式
        conn = get_db()
        cur = conn.cursor()

        # 查询所有肽数据
        query = '''
            SELECT 
                ID, Species, Sequence, Protein_Source, Method___Extraction,
                Function1, Function2, Function3, Function4, Function5,
                Function_Description, IC50, Validation,
                Length, Molecular_Weight, Isoelectric_Point,
                Charge, Gravy, Instability_Index, Aromaticity,
                SMILES, InChI
            FROM peptides
        '''
        cur.execute(query)
        rows = cur.fetchall()
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')

        if format_type == 'json':
            # 准备JSON格式
            data = []
            for row in rows:  # 使用已获取的rows
                peptide = dict(row)  # 转换为字典
                data.append(peptide)

            filename = f'peptides_{timestamp}.json'
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)

            with open(filepath, 'w', encoding='utf-8') as f:
                json.dump(data, f, ensure_ascii=False, indent=2)

        else:  # CSV格式
            filename = f'peptides_{timestamp}.csv'
            filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)

            with open(filepath, 'w', newline='', encoding='utf-8') as csvfile:
                writer = csv.writer(csvfile)
                # 写入表头
                writer.writerow([description[0] for description in cur.description])
                # 写入数据
                writer.writerows(rows)

        return send_file(
            filepath,
            as_attachment=True,
            download_name=filename,
            mimetype='text/csv' if format_type == 'csv' else 'application/json'
        )

    except Exception as e:
        logger.error(f"Download error: {str(e)}")
        return jsonify({'error': 'Download failed', 'message': str(e)}), 500
    finally:
        if conn:
            conn.close()
        # 清理临时文件
        if filepath and os.path.exists(filepath):
            try:
                os.remove(filepath)
            except Exception as e:
                logger.error(f"Error removing temp file: {str(e)}")


@app.route('/api/peptides/submit', methods=['POST'])
def submit_peptide():
    conn = None
    try:
        # 获取JSON数据
        data = request.get_json()
        if not data:
            return jsonify({'error': 'No data provided'}), 400

        # 验证数据
        try:
            validate_peptide_data(data)
        except ValueError as e:
            return jsonify({'error': 'Validation failed', 'message': str(e)}), 400

        conn = get_db()
        cur = conn.cursor()

        # 检查是否存在相同序列
        cur.execute('SELECT ID FROM peptides WHERE Sequence = ?', (data['sequence'],))
        if cur.fetchone():
            return jsonify({
                'error': 'Duplicate sequence',
                'message': 'A peptide with this sequence already exists'
            }), 409

        # 修改插入语句，只包含必需字段和可选字段
        insert_query = '''
            INSERT INTO Submission (
                Sequence, Link, Submitted_By,
                Name, Source, Method___Extraction,
                Function_Description, IC50, Validation,
                Submitted_Date, Status
            ) VALUES (
                ?, ?, ?,
                ?, ?, ?,
                ?, ?, ?,
                ?, ?
            )
        '''

        params = (
            data['sequence'],  # 必需字段
            data['link'],  # 必需字段
            data['submittedBy'],  # 必需字段
            data.get('name', ''),  # 可选字段
            data.get('source', ''),  # 可选字段
            data.get('extractionMethod', ''),  # 可选字段
            data.get('functionDescription', ''),  # 可选字段
            data.get('ic50', ''),  # 可选字段
            data.get('validation', ''),  # 可选字段
            datetime.now().isoformat(),  # 提交时间
            'pending'  # 状态
        )

        cur.execute(insert_query, params)
        conn.commit()

        return jsonify({
            'message': 'Peptide submitted successfully',
            'id': cur.lastrowid
        })

    except sqlite3.Error as e:
        logger.error(f"Database error: {str(e)}")
        return jsonify({'error': 'Database error', 'message': str(e)}), 500
    except Exception as e:
        logger.error(f"Submission error: {str(e)}")
        return jsonify({'error': 'Submission failed', 'message': str(e)}), 500
    finally:
        if conn:
            conn.close()


if __name__ == '__main__':
    app.run(debug=True)
