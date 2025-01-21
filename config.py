import os
from dotenv import load_dotenv

load_dotenv()

class Config:
    DATABASE = os.getenv('DATABASE_PATH', 'D:/HYPLPDB/data/peptides.db')
    UPLOAD_FOLDER = os.getenv('UPLOAD_FOLDER', 'D:/HYPLPDB/temp')
    MAX_CONTENT_LENGTH = 16 * 1024 * 1024  # 16MB
    ALLOWED_EXTENSIONS = {'csv', 'xlsx', 'mol2', 'sdf', 'pdb'}
    ALLOWED_ORIGIN = 'http://localhost:5173','http://localhost:5174'
