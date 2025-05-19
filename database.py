import sqlite3
from datetime import datetime

DB_PATH = 'predictions.db'

def init_db():
    with sqlite3.connect(DB_PATH) as conn:
        conn.execute('''
            CREATE TABLE IF NOT EXISTS predictions (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                email TEXT,
                smiles TEXT,
                prediction TEXT,
                toxicity_score REAL,
                model TEXT,
                timestamp TEXT
            )
        ''')

def log_prediction(email, smiles, prediction, score, model):
    with sqlite3.connect(DB_PATH) as conn:
        conn.execute('''
            INSERT INTO predictions (email, smiles, prediction, toxicity_score, model, timestamp)
            VALUES (?, ?, ?, ?, ?, ?)
        ''', (email, smiles, prediction, score, model, datetime.now().isoformat()))

def get_user_predictions(email):
    with sqlite3.connect(DB_PATH) as conn:
        cursor = conn.execute('''
            SELECT smiles, prediction, toxicity_score, timestamp
            FROM predictions WHERE email = ? ORDER BY timestamp DESC
        ''', (email,))
        return [
            {
                'smiles': row[0],
                'prediction': row[1],
                'toxicity_score': row[2],
                'timestamp': row[3],
                'is_toxic': row[1].lower() == 'toxic'
            } for row in cursor.fetchall()
        ]
