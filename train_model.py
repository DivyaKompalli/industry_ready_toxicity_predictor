# train_model.py

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import roc_auc_score, accuracy_score
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
import joblib
import os
from pathlib import Path
import argparse

# Configuration
MODEL_NAME = 'rf'
MODEL_DIR = Path('models')
DATA_DIR = Path('data/tox21-challenge')

def clean_tox21_scores(score_series):
    """Convert Tox21 scores (0,1,x) to numeric (0,1,NA)"""
    return pd.to_numeric(
        score_series.replace('x', np.nan),
        errors='coerce'
    )

def load_tox21_data(assay_name='NR-AR'):
    """Load and preprocess Tox21 data"""
    try:
        smiles_path = DATA_DIR / 'tox21_10k_challenge_score.smiles'
        scores_path = DATA_DIR / 'tox21_10k_challenge_score.txt'

        smiles_df = pd.read_csv(
            smiles_path,
            sep='\t',
            header=None,
            names=['smiles', 'sample_id'],
            skiprows=1
        )

        scores_df = pd.read_csv(
            scores_path,
            sep='\t',
            header=0
        )

        # Assay mapping
        assay_columns = {
            'NR-AR': 'NR-AR', 'SR-ARE': 'SR-ARE', 'NR-AR-LBD': 'NR-AR-LBD',
            'AAR': 'NR-AhR', 'AR': 'NR-AR', 'AR-LBD': 'NR-AR-LBD',
            'ARE': 'SR-ARE', 'aromatase': 'NR-Aromatase',
            'ATAD5': 'SR-ATAD5', 'ER': 'NR-ER',
            'ER-LBD': 'NR-ER-LBD', 'HSE': 'SR-HSE',
            'MMP': 'SR-MMP', 'p53': 'SR-p53',
            'PPAR-gamma': 'NR-PPAR-gamma'
        }

        if assay_name not in assay_columns:
            raise ValueError(f"Assay {assay_name} not found")

        # Combine and clean
        df = smiles_df.copy()
        score_column = assay_columns[assay_name]
        df['activity'] = clean_tox21_scores(scores_df[score_column])
        df = df[df['activity'].notna()]
        df['toxicity'] = (df['activity'] == 1).astype(int)

        # Validate SMILES
        df['valid'] = df['smiles'].apply(lambda x: Chem.MolFromSmiles(str(x)) is not None)
        df = df[df['valid']]

        if len(df) == 0:
            raise ValueError("No valid compounds remaining")

        print(f"Loaded {len(df)} valid compounds for {assay_name}")
        return df

    except Exception as e:
        print(f"Error loading data: {str(e)}")
        raise

def train_toxicity_model(assay_name='NR-AR'):
    """Train and save model"""
    try:
        MODEL_DIR.mkdir(parents=True, exist_ok=True)

        print(f"\nTraining model for {assay_name}")
        df = load_tox21_data(assay_name)

        # Featurization
        X = []
        y = []
        for smiles, tox in zip(df['smiles'], df['toxicity']):
            mol = Chem.MolFromSmiles(str(smiles))
            X.append([
                Descriptors.MolWt(mol),
                Descriptors.NumHAcceptors(mol),
                Descriptors.NumHDonors(mol),
                Descriptors.TPSA(mol),
                Descriptors.MolLogP(mol),
                rdMolDescriptors.CalcNumRotatableBonds(mol),
                rdMolDescriptors.CalcNumRings(mol)
            ])
            y.append(tox)

        # Split
        X_train, X_test, y_train, y_test = train_test_split(
            np.array(X), np.array(y), test_size=0.2, random_state=42
        )

        # Train
        model = RandomForestClassifier(
            n_estimators=100,
            max_depth=10,
            random_state=42,
            class_weight='balanced'
        )
        model.fit(X_train, y_train)

        # Evaluate
        test_pred = model.predict_proba(X_test)[:, 1]
        print(f"\nTest Performance:")
        print(f"AUC: {roc_auc_score(y_test, test_pred):.3f}")
        print(f"Accuracy: {accuracy_score(y_test, (test_pred > 0.5).astype(int)):.3f}")

        # Save model
        model_path = MODEL_DIR / f'tox21_{assay_name}_{MODEL_NAME}.joblib'
        joblib.dump(model, str(model_path))
        print(f"Model saved to {model_path}")

    except Exception as e:
        print(f"Training failed: {str(e)}")
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Train Tox21 prediction model')
    parser.add_argument('--assay', default='NR-AR', help='Tox21 assay name (default: NR-AR)')
    args = parser.parse_args()

    print("Starting model training...")
    train_toxicity_model(args.assay)
