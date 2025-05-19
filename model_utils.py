import joblib
import os

def load_model(model_name='rf'):
    path = f'models/tox21_{model_name}_model.joblib'
    if not os.path.exists(path):
        raise FileNotFoundError(f"Model file not found: {path}")
    print(f"Loaded model from {path}")
    return joblib.load(path)

def detect_data_drift(reference_data, new_data):
    import numpy as np
    drift_score = np.linalg.norm(np.mean(reference_data, axis=0) - np.mean(new_data, axis=0))
    return drift_score > 1.0
