import joblib
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED, rdMolDescriptors
import os

class ToxicityPredictor:
    def __init__(self):
        try:
            model_dir = os.path.join(os.path.dirname(__file__), 'models')
            self.model_ar = joblib.load(os.path.join(model_dir, 'tox21_NR-AR_rf.joblib'))
            self.model_general = joblib.load(os.path.join(model_dir, 'tox21_rf_model.joblib'))
        except Exception as e:
            raise RuntimeError(f"Failed to load models: {str(e)}")

    def predict_toxicity(self, smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return {"error": "Invalid SMILES string"}

            features = [[
                Descriptors.MolWt(mol),
                Descriptors.NumHAcceptors(mol),
                Descriptors.NumHDonors(mol),
                Descriptors.TPSA(mol),
                Descriptors.MolLogP(mol),
                rdMolDescriptors.CalcNumRotatableBonds(mol),
                rdMolDescriptors.CalcNumRings(mol)
            ]]

            return {
                "NR-AR": bool(self.model_ar.predict(features)[0]),
                "General_Toxicity": bool(self.model_general.predict(features)[0]),
                "error": None
            }
        except Exception as e:
            return {"error": str(e)}

class LeadOptimizer:
    def __init__(self):
        self.tox_predictor = ToxicityPredictor()
    
    def optimize(self, smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return {"error": "Invalid SMILES string"}
            
            return {
                "toxicity": self.tox_predictor.predict_toxicity(smiles),
                "properties": {
                    "Molecular_Weight": Descriptors.MolWt(mol),
                    "LogP": Descriptors.MolLogP(mol),
                    "H_Donors": Lipinski.NumHDonors(mol),
                    "H_Acceptors": Lipinski.NumHAcceptors(mol),
                    "TPSA": Descriptors.TPSA(mol)
                },
                "error": None
            }
        except Exception as e:
            return {"error": str(e)}

class ADMETPredictor:
    def predict_admet(self, smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return {"error": "Invalid SMILES string"}
            
            return {
                "Absorption": 0.8,
                "Distribution": 0.6,
                "Metabolism": 0.7,
                "Excretion": 0.5,
                "Toxicity": 0.3,
                "error": None
            }
        except Exception as e:
            return {"error": str(e)}

class TargetPredictor:
    def predict_targets(self, smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return {"error": "Invalid SMILES string"}
            
            return {
                "targets": ["Protein_A", "Protein_B"],
                "scores": [0.9, 0.7],
                "error": None
            }
        except Exception as e:
            return {"error": str(e)}

class DrugLikeness:
    def evaluate(self, smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                return {"error": "Invalid SMILES string"}
            
            lipinski = {
                "Molecular_Weight": Descriptors.MolWt(mol) <= 500,
                "LogP": Descriptors.MolLogP(mol) <= 5,
                "H_Donors": Lipinski.NumHDonors(mol) <= 5,
                "H_Acceptors": Lipinski.NumHAcceptors(mol) <= 10
            }
            
            return {
                "QED": round(QED.qed(mol), 2),
                "Lipinski": lipinski,
                "Passes_Lipinski": all(lipinski.values()),
                "error": None
            }
        except Exception as e:
            return {"error": str(e)}
