from flask import Flask, request, jsonify, send_from_directory, render_template, redirect, url_for, session
from predictor import ToxicityPredictor, LeadOptimizer, ADMETPredictor, TargetPredictor, DrugLikeness
from auth import auth_blueprint
from database import init_db, log_prediction, get_user_predictions
from i18n import init_i18n
import os
import logging
from rdkit import Chem
from rdkit.Chem import Draw , Descriptors
from io import BytesIO
import base64
import pandas as pd


app = Flask(__name__, static_folder='static', template_folder='templates')
app.secret_key = 'your_secret_key_here'
app.register_blueprint(auth_blueprint)
init_db()
init_i18n(app)

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# Predictors
tox_predictor = ToxicityPredictor()
optimizer = LeadOptimizer()
admet_predictor = ADMETPredictor()
target_predictor = TargetPredictor()
drug_evaluator = DrugLikeness()

# Image generator

def generate_molecule_image(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol)
        buffered = BytesIO()
        img.save(buffered, format="PNG")
        return base64.b64encode(buffered.getvalue()).decode('utf-8')
    return None

@app.route('/')
def home():
    return render_template('index.html')

@app.route('/dashboard')
def dashboard():
    if 'user' not in session:
        return redirect('/login')
    history = get_user_predictions(session['user'])
    return render_template('dashboard.html', history=history)

@app.route('/predict', methods=['POST'])
def predict():
    try:
        data = request.get_json() if request.is_json else request.form
        smiles = data.get('smiles')
        model_choice = data.get('model', 'rf')

        if not smiles:
            return render_template('error.html', message="Missing SMILES string")

        toxicity_result = tox_predictor.predict_toxicity(smiles)
        druglikeness_result = drug_evaluator.evaluate(smiles)
        admet_result = admet_predictor.predict_admet(smiles)
        mol_image = generate_molecule_image(smiles)

        is_toxic = toxicity_result["General_Toxicity"]
        score = 0.75 if is_toxic else 0.25
        prediction_label = "Toxic" if is_toxic else "Non-toxic"

        if 'user' in session:
            log_prediction(session['user'], smiles, prediction_label, score, model_choice)

        result = {
            "prediction": prediction_label,
            "toxicity_score": score,
            "is_toxic": is_toxic,
            "toxicity": toxicity_result,
            "druglikeness": druglikeness_result,
            "admet": admet_result,
            "mol_image": mol_image
        }

        explanation = {
            "features": ["MolWt", "LogP"],
            "values": [0.3, -0.2]
        }

        return render_template(
            'result.html',
            result=result,
            admet=admet_result,
            druglikeness=result["druglikeness"],
            explanation=explanation
        )
    except Exception as e:
        logger.error(f"Prediction error: {str(e)}")
        return render_template('error.html', message=str(e))

@app.route('/batch_predict', methods=['POST'])
def batch_predict():
    file = request.files.get('file')
    if not file:
        return render_template("error.html", message="No file uploaded")

    df = pd.read_csv(file)
    if 'smiles' not in df.columns:
        return render_template("error.html", message="CSV must contain a 'smiles' column")

    predictions = []
    for smiles in df['smiles']:
        tox = tox_predictor.predict_toxicity(smiles)
        predictions.append({
            "smiles": smiles,
            "is_toxic": tox["General_Toxicity"],
            "prediction": "Toxic" if tox["General_Toxicity"] else "Non-toxic",
            "toxicity_score": 0.75 if tox["General_Toxicity"] else 0.25
        })

    return render_template("dashboard.html", history=predictions)

@app.route('/static/<path:filename>')
def serve_static(filename):
    return send_from_directory(app.static_folder, filename)

if __name__ == '__main__':
    for d in ['static/css', 'static/js', 'templates', 'models']:
        os.makedirs(d, exist_ok=True)
    app.run(host='0.0.0.0', port=5000, debug=True)

@app.route('/optimize', methods=['GET', 'POST'])
def optimize():
    if request.method == 'POST':
        smiles = request.form.get('smiles')
        max_tox = float(request.form.get('max_toxicity', 0.3))
        target_logp = float(request.form.get('target_logp', 2.5))

        if not smiles:
            return render_template('error.html', message="Missing SMILES for optimization")

        results = []
        try:
            mol = Chem.MolFromSmiles(smiles)
            if not mol:
                raise ValueError("Invalid SMILES input")

            # Mutate the molecule (simple approach: remove H, add methyl)
            for i in range(5):
                new_mol = Chem.RWMol(mol)
                atom = new_mol.GetAtomWithIdx(i % mol.GetNumAtoms())
                atom.SetAtomicNum((atom.GetAtomicNum() % 10) + 1)  # tweak atom type
                Chem.SanitizeMol(new_mol)
                new_smiles = Chem.MolToSmiles(new_mol)

                tox = tox_predictor.predict_toxicity(new_smiles)
                logp = Descriptors.MolLogP(new_mol)
                score = 0.75 if tox['General_Toxicity'] else 0.25
                if score <= max_tox and logp <= target_logp:
                    results.append({
                        'smiles': new_smiles,
                        'toxicity': score,
                        'logp': logp,
                        'image': generate_molecule_image(new_smiles)
                    })
        except Exception as e:
            return render_template("error.html", message=str(e))

        return render_template("optimise.html", results=results, original=smiles)
    return render_template("optimise.html", results=[], original="")

