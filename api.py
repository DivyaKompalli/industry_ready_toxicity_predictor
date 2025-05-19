from flask import Blueprint, request, jsonify
from predictor import ToxicityPredictor

tox_predictor = ToxicityPredictor()
api_blueprint = Blueprint('api', __name__)

API_KEYS = {
    "test-key-123": "user@example.com",
    "admin-key-456": "admin@example.com"
}

@api_blueprint.route('/api/predict', methods=['POST'])
def api_predict():
    key = request.headers.get('X-API-Key')
    if key not in API_KEYS:
        return jsonify({"error": "Unauthorized - Invalid API key"}), 401

    data = request.get_json()
    smiles = data.get("smiles")
    if not smiles:
        return jsonify({"error": "Missing SMILES string"}), 400

    result = tox_predictor.predict_toxicity(smiles)
    return jsonify({"user": API_KEYS[key], "smiles": smiles, "toxicity": result})

