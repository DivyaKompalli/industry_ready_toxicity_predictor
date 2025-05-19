from flask import Blueprint, request, send_file, render_template
import io
import pandas as pd
from database import get_user_predictions
from fpdf import FPDF

export_blueprint = Blueprint('export', __name__)

@export_blueprint.route('/export/excel')
def export_excel():
    from flask import session
    if 'user' not in session:
        return "Unauthorized", 401

    data = get_user_predictions(session['user'])
    df = pd.DataFrame(data)
    buf = io.BytesIO()
    df.to_excel(buf, index=False)
    buf.seek(0)
    return send_file(buf, download_name="tox_results.xlsx", as_attachment=True)

@export_blueprint.route('/export/pdf')
def export_pdf():
    from flask import session
    if 'user' not in session:
        return "Unauthorized", 401

    data = get_user_predictions(session['user'])
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial", size=12)
    pdf.cell(200, 10, txt="Toxicity Predictions", ln=True, align="C")

    for row in data:
        line = f"{row['timestamp']} - {row['smiles']} - {row['prediction']} ({row['toxicity_score']})"
        pdf.cell(200, 10, txt=line, ln=True)

    buf = io.BytesIO()
    pdf.output(buf)
    buf.seek(0)
    return send_file(buf, download_name="tox_results.pdf", as_attachment=True)

