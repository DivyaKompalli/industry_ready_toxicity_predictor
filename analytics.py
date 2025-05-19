import sqlite3
from flask import Blueprint, render_template

analytics_blueprint = Blueprint('analytics', __name__)

@analytics_blueprint.route('/admin/analytics')
def analytics():
    conn = sqlite3.connect('predictions.db')
    cursor = conn.cursor()
    cursor.execute("SELECT email, COUNT(*) FROM predictions GROUP BY email")
    usage_by_user = cursor.fetchall()

    cursor.execute("SELECT model, COUNT(*) FROM predictions GROUP BY model")
    model_usage = cursor.fetchall()

    cursor.execute("SELECT COUNT(*) FROM predictions")
    total = cursor.fetchone()[0]

    conn.close()

    return render_template("admin_dashboard.html", total=total, user_stats=usage_by_user, model_stats=model_usage)

