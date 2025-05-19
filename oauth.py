from flask import Blueprint, redirect, url_for, session
from flask_dance.contrib.google import make_google_blueprint, google
import os

os.environ['OAUTHLIB_INSECURE_TRANSPORT'] = '1'

oauth_blueprint = Blueprint('oauth', __name__)
google_bp = make_google_blueprint(
    client_id="GOOGLE_CLIENT_ID",
    client_secret="GOOGLE_CLIENT_SECRET",
    scope=["profile", "email"]
)

oauth_blueprint.register_blueprint(google_bp, url_prefix="/login")

@oauth_blueprint.route("/google")
def google_login():
    if not google.authorized:
        return redirect(url_for("google.login"))
    resp = google.get("/oauth2/v2/userinfo")
    email = resp.json()["email"]
    session['user'] = email
    return redirect("/")
