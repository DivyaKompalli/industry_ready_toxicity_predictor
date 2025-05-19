# Authentication logic with role-based access
from flask import Blueprint, render_template, request, redirect, session

auth_blueprint = Blueprint('auth', __name__)

# Simulated user store (replace with real DB in production)
USERS = {
    "admin@example.com": {"password": "adminpass", "role": "admin"},
    "user@example.com": {"password": "userpass", "role": "user"}
}

@auth_blueprint.route('/login', methods=['GET', 'POST'])
def login():
    if request.method == 'POST':
        email = request.form['email']
        password = request.form['password']
        user = USERS.get(email)
        if user and user['password'] == password:
            session['user'] = email
            return redirect('/')
        else:
            return render_template('login.html', error="Invalid credentials")
    return render_template('login.html')

@auth_blueprint.route('/signup', methods=['GET', 'POST'])
def signup():
    if request.method == 'POST':
        email = request.form['email']
        password = request.form['password']
        if email in USERS:
            return render_template('signup.html', error="User already exists")
        USERS[email] = {"password": password, "role": "user"}
        return redirect('/login')
    return render_template('signup.html')

@auth_blueprint.route('/logout')
def logout():
    session.clear()
    return redirect('/login')

def login_required(f):
    from functools import wraps
    @wraps(f)
    def wrapper(*args, **kwargs):
        if 'user' not in session:
            return redirect('/login')
        return f(*args, **kwargs)
    return wrapper

def get_user_role(email):
    return USERS.get(email, {}).get("role", "guest")
