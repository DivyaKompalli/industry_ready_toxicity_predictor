from flask_babel import Babel
from flask import request

def init_i18n(app):
    def get_locale():
        return request.accept_languages.best_match(['en', 'es'])

    Babel(app, locale_selector=get_locale)
