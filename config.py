"""Flask config."""
from os import environ, path

from dotenv import load_dotenv

BASE_DIR = path.abspath(path.dirname(__file__))
load_dotenv(path.join(BASE_DIR, ".env"))


class Config:
    """Flask configuration variables."""

    # General Config
    FLASK_APP = "wsgi.py"
    FLASK_ENV = environ.get("FLASK_ENV")
    SECRET_KEY = environ.get("SECRET_KEY")
    BASIC_AUTH_FORCE = True
    BASIC_AUTH_USERNAME = environ.get("BASIC_AUTH_USERNAME")
    BASIC_AUTH_PASSWORD = environ.get("BASIC_AUTH_PASSWORD")
    # Assets
    LESS_BIN = environ.get("LESS_BIN")
    ASSETS_DEBUG = environ.get("ASSETS_DEBUG")
    LESS_RUN_IN_DEBUG = environ.get("LESS_RUN_IN_DEBUG")

    # Static Assets
    STATIC_FOLDER = "static"
    TEMPLATES_FOLDER = "templates"
    COMPRESSOR_DEBUG = environ.get("COMPRESSOR_DEBUG")

    # Gaussian
    GAUSSIAN_CMD = environ.get("GAUSSIAN_CMD")
    GAUSSIAN_N = environ.get("GAUSSIAN_N")
    GAUSSIAN_M = environ.get("GAUSSIAN_M")