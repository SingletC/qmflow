"""Routes for parent Flask app."""
import os
import pathlib

import flask
from flask import current_app as app
from flask import render_template
from rdkit import Chem
from rdkit.Chem import AllChem


@app.route("/")
def home():
    """Landing page."""
    return render_template(
        "index.jinja2",
        title="Dash",
        description="",
        template="home-template",
        body="This is a homepage served with Flask.",
    )
@app.route('/downloads/<path:path>')
def serve_static(path):
    if path[-5:] != '.fchk':
        return ''
    root_dir = os.getcwd()
    return flask.send_from_directory(
        os.path.join(root_dir), path
    )

