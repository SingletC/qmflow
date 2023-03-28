"""Routes for parent Flask app."""
import os
import pathlib

import flask
from flask import current_app as app
from flask import render_template
from rdkit import Chem
from rdkit.Chem import AllChem

from calc.utils import rdkit_2_base64png, rdkit2png


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

@app.route('/smiles/<string:smiles>')
def serve_smiles(smiles):
    rdkit_mol = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(rdkit_mol)
    img = rdkit2png(rdkit_mol)
    return flask.send_file(img,mimetype='image/png')