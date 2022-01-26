# This script will create a Python virtual environment (venv) as well as the Python packages PEtab & PEtab-select which are necessary to use arPEtabSelect i$
# Execute by navigating to this directory in the terminal and typing "bash prepare_python_venv.sh"
# Assumes a symlink named "python3.7" that calls Python 3.7

# create python venv in ~/ directory
python3.7 -m venv ~/_d2d_python_venv

# activate virtual environment
source ~/_d2d_python_venv/bin/activate

# install PEtab and PEtab-select
pip3 install --upgrade pip
pip3 install PEtab --upgrade
pip3 install PEtab-select --upgrade
