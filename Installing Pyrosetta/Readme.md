```
# Create an environment
conda create -n pyrosetta python=3.11 -y

# Activate it
conda activate pyrosetta

#Install pyrosetta installer
pip install pyrosetta-installer

#install PyRosetta, after installing package, run:
python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()

#Check if pyrosetta is successfully installed
python -c "from pyrosetta import *; init(); print('PyRosetta initialized successfully')"
```
