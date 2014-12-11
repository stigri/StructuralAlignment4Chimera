import chimera
import wurstl
from chimera.tkgui import app


def startup():
    # Store all currently open models in the variable called "models"
    models = chimera.openModels.list()

    model = []

    # Loop over all atoms of all open models and print the atom's coordinates
    for p in models:
	model.append(p.name)

	# convert all atoms of the molecule into a list of atoms
        atoms = []
        for a in p.atoms:
		atom = []
		atom.append(a.name) # str
		atom.append(a.idatmType) # str
		atom.append(a.residue.type) # str
		atom.append(a.residue.id.insertionCode) # str
		atom.append(a.residue.id.chainId) # str
		atom.append(a.residue.id.position) # int
		atom.append(a.coord().x) # float
		atom.append(a.coord().y) # float
		atom.append(a.coord().z) # float

		# stuff list of atom properties into list of atoms
		atoms.append(atom)

	model.append(atoms)

    print wurstl.structural_alignment(model)

# If startup is immediately executed, it will miss the models that are being
# opened (Its likelyhood is significantly higher in commandline mode, which
# is used for development).
# Inspired by http://www.cgl.ucsf.edu/pipermail/chimera-users/2008-October/003226.html
# Btw. it's unfortunate that it creates a depedency towards tkgui. It means that
# this script cannot run in '--nogui' mode any longer.
app.after(100, startup)
