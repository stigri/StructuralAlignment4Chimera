import chimera

# Store all currently open models in the variable called "models"
models = chimera.openModels.list(modelTypes = [chimera.Molecule])

# Loop over all atoms of all open models and print the atom's coordinates
for p in models:
    print p, "\n", "#####"
    for a in p.atoms:
	print "Name:", a.name
	print "Coord:", a.coord()
	print "xform:", a.xformCoord()
	print "ID ATM Type:", a.idatmType
	print "Residue:",  a.residue.type, a.residue.id
	print "--------"

