import chimera
import wurstl
from chimera import replyobj

def startup():
    # Store all currently open Chimera models in the variable called "models"
    models = chimera.openModels.list()

    # Need at least to models to work with (technically we need exactly two).
    if len(models) < 2:
	replyobj.error("Cannot create alignment out of less than two models.")
	return

    # Loop over all atoms of all open models and convert Chimera's
    # data model into a user-type (model, atom, residue) independent
    # representation solely based on Python primitives. This can be
    # passed to wurstl and further converted to C data types.
    #
    # The resulting model data structure looks like:
    #
    # ["compnd first model",
    #         [("0", "idamtype", ..., x, y, z),
    #          (second atom), 
    #          ...,
    #          (n-th atom)
    #         ],
    #  "compnd second model",
    #         [(first atom of second model),
    #          (...),
    #          (n-th atom of second model)
    #         ]
    # ]
    #
    model = []
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

    # wurstl returns a list of pairs. Each pair has the indices
    # of an atom in the first model aligned to an atom in the
    # second model. alignments is therefore a mapping between
    # atoms in both models. Note that we expect to only receive
    # alignments for two models even though that the loop above
    # might iterate over more models. However, wurstl as of
    # today cannot align more than two models.
    alignments = wurstl.structural_alignment(model)

    # mobileAtoms and referenceAtoms are logically mapped lists
    # of aligned atoms. Lets use the pairs returned by wurstl
    # to lookup the pairs in the first and second model.
    mobileAtoms = []
    referenceAtoms = []
    for alignment in alignments:
		idxA = models[0].atoms[alignment[0]]
		idxB = models[1].atoms[alignment[1]]
		if idxA != -1 and idxB != -1:
			mobileAtoms.append(idxA)
			referenceAtoms.append(idxB)

    # Use Midas' match command to visually align
    # both models.
    # see https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/match.html
    from Midas import match, TooFewAtomsError
    try:
	match(mobileAtoms, referenceAtoms, None)
    except TooFewAtomsError:
	replyobj.error('Too few corresponding atoms (<4) to '
		'match both models\n')
