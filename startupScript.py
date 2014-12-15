from __init__ import startup
from chimera.tkgui import app

# If startup is immediately executed, it will miss the models that are being
# opened concurrently (Its likelyhood is significantly higher in commandline mode, which
# is used for development).
# Inspired by http://www.cgl.ucsf.edu/pipermail/chimera-users/2008-October/003226.html
# Btw. it's unfortunate that it creates a depedency towards tkgui. It means that
# this script cannot run in '--nogui' mode any longer, but life is too short to deal with
# this right now.
app.after(100, startup)
