Debug Chimera
=============
	```bash
	# Other potentially necessary environment variables can be extracted from Chimera's startup script in /home/stine/.local/UCSF-Chimera64-1.10/bin/chimera
	export LD_LIBRARY_PATH=/home/stine/.local/UCSF-Chimera64-1.10/lib/:$LD_LIBRARY_PATH
	export PYTHONPATH=/home/stine/.local/UCSF-Chimera64-1.10/lib/python2.7/site-packages/:/home/stine/.local/UCSF-Chimera64-1.10/lib/python2.7/lib-dynload
	export TCL_LIBRARY=/home/stine/.local/UCSF-Chimera64-1.10/lib/tcl8.6/
	python -m pudb.run /home/stine/.local/UCSF-Chimera64-1.10/share/__main__.py --
	
	# Run a startup script (the pdb file is a model)
	python -m pudb.run /home/stine/.local/UCSF-Chimera64-1.10/share/__main__.py --nogui --debug --script startupScript.py models/1BIK.pdb models/1knt.pdb

Random Notes
============
* Chimera's C++ Atom.h is in  include/_molecule/Atom.h
* [Getting started with Chimera and Python](https://www.student.cs.uwaterloo.ca/~cs483/Getting_Started_with_Chimera_and_Python_H.pdf)
* ["getcrd" command lists an atoms coordinates](http://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/getcrd.html)

Markdown itself
===============
[Markdown's syntax](http://daringfireball.net/projects/markdown/syntax)