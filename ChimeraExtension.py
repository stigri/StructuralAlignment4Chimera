import chimera.extension
from chimera import tasks

class SA4CEMO(chimera.extension.EMO):
	# This class is derived from chimera.extension.EMO
	# Chimera uses this class to integrate the extension
	
	# Return the actual name of the extension.
	def name(self):
		return "Create/View structural Alignment"
	
	# Return the short description that typically appears as
	# balloon help or in the Tools preference category.
	def description(self):
		return "Align structures with WurstL and display their superposition"
	
	# Return the categories in which this extension should appear.
	# It is either a list or a dictionary.  If it is a dictionary
	# then the keys are the category names and the values are
	# category-specific descriptions (and the description() method
	# above is ignored).
	def categories(self):
		return ["Structural Alignment"]
	
	# Return the name of a file containing an icon that may be used
	# on the tool bar to provide a shortcut for launching the extension.
	def icon(self):
		# TODO need a tiff
		return self.path("SA4C.tiff")
	
	# Invoke the extension.  Note that when this method is called,
	# the content of "__init__.py" is not available.  To simplify
	# calling functions, the 'EMO' provides a 'module' method that
	# locates modules in the extension package by name; if no name
	# is supplied, the "__init__.py" module is returned.
	def activate(self):
		BackgroundTask(self.module("__init__").startup)

class BackgroundTask:
	def __init__(self, func):
		self.checkCount = 0
		self.task = tasks.Task("Structural Alignment", self.cancelCB, self.statusCB)
		# Execute func in a thread different from the main (eventloop) thread
		# to prevent func from freezing the UI (which is re-drawn by the eventloop).
		from thread import start_new_thread
		start_new_thread(self.execute,(func,))

	def cancelCB(self):
		# TODO cancellation doesn't work yet
		self.finished()

	def statusCB(self):
		# TODO Provide status beyond an incrementing counter
		self.checkCount += 1
		self.task.updateStatus("count=%d" % self.checkCount)

	def execute(self, func):
		# This method executes in its own thread. Once
		# func returns, we signal completion to the 
		# BackgroundTask.
		func()
		self.task.finished()

# Register an instance of 'PscViewEMO' with the Chimera extension manager
chimera.extension.manager.registerExtension(SA4CEMO(__file__))
