
class AcceptanceStats:
	"""Utility class to keep track of acceptance ratios"""
	def __init__(self):
		self.accepts = 0
		self.rejects = 0

	def addAccept(self):
		self.accepts += 1

	def addReject(self):
		self.rejects += 1
	
	def getAcceptanceRatio(self):
		return float(self.accepts) / float(self.accepts + self.rejects)

	def reset(self):
		self.accepts = 0
		self.rejects = 0
