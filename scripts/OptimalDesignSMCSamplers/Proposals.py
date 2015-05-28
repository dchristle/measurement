from __future__ import division
from scipy import stats

class OutcomeProposal:
	def __init__( self ):
		self.dist = stats.t(4.0, loc = 0.0, scale = 10.0)
	def propose(self,d):
		return d + self.dist.rvs()
	def evaluate(self, d1, d2):
		return self.dist.pdf(d2-d1)

class OutcomeLocalProposal:
	def __init__( self ):
		self.dist = stats.norm(loc = 0.0, scale = 5.0)
	def propose(self,d):
		return d + self.dist.rvs()
	def evaluate(self, d1, d2):
		return self.dist.pdf(d2-d1)

class OutcomeGlobalProposal:
	def __init__( self ):
		self.dist = stats.norm(loc = 0.0, scale = 15.0)
	def propose(self,d):
		return self.dist.rvs()
	def evaluate(self, d1, d2):
		return self.dist.pdf(d2)

class OutcomeMixtureProposal:
	"""Mixture of local and global proposal moves for experiment outcomes"""
	def __init__( self ):
		self.localProp = OutcomeLocalProposal()
		self.globalProp = OutcomeGlobalProposal()
		self.mixtureWeight = 0.10
	def propose(self,d):
		if random.random() > self.mixtureWeight:
			return self.localProp.propose(d)
		else:
			return self.globalProp.propose(d)
	def evaluate(self, d1, d2):
		return self.mixtureWeight * self.globalProp.evaluate(d1,d2) + \
			   (1.0 - self.mixtureWeight) * self.localProp.evaluate(d1,d2)

# Theses proposals are used for proposing a new design point (x location) based on the 
# previous location
class DesignLocalProposal:
	def __init__( self ):
		self.dist = stats.t(4,loc = 0, scale = 0.15)
	def propose(self,d):
		return d + self.dist.rvs()
	def evaluate(self, d1, d2):
		return self.dist.pdf(d1-d2)

class DesignGlobalProposal:
	def __init__( self ):
		self.dist = stats.uniform(loc = -2, scale = 4)
	def propose(self,d):
		return self.dist.rvs()
	def evaluate(self, d1, d2):
		return self.dist.pdf(d2)

class DesignMixtureProposal:
	def __init__( self ):
		self.localProp = DesignLocalProposal()
		self.globalProp = DesignGlobalProposal()
		self.mixtureWeight = 0.10
	def propose(self,d):
		if random.random() > self.mixtureWeight:
			return self.localProp.propose(d)
		else:
			return self.globalProp.propose(d)
	def evaluate(self, d1, d2):
		return self.mixtureWeight * self.globalProp.evaluate(d1,d2) + \
			   (1.0 - self.mixtureWeight) * self.localProp.evaluate(d1,d2)
