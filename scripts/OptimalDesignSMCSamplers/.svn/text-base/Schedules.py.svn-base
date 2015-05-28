from __future__ import division

class RealValuedLinearSchedule:
	"""starts at a given start value, increases output value by constant stepsize at each iteration"""
	def __init__(self, startvalue, stepsize):
		self.startvalue = startvalue
		self.stepsize = stepsize

	def __call__(self,	iteration ):
		return self.startvalue + iteration * self.stepsize
		
class Interpolate0to1Schedule:
	"""starts at 0, ends up at 1 after the given number of steps"""
	def __init__(self, numberSteps):
		self.numberSteps = numberSteps

	def __call__(self,	iteration ):
		return iteration / (self.numberSteps - 1)


class IntegerStepsSchedule:
	def __init__(self, subIntegerSteps):
		self.subIntegerSteps = subIntegerSteps
	
	def __call__(self, iteration):
		"""Returns n(t) and nu(t) for a given iteration t
		   n(t) is the number of outcomes
		   nu(t) is the annealing factor for the last outcome"""
		n = iteration // self.subIntegerSteps + 1
		return (n,1)


class LinearFractionalSchedule:
	"""This schedule is designed to exactly hit integer values. It starts with n = 1, nu = 1.0 for iteration 0 """
	def __init__(self, subIntegerSteps):
		"""The parameter subIntegerSteps specifies the number of steps """
		self.subIntegerSteps = subIntegerSteps
		self.startIteration = 2 * self.subIntegerSteps - 1
	
	def __call__(self, iteration):
		"""Returns n(t) and nu(t) for a given iteration t
		   n(t) is the number of outcomes
		   nu(t) is the annealing factor for the last outcome (slowly fading in to 1.0)"""
		n = (iteration + self.startIteration) // self.subIntegerSteps
		nu = (((iteration + self.startIteration) % self.subIntegerSteps) + 1) / self.subIntegerSteps
		return (n,nu)

class Fractional0to1Schedule (object):
	"""This schedule is designed to exactly hit integer values. It starts with n = 1, nu = 1.0 for iteration 0.
	   The difference to LinearFractionalSchedule is that after n, nu=1.0 it goes to n+1, nu=0.0 """
	def __init__(self, subIntegerSteps):
		"""The parameter subIntegerSteps specifies the number of steps """
		self.subIntegerSteps = subIntegerSteps
		self.startIteration =  self.subIntegerSteps 
	
	def __call__(self, iteration):
		"""Returns n(t) and nu(t) for a given iteration t
		   n(t) is the number of outcomes
		   nu(t) is the annealing factor for the last outcome (slowly fading in to 1.0)"""
		n = (iteration + self.startIteration) // self.subIntegerSteps
		nu = ((iteration + self.startIteration) % self.subIntegerSteps)  / (self.subIntegerSteps - 1)
		return (n,nu)

if __name__ == '__main__':
	sched = IntegerStepsSchedule(3)
	for i in xrange(10):
		print "sched(",i,") = ", sched(i)
	
