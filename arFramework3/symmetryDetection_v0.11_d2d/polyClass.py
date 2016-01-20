# Author: Benjamin Merkt, Physikalisches Institut, Universitaet Freiburg

import sympy as spy
import numpy as np

from copy import deepcopy

### efficient class for polynomial calculations
class Apoly:
	def __init__(self, expr, variables, rs):
		if expr is None:
			self.coefs = []
			self.exps = []
			self.vars = variables
			self.rs = rs
		else:
			poly = spy.Poly(expr, variables).as_dict()
			#extract coefficients from polynomial
			if rs is None:
				self.coefs = poly.values()
			else:
				coefsTmp = poly.values()
				self.coefs = [0]*len(coefsTmp)
				for i in xrange(len(coefsTmp)):
					self.coefs[i] = np.zeros(len(rs))
					for j, r in enumerate(rs):
						if coefsTmp[i].has(r):
							self.coefs[i][j] = spy.diff(coefsTmp[i], r)
			#extract exponents from polynomial
			self.exps = poly.keys()
			for i in xrange(len(self.exps)):
				self.exps[i] = np.array(self.exps[i])
			self.vars = variables
			self.rs = rs

	def __repr__(self):
		return str(self.coefs) + '\n' + str(self.exps)
	def __str__(self):
		return str(self.coefs) + '\n' + str(self.exps)

	### return a copy of self
	def getCopy(self):
		newPoly  = Apoly(None,self.vars, self.rs)
		newPoly.coefs = deepcopy(self.coefs)
		newPoly.exps = deepcopy(self.exps)
		
		return newPoly

	### add a second polynomial
	### self is overwritten with result
	def add(self, otherPoly):
		for i in xrange(len(otherPoly.exps)):
			for j in xrange(len(self.exps)):
				if np.array_equal(otherPoly.exps[i], self.exps[j]):
					self.coefs[j] = self.coefs[j] + otherPoly.coefs[i]
					if not np.any(self.coefs[j]):
						self.coefs.pop(j)
						self.exps.pop(j)
					break
			else:
				self.coefs.append(otherPoly.coefs[i])
				self.exps.append(otherPoly.exps[i])

	### substract a second polynomial
	### self is overwritten with result
	def sub(self, otherPoly):
		for i in xrange(len(otherPoly.exps)):
			for j in xrange(len(self.exps)):
				if np.array_equal(otherPoly.exps[i], self.exps[j]):
					self.coefs[j] = self.coefs[j] - otherPoly.coefs[i]
					if not np.any(self.coefs[j]):
						self.coefs.pop(j)
						self.exps.pop(j)
					break
			else:
				self.coefs.append(-1*otherPoly.coefs[i])
				self.exps.append(otherPoly.exps[i])

	### multiply with a second polynomial
	### a new Apoly is created and returned. self remains unchanged
	def mul(self, otherPoly):
		newPoly = Apoly(None, self.vars, self.rs)
		newPoly.coefs = [0]*(len(self.coefs)*len(otherPoly.coefs))
		newPoly.exps = [0]*(len(self.coefs)*len(otherPoly.coefs))

		k = 0
		for i in xrange(len(otherPoly.exps)):
			for j in xrange(len(self.exps)):
				newPoly.coefs[k] = otherPoly.coefs[i] * self.coefs[j] #works only because only one poly has rs
				newPoly.exps[k] = otherPoly.exps[i] + self.exps[j]
				k += 1

		i = 0
		while i < len(newPoly.coefs):
			j = i+1
			while j <len(newPoly.coefs):
				if np.array_equal(newPoly.exps[i], newPoly.exps[j]):
					newPoly.exps.pop(j)
					newPoly.coefs[i] = newPoly.coefs[i] + newPoly.coefs.pop(j)
				else:
					j += 1
			i += 1

		return newPoly

	### differentiate the polynomial
	### a new Apoly is created and returned. self remains unchanged
	def diff(self, j):
		newPoly = self.getCopy()
		i = 0
		while i < len(newPoly.exps):
			if newPoly.exps[i][j] != 0:
				newPoly.coefs[i] = newPoly.coefs[i]*newPoly.exps[i][j]
				newPoly.exps[i][j] -= 1
				i += 1
			else:
				newPoly.coefs.pop(i)
				newPoly.exps.pop(i)
		return newPoly

	### transform polynomial to regular sympy expression
	def as_expr(self):
		expr = 0
		for i in range(len(self.coefs)):
			fact = 1
			for j in range(len(self.vars)):
				fact = fact*self.vars[j]**self.exps[i][j]
			if self.rs is None:
				expr += self.coefs[i]*fact
			else:
				coef = 0
				for j in range(len(self.rs)):
					coef += self.rs[j]*self.coefs[i][j]
				expr += coef*fact
		return spy.nsimplify(expr)
