# Author: Benjamin Merkt, Physikalisches Institut, Universitaet Freiburg

import sys

import sympy as spy
import numpy as np

from multiprocessing import Queue, Queue, Process

# try/except necessary for R interface (imports automatically and does not find other files)
try:
	from functions import *
	from polyClass import *
except:
	pass

### calculate conditions for a differential equation
def doEquation(k, numerators, denominators, derivativesNum, infis,
		diffInfis, allVariables, rs, ansatz, queue):

	n = len(allVariables)
	m = len(numerators)
	
	polynomial = Apoly(None, allVariables, rs)
	if ansatz == 'uni' or ansatz == 'par':
		#calculate polynomial
		polynomial.add(diffInfis[0][k].mul(denominators[k]).mul(numerators[k]))
		for i in range(n):
			polynomial.sub(infis[i].mul(derivativesNum[k][i]))

	elif ansatz == 'multi':
		for j in range(m):
			summand = diffInfis[k][j].mul(denominators[k]).mul(numerators[j])
			for l in range(m):
				if l != j:
					summand = summand.mul(denominators[l])
			polynomial.add(summand)
		for i in range(n):
			summand = infis[i].mul(derivativesNum[k][i])
			for l in range(m):
				if l != k:
					summand = summand.mul(denominators[l])
			polynomial.sub(summand)
	#determine rSystem such that the coefficients vanish
	lgs = np.empty([len(polynomial.coefs), len(rs)])
	for i in range(len(polynomial.coefs)):
		lgs[i,:] = polynomial.coefs[i]

	queue.put(lgs)

### calculate conditions for an observation equation
def doObsEquation(k, obsDerivativesNum, infis, allVariables, rs, queue):
	n = len(allVariables)

	#calculate polynomial
	polynomial = Apoly(None, allVariables, rs)
	for l in range(n):
		polynomial.add(infis[l].mul(obsDerivativesNum[k][l]))

	#determine rSystem such that the coefficients vanish
	lgs = np.empty([len(polynomial.coefs), len(rs)])
	for i in range(len(polynomial.coefs)):
		lgs[i,:] = polynomial.coefs[i]

	queue.put(lgs)

### calculate conditions for an initial equation
def doInitEquation(k, initDenominators, initDerivativesNum,
			initFunctions, infis, allVariables, rs, queue):
	n = len(allVariables)
	m = len(initFunctions)

	#calculate polynomial
	polynomial = infis[k].mul(initDenominators[k]).mul(initDenominators[k])
	for i in range(n):
		polynomial.sub(infis[i].mul(initDerivativesNum[k][i]))

	#substitute initial Functions into conditions
	polynomial = polynomial.as_expr()
	for i in range(m):
		if polynomial.has(allVariables[i]):
			polynomial = polynomial.subs(allVariables[i], initFunctions[i])

	#determine rSystem such that the coefficients vanish
	polynomial = Apoly(polynomial, allVariables, rs)
	lgs = np.empty([len(polynomial.coefs), len(rs)])
	for i in range(len(polynomial.coefs)):
		lgs[i,:] = polynomial.coefs[i]

	queue.put(lgs)

def buildSystem(numerators, denominators, derivativesNum, obsDerivativesNum,
			initDenominators, initDerivativesNum, initFunctions, 
			infis, diffInfis, allVariables, rs, nProc, ansatz):
	if nProc>1:
		from multiprocessing import Queue, Process
	else:
		from multiprocessing import Queue

	n = len(allVariables)
	m = len(numerators)
	h = len(obsDerivativesNum)
	o = len(initFunctions)

	### start the calculations for the first equations
	ns = 0
	queue = Queue()
	while ns < min([m+h+o, nProc]):
		if ns < m:
			if nProc>1: p = Process(target=doEquation, args=(ns, numerators, denominators, derivativesNum, infis,
								diffInfis, allVariables, rs, ansatz, queue))
			else: doEquation(ns, numerators, denominators, derivativesNum, infis,
								diffInfis, allVariables, rs, ansatz, queue)
		elif ns < m+h:
			if nProc>1: p = Process(target=doObsEquation, args=(ns-m, obsDerivativesNum, infis, allVariables, rs, queue))
			else: doObsEquation(ns-m, obsDerivativesNum, infis, allVariables, rs, queue)
		else:
			if nProc>1: p = Process(target=doInitEquation, args=(ns-m-h, initDenominators, initDerivativesNum,
								initFunctions, infis, allVariables, rs, queue))
			else: doInitEquation(ns-m-h, initDenominators, initDerivativesNum,
								initFunctions, infis, allVariables, rs, queue)
		if nProc>1: p.start()
		ns += 1

	sys.stdout.write("\rBuilding system...0%")
	sys.stdout.flush()

	### wait till a process has finished and start the calculation for a new equation
	lgsList = []
	lgsSize = 0
	finished = 0
	while ns < m+h+o:
		lgs = queue.get()
		if ns < m:
			if nProc>1: p = Process(target=doEquation, args=(ns,numerators, denominators, derivativesNum, infis,
								diffInfis, allVariables, rs, ansatz, queue))
			else: doEquation(ns,numerators, denominators, derivativesNum, infis,
								diffInfis, allVariables, rs, ansatz, queue)
		elif ns < m+h:
			if nProc>1: p = Process(target=doObsEquation, args=(ns-m, obsDerivativesNum, infis, allVariables, rs, queue))
			else: doObsEquation(ns-m, obsDerivativesNum, infis, allVariables, rs, queue)
		else:
			if nProc>1: p = Process(target=doInitEquation, args=(ns-m-h, initDenominators, initDerivativesNum,
								initFunctions, infis, allVariables, rs, queue))
			else: doInitEquation(ns-m-h, initDenominators, initDerivativesNum,
								initFunctions, infis, allVariables, rs, queue)
		if nProc>1: p.start()
		ns += 1

		lgsList.append(lgs)
		lgsSize += lgs.shape[0]	
		finished += 1

		prog = int(float(finished)/(m+h+o)*100)
		sys.stdout.write("\rBuilding system...%d%%" %prog)
		sys.stdout.flush()

	### wait for all processes to finish
	while finished < m+h+o:
		lgs = queue.get()

		lgsList.append(lgs)
		lgsSize += lgs.shape[0]		
		finished += 1

		prog = int(float(finished)/(m+h+o)*100)
		sys.stdout.write("\rBuilding system...%d%%" %prog)
		sys.stdout.flush()
	
	sys.stdout.write("\nCombining system...")
	sys.stdout.flush()	

	### combine all conditions into one matrix
	rSystem = np.empty([lgsSize, len(rs)])
	pos = 0
	for lgs in lgsList:
		rSystem[pos:(pos+lgs.shape[0]), :] = lgs
		pos += lgs.shape[0]

	return rSystem
		


