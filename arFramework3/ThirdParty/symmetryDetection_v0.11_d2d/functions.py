# Author: Benjamin Merkt, Physikalisches Institut, Universitaet Freiburg

import sys
import time

import numpy as np
import sympy as spy

from sympy.parsing.sympy_parser import parse_expr

# try/except necessary for R interface (imports automatically and does not find other files)
try:
	from polyClass import *
except:
	pass

# readline might not be available
try:
	import readline
	readlineAvailable = True
except:
	readlineAvailable = False
	
extension_str = '_93502158393064762'

class bcolors:
	purple = '\033[95m'
	blue = '\033[94m'
	green = '\033[92m'
	yellow = '\033[93m'
	red = '\033[91m'
	endc = '\033[0m'
	bold = '\033[1m'
	underline = '\033[4m'
	
# wrapper on spy.var for renaming QCOSINE variables
def giveVar(expr):
	if expr == 'epsilon':
		raise(UserWarning("Transformation parameter 'epsilon' not allowed in any input"))
	if expr == 'lamda':
		raise(UserWarning("Transformation parameter 'lamda' not allowed in any input"))
	if expr == 'T':
		raise(UserWarning("Variable 'T' reserved for time"))
		
	for v in ['Q', 'C', 'O', 'S', 'I', 'N', 'E']:
		expr = expr.replace(v, v + extension_str)
		
	return spy.var(expr)
	
# wrapper on sympy.parsing.sympy_parser.parse_expr for renaming QCOSINE variables
def giveParsed(expr):
	for v in ['Q', 'C', 'O', 'S', 'I', 'N', 'E']:
		expr = expr.replace(v, v + extension_str)
		
	return parse_expr(expr)
	
	
# recursive function to construct a multidimensional polynomial
# vars: variables, i: position in vars, p: degree left for other variables
# summand: current monom under construction, poly: full polynomial
# num: umber of coefficients, k: ansatz for which variable, rs: list of coefficiets
def giveDegree(vars, i, p, summand, poly, num, k, rs):

	if i == len(vars)-1:
		rs.append(giveVar('r_'+str(vars[k])+'_'+str(num)))
		poly += rs[-1]*summand*vars[i]**p
		return poly, num+1
	else:
		for j in range(p+1):
			poly, num = giveDegree(vars, i+1, p-j, summand*vars[i]**j, poly, num, k, rs)

	return poly, num

# make infinitesimal ansatz
def makeAnsatz(ansatz, allVariables, m, q, pMax, fixed):
	n = len(allVariables)

	if ansatz == 'uni':
		#construct polynomial
		rs = []
		infis = []
		for k in range(n):
			infis.append(spy.sympify(0))
			if allVariables[k] in fixed: continue #if in fixed, ansatz is 0
			for p in range(pMax+1):
				rs.append(giveVar('r_'+str(allVariables[k])+'_'+str(p)))
				infis[-1] += rs[-1] * allVariables[k]**p

		#calculate derivatives
		diffInfis = [[0]*n]
		for i in range(n):
			diffInfis[0][i] = spy.diff(infis[i],allVariables[i])

	elif ansatz == 'par':
		rs = []
		infis = []
		for k in range(n):
			infis.append(spy.sympify(0))
			if allVariables[k] in fixed: continue #if in fixed, ansatz is 0
			num = 0
			for p in range(pMax+1): #for every degree for 0 to pMax
				vari = allVariables[m+q:] #all parameters
				if k < (m+q): #if ansatz is not for a 
					vari.append(allVariables[k])
					kp = len(vari)-1
				else:	
					kp = k-(m+q)
				degree, num = giveDegree(vari, 0, p, 1, 0, num, kp, rs)
				infis[-1] += degree

		#calculate derivatives
		diffInfis = [[0]*n]
		for i in range(n):
			diffInfis[0][i] = spy.diff(infis[i], allVariables[i])

	elif ansatz == 'multi':
		rs = []
		infis = []
		for k in range(n):
			infis.append(spy.sympify(0))
			if allVariables[k] in fixed: continue #if in fixed, ansatz is 0
			num = 0
			for p in range(pMax+1): #for every degree for 0 to pMax
				if k < m:  #if ansatz is for a dynamic variable
					vari = allVariables[:m] + allVariables[m+q:]
					kp = k
				elif k < m+q:  #if ansatz is for an input
					vari = allVariables[:]
					kp = k
				else:  #if ansatz is for a parameter	
					vari = allVariables[m+q:]  #all parameters
					kp = k-(m+q)
				degree, num = giveDegree(vari, 0, p, 1, 0, num, kp, rs)
				infis[-1] += degree

		#calculate derivatives
		diffInfis = [[0]*n for i in range(n)]
		for i in range(n):
			for j in range(n):
				diffInfis[i][j] = spy.diff(infis[i],allVariables[j])

	return infis, diffInfis, rs

def transformExprToPoly(diff, i, infis, queue, allVariables, rs):
	if diff:
		queue.put((Apoly(infis[i[0]][i[1]], allVariables, rs), diff, i))
	else:
		queue.put((Apoly(infis[i], allVariables, rs), diff, i))

def transformInfisToPoly(infis, diffInfis, allVariables, rs, nProc, ansatz):
	if nProc > 1:
		from multiprocessing import Queue, Process
	else:
		from multiprocessing import Queue

	n = len(infis)
	k = len(diffInfis)

	ns = 0	
	queue = Queue()
	### start the transformation for the first equations
	while ns < min([n+k*n, nProc]):
		if ns < n:
			if nProc > 1: p = Process(target=transformExprToPoly, args=(False, ns, infis, queue, 
										allVariables, rs))
			else: transformExprToPoly(False, ns, infis, queue, allVariables, rs)
		else:
			if ansatz == 'multi': i = divmod(ns-n,n)		
			else: i = (0, ns-n)
			if nProc > 1: p = Process(target=transformExprToPoly, args=(True, i, diffInfis, queue, 
										allVariables, rs))
			else: transformExprToPoly(True, i, diffInfis, queue, allVariables, rs)
		if nProc > 1: p.start()
		ns += 1

	sys.stdout.write("\rPreparing equations...0%")
	sys.stdout.flush()

	### wait till a process has finished and start the transformation for a new equation
	infisPoly = [0]*n
	diffInfisPoly = [0]*k
	for i in range(k):
		diffInfisPoly[i] = [0]*n
	finished = 0
	while ns < n+k*n:
		#if mp:
		poly, diff, i = queue.get()
		if diff: diffInfisPoly[i[0]][i[1]] = poly
		else: infisPoly[i] = poly
		finished += 1

		if ns < n:
			if nProc > 1: p = Process(target=transformExprToPoly, args=(False, ns, infis, queue, 
										allVariables, rs))
			else: transformExprToPoly(False, ns, infis, queue, allVariables, rs)
		else:
			if ansatz == 'multi':  i = divmod(ns-n,n)		
			else: i = (0, ns-n)
			if nProc > 1: p = Process(target=transformExprToPoly, args=(True, i, diffInfis, queue, 
										allVariables, rs))
			else: transformExprToPoly(True, i, diffInfis, queue, allVariables, rs)
		if nProc > 1: p.start()
		ns += 1

		prog = int(float(finished)/(n+k*n)*100)
		sys.stdout.write("\rPreparing equations...%d%%" %prog)
		sys.stdout.flush()

	### wait for all processes to finish
	while finished < n+k*n:
		poly, diff, i = queue.get()
		if diff: diffInfisPoly[i[0]][i[1]] = poly
		else: infisPoly[i] = poly
		finished += 1

		prog = int(float(finished)/(n+k*n)*100)
		sys.stdout.write("\rPreparing equations...%d%%" %prog)
		sys.stdout.flush()

	return infisPoly, diffInfisPoly


### calculate rref from a upper triangular matrix
def getrref(rSystem):
	pivots = []
	pivotLines = []
	i = -1
	for j in xrange(rSystem.shape[1]):
		if rSystem[j,j] == 0:
			k = 1
			while j-k > i:
				if rSystem[j-k,j] != 0:
					i = j-k
					break
				k += 1
			else:
				k = i-1
				while k >= 0:
					if rSystem[k,j] != 0 and (not k in pivotLines):
						rSystem[[j,k],:] = rSystem[[k,j],:]
						i = j
						break
					k -= 1
				else: continue
		else:
			i = j

		pivots.append(j)
		pivotLines.append(i)

		coeff = rSystem[i,j]
		rSystem[i,:] = rSystem[i,:]/coeff

		for k in xrange(i):
			coeff = rSystem[k,j]
			if coeff != 0:
				rSystem[k,:] = rSystem[k,:] - coeff*rSystem[i,:]

	return rSystem[pivotLines,:], pivots


### returns a matrix of base vectors of the null space given a matrix in rref
### the base vectors are the columns of the matrix
def nullSpace(matrix, pivots):
	m = matrix.shape[1]

	notPivots = []
	solutions = np.zeros((m, m-len(pivots)))

	i, k, l = m-1, 0, matrix.shape[0]-1
	while i >= 0:
		if i in pivots:
			for h in range(len(notPivots)):
				solutions[i,h] = - matrix[l,notPivots[h]]
			l -= 1
		else:	
			notPivots.append(i)
			solutions[i,k] = 1
			k += 1
		i -= 1
	
	return solutions

def checkForCommonFactor(infisTmp, allVariables, m):
	spy.var('epsilon')
	# extract all factors from first infinitesimal
	for i in range(len(infisTmp)):
		if infisTmp[i] != 0:
			fac = spy.factor(infisTmp[i])
			if type(fac) == type(epsilon+1):
				factors = [infisTmp[i]]
			elif type(fac) == type(epsilon):
				factors = [fac]
			else:
				factors = list(fac.args)
		
			break

	# remove all non-parameters from factors list
	i = 0
	while i < len(factors):
		if factors[i].is_number:
			factors.pop(i)
		elif factors[i] in allVariables[:m]:
			factors.pop(i)
		elif type(factors[i]) == type(epsilon+1):
			factors.pop(i)				
		elif type(factors[i]) == type(epsilon**2):
			if type(factors[i].args[0]) != type(epsilon+1):
				factors[i] = factors[i].args[0]
				i += 1
			else:
				factors.pop(i)					
		else:
			i += 1

	#check which of the factors is in all other infinitesimals
	for i in range(1,len(infisTmp)):
		if infisTmp[i] == 0: continue
		fac = spy.factor(infisTmp[i])
		if type(fac) == type(epsilon+1):
			factorsTmp = [fac]
		elif type(fac) == type(epsilon):
			factorsTmp = [fac]
		else:
			factorsTmp = list(fac.args)

		j = 0
		while j < len(factors):
			k = 0
			while k < len(factorsTmp):
				if factorsTmp[k].is_number:
					factorsTmp.pop(k)
				elif factorsTmp[k] in allVariables[:m]:
					factorsTmp.pop(k)
				elif type(factorsTmp[k]) == type(epsilon+1):
					factorsTmp.pop(k)				
				elif type(factorsTmp[k]) == type(epsilon**2):
					if type(factorsTmp[k].args[0]) != type(epsilon+1):
						factorsTmp[k] = factorsTmp[k].args[0]
						k += 1
					else:
						factorsTmp.pop(k)					
				else:
					k += 1
			if factors[j] in factorsTmp:
				j += 1
				continue
			else:
				factors.pop(j)

		if len(factors) != 0:
			continue	#if potential common factors are left, try next ifinitesimal
		else:
			break		#otherwise treat next solution

	if len(factors) == 0:
		return False
	else:
		return True

### determine known transformations from infinitesimals
def buildTransformation(infis, allVariables):
	n = len(infis)
	spy.var('epsilon')
	spy.var('lamda')
	
	transformations = [0]*n
	tType = np.zeros(6,dtype=bool)  #0: unknown, 
									#1: scaling, 
									#2: translation, 
									#3: MM-like, 
									#4: p>2, 
									#5: generalized translation
	for i in range(n):
		if infis[i] == 0: # not transformed
			transformations[i] = allVariables[i]
		else:
			poly = spy.Poly(infis[i], allVariables).as_dict()
			monomials = poly.keys()
			coefs = poly.values()
			if len(monomials) == 1:
				p = None
				for j in range(n):
					if monomials[0][j] != 0:
						if j == i and p == None: # p Symmetry
							p = monomials[0][i]
						elif p == None and monomials[0][j] == 1: # 
							p = -1-j							
						else:
							transformations[i] = '-?-'
							tType[0] = True
							break
				else:
					if p == None: # translation
						transformations[i] = allVariables[i] + epsilon*coefs[0]
						tType[2] = True
					elif p <= 0: # generalized translation
						transformations[i] = allVariables[i] + epsilon*coefs[0] * \
												allVariables[-p-1]
						tType[5] = True
					elif p == 1: # scaling
						#transformations[i] = spy.exp(epsilon*coefs[0]) * allVariables[i]
						transformations[i] = lamda**coefs[0] * allVariables[i]
						tType[1] = True
					else: # p Symmetry
						transformations[i] = spy.simplify(allVariables[i]/(1-(p-1) * epsilon * \
												allVariables[i]**(p-1))**(spy.sympify(1)/(p-1)))
						if p == 2: tType[3] = True
						else: tType[4] = True
			else:	
				transformations[i] = '-?-'
				tType[0] = True

	if tType[1] and not np.all(tType == np.array([0,1,0,0,0,0], dtype=np.bool)):
		transformations = map(lambda trans: trans.subs(lamda, spy.exp(epsilon)), transformations)

	string = 'Type: '
	if tType[0]: string += 'unknown, '
	if tType[1]: string += 'scaling, '
	if tType[2]: string += 'translation, '
	if tType[3]: string += 'MM-like, '
	if tType[4]: string += 'p>2, '
	if tType[5]: string += 'gen. tanslation, '

	string = string[0:(len(string)-2)]

	return transformations, string

### print found transformations
def printTransformations(infisAll, allVariables, pretty, suffix):
	st = str(len(infisAll)) + ' transformation(s) found:\n'
 	if pretty:
 		st = bcolors.purple + bcolors.bold + st + bcolors.endc
 	print st

	n = len(infisAll[0])

	fhandle = open( suffix + '_result.txt', 'w');
	
	length1 = 8
	length2 = 13
	length3 = 14
	transformations = [0]*len(infisAll)
	types = [0]*len(infisAll)
	outputs = []
	for l in range(len(infisAll)):
		for i in range(n):
			infisAll[l][i] = spy.nsimplify(infisAll[l][i])
		transformations[l], types[l] = buildTransformation(infisAll[l], allVariables)
		
		outputs.append([])
		for i in range(n):
			if infisAll[l][i] != 0:
				# get stuff for output line
				t1, t2, t3 = allVariables[i], infisAll[l][i], transformations[l][i]
				for u in range(len(allVariables)):
					v_str = str(allVariables[u])
					for j in ['Q', 'C', 'O', 'S', 'I', 'N', 'E']:
						v_str = v_str.replace(j + extension_str, j)
					v = spy.var(v_str)
					t1 = t1.subs(allVariables[u], v)
					t2 = t2.subs(allVariables[u], v)
					if not isinstance(t3, basestring): 
						t3 = t3.subs(allVariables[u], v)
				if pretty:
					t1, t2, t3 = spy.pretty(t1), spy.pretty(t2), spy.pretty(t3)
				else:
					t1, t2, t3 = str(t1), str(t2), str(t3)
				outputs[-1].append([t1, t2, t3])

				# remove string extension
				#for v in ['Q', 'C', 'O', 'S', 'I', 'N', 'E']:
					#outputs[-1][-1][0] = outputs[-1][-1][0].replace(v + extension_str, v)
					#outputs[-1][-1][1] = outputs[-1][-1][1].replace(v + extension_str, v)
					#outputs[-1][-1][2] = outputs[-1][-1][2].replace(v + extension_str, v)
				
				# search for longest string
				if len(outputs[-1][-1][0]) > length1:
					length1 = len(outputs[-1][-1][0])					
				if len(outputs[-1][-1][1]) > length2:
					length2 = len(outputs[-1][-1][1])					
				if len(outputs[-1][-1][2]) > length3:
					length3 = len(outputs[-1][-1][2])

	# print all stuff
	if pretty:
		print bcolors.bold + bcolors.purple + ('{0:'+str(length1)+'s}').format('Variable') + bcolors.endc + ' : '\
			+ bcolors.bold + bcolors.purple + ('{0:'+str(length2)+'s}').format('Infinitesimal') + bcolors.endc + ' : '\
			+ bcolors.bold + bcolors.purple + str('Transformation') + bcolors.endc

		for l in range(len(infisAll)):
			print '='*(length1+length2+14+6)
			print bcolors.blue + bcolors.bold + '#' + str(l+1) + ': ' + types[l] + bcolors.endc#+'\n'
			
			for i, lst in enumerate(outputs[l]):
				print '-'*(length1+length2+14+6)
				str1 = lst[0].splitlines()
				str2 = lst[1].splitlines()
				str3 = lst[2].splitlines()
				lines = max([len(str1), len(str2), len(str3)])

				for j in range(lines):
					if j >= len(str1):	s1 = ''
					else:				s1 = str1[j]
					if j >= len(str2):	s2 = ''
					else:				s2 = str2[j]
					if j >= len(str3):	s3 = ''
					else:				s3 = str3[j]
					print bcolors.bold + s1 + ' '*(length1-len(s1)) + bcolors.endc + ' : ' + \
							s2 + ' '*(length2-len(s2)) + ' : ' +\
							s3 + ' '*(length3-len(s3))


				# print bcolors.bold + lst[0] + ' '*(length1-len(lst[0])) + bcolors.endc + ' : '\
				# 		+ lst[1] + ' '*(length2-len(lst[1])) + ' : '\
				# 		+ lst[2].replace('\n','\n'+' '*length1 + ' : ' + ' '*length2 + ' : ')
				#print ' '*length1 + ' : ' + ' '*length2 + ' : '

	else:
		print ('{0:'+str(length1)+'s} : ').format('variable') \
			+ ('{0:'+str(length2)+'s} : ').format('infinitesimal')\
			+ str('transformation')

		for l in range(len(infisAll)):
			print '-'*(length1+length2+length3+6)
			print '#' + str(l+1) + ': ' + types[l]
			
			for lst in outputs[l]:
				print ('{0:'+str(length1)+'s} : ').format(lst[0]) \
						+ ('{0:'+str(length2)+'s} : ').format(lst[1])\
						+ lst[2]

			fhandle.write('#' + str(l+1) + '\t' + types[l] + '\t')
			fhandle.write(str(lst[0]))
			fhandle.write('\t')
			fhandle.write(str(lst[2]))
			fhandle.write('\n')

	fhandle.close()						


class bcolors:
	purple = '\033[95m'
	blue = '\033[94m'
	green = '\033[92m'
	yellow = '\033[93m'
	red = '\033[91m'
	endc = '\033[0m'
	bold = '\033[1m'
	underline = '\033[4m'