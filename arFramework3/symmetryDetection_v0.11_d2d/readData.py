# Author: Benjamin Merkt, Physikalisches Institut, Universitaet Freiburg

import csv
import sys

import sympy as spy

#from sympy.parsing.sympy_parser import parse_expr
from sympy.parsing.sympy_tokenize import tokenize

# try/except necessary for R interface (imports automatically and does not find other files)
try:
	from functions import *
except:
	pass

def readModel(fileName, delimT):
	if delimT == 't':
		delim = '\t'
	else:
		delim = delimT
	variables = []
	parameters = []
	flows = []
	stoichiometry = []

	global l; l = -1
	
	with open(fileName, 'rb') as defFile:
		reader = csv.reader(defFile, delimiter=delim, quoting=csv.QUOTE_NONE)
		
		row = reader.next()
		for i in range(2,len(row)):
			row[i] = row[i].replace('"','')
			variables.append(giveVar(row[i]))
		
		lines = 0 
		stoichiometryList = []
		for row in reader:
			row[1] = row[1].replace('"','')
			row[1] = row[1].replace('^','**')
			flows.append(row[1])
			
			lines += 1
			for i in range(2,len(row)):
				if row[i] == '': num = 0
				else: 
					row[i] = row[i].replace('"','')
					if row[i] == '':
						row[i] = 0
					num = int(row[i])
				stoichiometryList.append(num)

		stoichiometryT = spy.Matrix(lines,len(variables),stoichiometryList)
		stoichiometry = stoichiometryT.transpose()

		def read():
			global l
			l += 1
			if l >= len(flows): raise StopIteration
			else: return flows[l]

		def useToken(key, value, Coord1, Coord2, fullLine):
			if key == 1:
				parameters.append(giveVar(value))

		tokenize(read,useToken)#get parameters from flows

	parameters = sorted(list(set(parameters)), key=spy.default_sort_key)

	for entry in variables:
		if entry in parameters:
			parameters.remove(entry)

	for f in range(len(flows)):
		flows[f] = giveParsed(flows[f])
		
	return variables, parameters, spy.Matrix(len(flows),1,flows), stoichiometry



def readEquations(equationSource):
	if isinstance(equationSource, basestring):
		eq_file = open(equationSource,'r')
		def read():
			line = eq_file.readline()
			line = line.replace('"','').replace(',','')
			return line
	else:
		global l
		l = 0
		def read():
			global l
			if l == len(equationSource): raise StopIteration
			line = equationSource[l]
			line = line.replace('"','').replace(',','').strip()
			l += 1
			return line + '\n'
	
	global newLine; newLine = True
	global variables; variables = []
	global functions; functions = []
	global parameters; parameters = []
	
	def useToken(key, value, Coord1, Coord2, fullLine):
		global newLine, variables, obsFunctions, parameters
		if key == 1: #1: NAME  2: NUMBER  51: OP   4: NEWLINE  0: ENDMARKER
			if newLine == True:
				variables.append(giveVar(value))
				functions.append(giveParsed(fullLine[(fullLine.find('=')+1):len(fullLine)]))
			else:
				parameters.append(giveVar(value))
			newLine = False
		elif key == 4:
			newLine = True
		
	tokenize(read,useToken)
	
	parameters = sorted(list(set(parameters)), key=spy.default_sort_key)
	
	for entry in variables:
		if entry in parameters:
			parameters.remove(entry)
			
	return variables, functions, parameters



def readObservation(observation_path, variables, parameters):

    observables, obsFunctions, obsParameters = readEquations(observation_path)

	#remove dynamic parameters and variables from observation Parameters
    for var in variables:
		if var in obsParameters:
			obsParameters.remove(var)
    for par in parameters:
		if par in obsParameters:
			obsParameters.remove(par)

    return observables, obsFunctions, parameters+obsParameters



def readInitialValues(initial_path, variables, parameters):

	initVars, initFunctions, initParameters = readEquations(initial_path)
	o = len(initVars)
	m = len(variables)

    #remove variables and other parameters for initParameters
	i = 0
	while i < len(initParameters):
		if initParameters[i] in variables+parameters:
			initParameters.pop(i)
		else:
			i += 1

    #if variabel not restricted, introduce inital value parameter
	for i in range(o):
		if initVars[i] == initFunctions[i]:
			initFunctions[i] = giveVar(str(initVars[i])+'_0')
			initParameters.append(initFunctions[i])

    #subsitute dependence of other variables
	substituted = True
	counter = 0
	while substituted:
		substituted = False
		for k in range(o):
			for j in range(o):
				if initFunctions[k].has(initVars[j]):
					initFunctions[k] = initFunctions[k].subs(initVars[j],initFunctions[j])
					substituted = True
		counter += 1
		if counter > 100:
			raise(UserWarning('There seems to be an infinite recursion in the initial value functions'))

    #order varaibels according to equations
	initFunctionsOrdered = [0]*m
	for i in range(m):
		try: initFunctionsOrdered[i] = initFunctions[initVars.index(variables[i])]
		except ValueError: #if not contained introduce new unconstrained parameter
 			initFunctionsOrdered[i] = giveVar(str(variables[i])+'_0')
			initParameters.append(initFunctionsOrdered[i])

	return initFunctionsOrdered, parameters+initParameters



def readPredictions(prediction_path, variables, parameters,):

    predictions, predFunctions, predParameters = readEquations(prediction_path)

    i = 0
    while i < len(predParameters):
	    if predParameters[i] in variables+parameters:
		    predParameters.pop(i)
	    else:
		    i += 1
    if len(predParameters) != 0:
	    raise(UserWarning('Error: New parameters occured in predictions: ' + str(predParameters)))

    return predictions, predFunctions

