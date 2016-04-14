# Author: Benjamin Merkt, Physikalisches Institut, Universitaet Freiburg
# Version: 0.11

import sys
import argparse
import time

import sympy as spy
import scipy.linalg

# try/except necessary for R interface which imports automatically after loading
try:
	from readData import *
	from functions import *
	from buildSystem import *
	from polyClass import *
	from checkPredictions import *
except:
	pass

t0 = time.time()

spy.var('epsilon')
spy.var('T')

def symmetryDetection(allVariables, diffEquations, observables, obsFunctions, initFunctions,
						predictions, predFunctions, ansatz = 'uni', pMax = 2, inputs = [], 
						fixed = [], parallel = 1, allTrafos = False, timeTrans = False,
						pretty = True, suffix=''):
	
	n = len(allVariables)
	m = len(diffEquations)
	h = len(observables)

	###########################################################################################
	#############################     prepare equations    ####################################
	###########################################################################################
	sys.stdout.write('Preparing equations...')
	sys.stdout.flush()

	# make infinitesimal ansatz
	infis, diffInfis, rs = makeAnsatz(ansatz, allVariables, m, len(inputs), pMax, fixed)

	# get infinitesimals of time transformation
	if timeTrans:
		rs.append(spy.var('r_T_1'))
		diffInfiT = rs[-1]
		allVariables += [T]
	else:
		diffInfiT = None

	# and convert to polynomial
	infis, diffInfis = transformInfisToPoly(infis, diffInfis, allVariables, rs, parallel, ansatz)
	
	diffInfiT = Apoly(diffInfiT, allVariables, rs)

	### extract numerator and denominator of equations
	#differential equations
	numerators = [0]*m
	denominators = [0]*m
	for k in range(m):
		rational = spy.together(diffEquations[k])
		numerators[k] = Apoly(spy.numer(rational), allVariables, None)
		denominators[k] = Apoly(spy.denom(rational), allVariables, None)

	#observation functions
	obsNumerators = [0]*h
	obsDenominatros = [0]*h
	for k in range(h):
		rational = spy.together(obsFunctions[k])
		obsNumerators[k] = Apoly(spy.numer(rational), allVariables, None)
		obsDenominatros[k] = Apoly(spy.denom(rational), allVariables, None)

	#initial functions
	if len(initFunctions) != 0:
		initNumerators = [0]*m
		initDenominatros = [0]*m
		for k in range(m):
			rational = spy.together(initFunctions[k])
			initNumerators[k] = Apoly(spy.numer(rational), allVariables, None)
			initDenominatros[k] = Apoly(spy.denom(rational), allVariables, None)
	else:
		initNumerators = []
		initDenominatros = []

	### calculate numerator of derivatives of equations
	#differential equatioins
	derivativesNum = [0]*m
	for i in range(m):
		derivativesNum[i] = [0]*n
	for k in range(m):
		for l in range(n):
			derivativesNum[k][l] = Apoly(None, allVariables, None)
			derivativesNum[k][l].add(numerators[k].diff(l).mul(denominators[k]))
			derivativesNum[k][l].sub(numerators[k].mul(denominators[k].diff(l)))

	#observation functions
	obsDerivativesNum = [0]*h
	for i in range(h):
		obsDerivativesNum[i] = [0]*n
	for k in range(h):
		for l in range(n):
			obsDerivativesNum[k][l] = Apoly(None, allVariables, None)
			obsDerivativesNum[k][l].add(obsNumerators[k].diff(l).mul(obsDenominatros[k]))
			obsDerivativesNum[k][l].sub(obsNumerators[k].mul(obsDenominatros[k].diff(l)))

	#initial functions
	if len(initFunctions) != 0:
		initDerivativesNum = [0]*len(initFunctions)
		for i in range(m):
			initDerivativesNum[i] = [0]*n
		for k in range(m):
			for l in range(n):
				initDerivativesNum[k][l] = Apoly(None, allVariables, None)
				initDerivativesNum[k][l].add(initNumerators[k].diff(l).mul(initDenominatros[k]))
				initDerivativesNum[k][l].sub(initNumerators[k].mul(initDenominatros[k].diff(l)))
	else:
		initDerivativesNum = []

	sys.stdout.write('\rPreparing equations...done\n')
	sys.stdout.flush()

	###########################################################################################
	############################     build linear system    ###################################
	###########################################################################################
	sys.stdout.write('\nBuilding system...')
	sys.stdout.flush()

	rSystem = buildSystem(numerators, denominators, derivativesNum, obsDerivativesNum,
				initDenominatros, initDerivativesNum, initFunctions, 
				infis, diffInfis, diffInfiT, allVariables, rs, parallel, ansatz)

	sys.stdout.write('done\n')
	sys.stdout.flush()

	
	###########################################################################################
	##############################     solve system    ########################################
	###########################################################################################
	sys.stdout.write('\nSolving system of size ' + str(rSystem.shape[0]) + 'x' +\
						str(rSystem.shape[1]) + '...')
	sys.stdout.flush()

	#get LU decomposition from scipy
	rSystem = scipy.linalg.lu(rSystem, permute_l=True)[1]

	#calculate reduced row echelon form
	rSystem, pivots = getrref(rSystem)

	sys.stdout.write('done\n')
	sys.stdout.flush()

	###########################################################################################
	#############################     process results    ######################################
	###########################################################################################
	sys.stdout.write('\nProcessing results...')
	sys.stdout.flush()

	# calculate solution space
	sys.stdout.write('\n  calculating solution space')
	sys.stdout.flush()
	baseMatrix = nullSpace(rSystem, pivots)

	#substitute solutions into infinitesimals
	#(and remove the ones with common parameter factors)
	sys.stdout.write('\n  substituting solutions')
	sys.stdout.flush()
	infisAll = []
	for l in range(baseMatrix.shape[1]):
		infisTmp = [0]*n
		for i in range(len(infis)):
			infisTmp[i] = infis[i].getCopy()
			infisTmp[i].rs = baseMatrix[:,l]
			infisTmp[i] = infisTmp[i].as_expr()
		if timeTrans:
			infisTmp.append(baseMatrix[-1,l] * T)

		if allTrafos:
			infisAll.append(infisTmp)
		else:
			if not checkForCommonFactor(infisTmp, allVariables, m):
				infisAll.append(infisTmp)
			
	print('')
	sys.stdout.write('done\n')
	sys.stdout.flush()

	# print transformations
	print('\n\n')
	if len(infisAll) != 0: printTransformations(infisAll, allVariables, pretty,suffix)

	###########################################################################################
	############################     check predictions    #####################################
	###########################################################################################
	if predictions != False:
		checkPredictions(predictions, predFunctions, infisAll, allVariables)

	print(time.strftime('\nTotal time: %Hh:%Mm:%Ss', time.gmtime(time.time()-t0)))


###########################################################################################
###################     main (start program from terminal)    #############################
###########################################################################################
def main():	

	# check if run with arguments (i.e. from terminal)
	try:
		sys.argv[0]
	except:
		return
		
	parser = argparse.ArgumentParser(usage='%(prog)s model_path observation_path [prediction_path] [options]', description='Detect symmetries in systems of ODEs.')
	parser.add_argument('model_path', help = 'model csv-file with path')
	parser.add_argument('observation_path', help = 'observation txt-file with path')
	parser.add_argument('prediction_path', nargs='?', default=False, 
						help = 'prediction txt-file with path (optional)')
	parser.add_argument('-I','--initial', nargs = 1, default=[False], 
						help = 'initial values txt-file with path')
	parser.add_argument('-d','--delim', nargs = 1, default = [','], 
						help = 'delimiter used in the model csv (default = ,)')
	parser.add_argument('-a','--ansatz', choices=['uni', 'par', 'multi'], default = 'uni', 
						help='ansatz made for infinitesimals (default = uni)')
	parser.add_argument('-p','--pMax', nargs = 1, default = [2], type = int, 
						help = 'maximal power used in the infinitesimal generator (default = 2)')
	parser.add_argument('-i','--input', nargs = '+', default = [], 
						help = 'input variables')
	parser.add_argument('-f','--fixed', nargs = '+', default = [], 
						help = 'variables to consider fixed')
	parser.add_argument('-P','--parallel', nargs = 1, default=[1], 
						help = 'maximal number of processes (default = 1)')
	parser.add_argument('-t','--timeTrans', action='store_true', default=False, 
						help = 'allow scaling transformations of time variable')
	parser.add_argument('-A','--allTrafos', action='store_true', default=False, 
						help = 'do not remove transformations with common parameter factors')
	parser.add_argument('--notPretty', action='store_false', default=True, 
						help = 'do not use pretty printing for output of transformations')

	args = parser.parse_args()
    
	inputs = args.input
    #read and print input and fixed variables
	if len(inputs) != 0:
		s = 'Input variables: '
		for v in range(len(inputs)):
			s = s + str(inputs[v]) + ', '
			inputs[v] = giveVar(inputs[v])
		sys.stdout.write(s[0:len(s)-2] + '\n')
		sys.stdout.flush()
	fixed = args.fixed
	if len(fixed) != 0:
		s = 'Fixed variables: '
		for v in range(len(fixed)):
			s = s + str(fixed[v]) + ', '
			fixed[v] = giveVar(fixed[v])
		sys.stdout.write(s[0:len(s)-2] + '\n')
		sys.stdout.flush()


	####### read data from files #######
	sys.stdout.write('\nReading files...')
	sys.stdout.flush()

	# read model
	variables, parameters, flows, stoichiometry = readModel(args.model_path, args.delim[0])
	suffix = str(args.model_path)
	print('%%%' + suffix)

	# read observation
	observables, obsFunctions, parameters = readObservation(args.observation_path, variables,
															parameters)


	# read initial values
	if args.initial[0] != False:
		initFunctions, parameters = readInitialValues(args.initial[0], variables, parameters)
	else:
		initFunctions = []

	# read predictions
	if args.prediction_path != False:
		predictions, predFunctions = readPredictions(args.prediction_path, variables, parameters)
	else:
		predictions, predFunctions = False, False

	# remove inputs from parameters
	for par in inputs:
		if par in parameters:
			parameters.remove(par)

	#define some stuff
	diffEquations = stoichiometry * flows
	allVariables = variables + inputs + parameters

	sys.stdout.write('done\n')
	sys.stdout.flush()
        
	symmetryDetection(allVariables, diffEquations, observables, obsFunctions, initFunctions,
						 predictions, predFunctions, args.ansatz, args.pMax[0], args.input,
						 args.fixed, int(args.parallel[0]), args.allTrafos, args.timeTrans,
						 args.notPretty, suffix)



###########################################################################################
#########################     start program from dMod   ###################################
###########################################################################################
def symmetryDetectiondMod(model, observation, prediction, initial, ansatz, pMax, inputs, fixed,
							parallel, allTrafos, timeTrans=False):
	'''start program from dMod'''

	if model == None:
		model = []
	elif isinstance(model, basestring):
		model = [model]
		
	if observation == None:
		observation = []
	elif isinstance(observation, basestring):
		observation = [observation]
		
	if prediction == None:
		prediction = []
	elif isinstance(prediction, basestring):
		prediction = [prediction]
		
	if initial == None:
		initial = []
	elif isinstance(initial, basestring):
		initial = [initial]	

	if fixed == None:
		fixed = []
	elif isinstance(fixed, basestring):
		fixed = [str(fixed)]
	if len(fixed) != 0:
		s = 'Fixed variables: '
		for v in range(len(fixed)):
			s = s + str(fixed[v]) + ', '
			fixed[v] = giveVar(fixed[v])
		sys.stdout.write(s[0:len(s)-2] + '\n')
		sys.stdout.flush()

	if inputs == None:
		inputs = []
	elif isinstance(inputs, basestring):
		inputs = [str(inputs)]
	if len(inputs) != 0:
		s = 'Input variables: '
		for v in range(len(inputs)):
			s = s + str(inputs[v]) + ', '
			inputs[v] = giveVar(inputs[v])
		sys.stdout.write(s[0:len(s)-2] + '\n')
		sys.stdout.flush()
	
	sys.stdout.write('\nReading input...')
	sys.stdout.flush()

	# read model
	variables, diffEquations, parameters = readEquations(model)

	# read observation
	observables, obsFunctions, parameters = readObservation(observation, variables, parameters)

	# read initial values
	if len(initial) != 0:
		initFunctions, parameters = readInitialValues(initial, variables, parameters)
	else:
		initFunctions = []

	# read predictions
	if len(prediction) != 0:
		predictions, predFunctions = readPredictions(prediction, variables, parameters)
	else:
		predictions, predFunctions = False, False

	# remove inputs from parameters
	for par in inputs:
		if par in parameters:
			parameters.remove(par)

	#define some stuff
	allVariables = variables + inputs + parameters
	
	sys.stdout.write('done\n')
	sys.stdout.flush()
	
	symmetryDetection(allVariables, diffEquations, observables, obsFunctions, initFunctions,
						predictions, predFunctions,	ansatz, pMax, inputs, fixed, parallel, 
						allTrafos, timeTrans=timeTrans, pretty=True)

if __name__ == "__main__":
	main()
		
