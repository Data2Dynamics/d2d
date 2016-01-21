# Author: Benjamin Merkt, Physikalisches Institut, Universitaet Freiburg

import sys

import sympy as spy

# try/except necessary for R interface (imports automatically and does not find other files)
try:
	from functions import extension_str
except:
	pass

def checkPredictions(predictions, predFunctions, infisAll, allVariables):
	n = len(allVariables)
	print '\nChecking predictions:'

	printStrings = []
	for i in range(len(predictions)):
		printStrings.append([])

		admits = True
		for j in range(len(infisAll)):
			infiPred = 0
			for k in range(n):
				if infisAll[j][k] != 0:
					infiPred += infisAll[j][k] * spy.diff(predFunctions[i], allVariables[k])
	
			infiPred = spy.simplify(infiPred)
			if infiPred != 0:
				admits = False
				p = spy.Wild('p',exclude=[0])
				c = spy.Wild('c')
				if infiPred.match(c*predFunctions[i]**p) != None:
					matches = infiPred.match(c*predFunctions[i]**p)
					printStrings[i].append([\
							str(predictions[i]),
							'#'+str(j+1),
							str((c*predictions[i]**p).subs(c,matches[c]).subs(p,matches[p]))])
				elif infiPred.match(c*(-1*predFunctions[i])**p) != None:
					matches = infiPred.match(c*(-1*predFunctions[i])**p)
					printStrings[i].append([\
							str(predictions[i]),							
							'#'+str(j+1),
							str((c*(-1)**p*predictions[i]**p).subs(c,matches[c]).subs(p,matches[p]))])
				else:
					printStrings[i].append([str(predictions[i]), '#'+str(j+1), str(infiPred)])

		if admits:
			printStrings[i] = True

	length0 = 10
	length1 = 10
	length2 = 13
	for i in range(len(printStrings)):
		tmp = str(predictions[i])
		for v in ['Q', 'C', 'O', 'S', 'I', 'N', 'E']:
			tmp = tmp.replace(v + extension_str, v)	
					
		if length0 <= len(tmp):
			length0 = len(tmp)
		if printStrings[i] == True: continue
		for j in range(len(printStrings[i])):
		
			for v in ['Q', 'C', 'O', 'S', 'I', 'N', 'E']:
					printStrings[i][j][0] = printStrings[i][j][0].replace(v + extension_str, v)	
					printStrings[i][j][2] = printStrings[i][j][2].replace(v + extension_str, v)		
		
			if length1 <= len(printStrings[i][j][1]):
				length1 = len(printStrings[i][j][1])
			if length2 <= len(printStrings[i][j][2]):
				length2 = len(printStrings[i][j][2])

	print ('{0:'+str(length0)+'s} : ').format('prediction') \
		+ ('{0:'+str(length1)+'s} : ').format('symmetry')\
		+ str('infinitesimal')
	for i in range(len(predictions)):
		print '-'*(length0+length1+length2+6)
		if printStrings[i] == True:
			print ('{0:'+str(length0)+'s} : ').format(tmp) \
				+ ('{0:'+str(length1)+'s} : ').format('admits all')\
				+ ('{0:'+str(length2)+'s}').format('      -      ')
			continue
			
		print ('{0:'+str(length0)+'s} : ').format(printStrings[i][0][0]) \
			+ ('{0:'+str(length1)+'s} : ').format(printStrings[i][0][1])\
			+ str(printStrings[i][0][2])
		for j in range(1,len(printStrings[i])):
			print ('{0:'+str(length0)+'s} : ').format('') \
				+ ('{0:'+str(length1)+'s} : ').format(printStrings[i][j][1])\
				+ str(printStrings[i][j][2])
