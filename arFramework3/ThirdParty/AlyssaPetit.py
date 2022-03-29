# AlyssaPetit version 1.0
# Use with python 3.x

import numpy
import sympy
from sympy import Matrix, simplify, expand, solve
from numpy import shape, zeros, concatenate
from numpy.linalg import matrix_rank
from sympy.parsing.sympy_parser import parse_expr
from sympy.matrices import *
from sympy.matrices import matrix_multiply_elementwise
import csv
import random
from random import shuffle

def LCS(s1, s2):
    m = [[0] * (1 + len(s2)) for i in range(1 + len(s1))]
    longest, x_longest = 0, 0
    for x in range(1, 1 + len(s1)):
        for y in range(1, 1 + len(s2)):
            if s1[x - 1] == s2[y - 1]:
                m[x][y] = m[x - 1][y - 1] + 1
                if m[x][y] > longest:
                    longest = m[x][y]
                    x_longest = x
            else:
                m[x][y] = 0
    return s1[x_longest - longest: x_longest]

def SolveSymbLES(A,b):
    dim=shape(A)[0]
    Asave=A[:]
    Asave=Matrix(dim, dim, Asave)
    #printmatrix(Asave)
    #print(b)
    determinant=Asave.det()
    if(determinant==0):
        #print('Determinant of LCL-calculation is zero! Try to specify LCLs yourself!')
        return([])
    result=[]
    for i in range(dim):
        A=Matrix(dim,dim,Asave)
        A.col_del(i)
        A=A.col_insert(i,b)
        result.append(simplify(A.det()/determinant))
    
    return(result)

def CutStringListatSymbol(liste, symbol):
    out=[]    
    for el in liste:
        if(symbol in el):
            add=el.split(symbol)
        else:
            add=[el]
        out=out+add
    return(out)

def FillwithRanNum(M):
    dimx=len(M.row(0))
    dimy=len(M.col(0))
    ranM=zeros(dimy, dimx)
    parlist=[]
    ranlist=[]
    for i in M[:]:
        if(i!=0):
            if(str(i)[0]=='-'):
                parlist.append(str(i)[1:])
            else:
                parlist.append(str(i))
    parlist=list(set(parlist))
    for symbol in [' - ', ' + ', '*', '/', '(',')']:
        parlist=CutStringListatSymbol(parlist,symbol)
    parlist=list(set(parlist))
    temp=[]    
    for i in parlist:
        if(i!=''):
            if(not is_number(i)):
                temp.append(i)
                ranlist.append(random.random())
    parlist=temp
    for i in range(dimy):
        for j in range(dimx):
            ranM[i,j]=M[i,j]
            if(ranM[i,j]!=0):
                for p in range(len(parlist)):
                   ranM[i,j]=ranM[i,j].subs(parse_expr(parlist[p]),ranlist[p])
    return(ranM)

def FindLinDep(M, tol=1e-12):
    ranM=FillwithRanNum(M)
    Q,R=numpy.linalg.qr(ranM)
    for i in range(shape(R)[0]):
        for j in range(shape(R)[1]):
            if(abs(R[i,j]) < tol):
                R[i,j]=0.0
                
    LinDepList=[]
    for i in range(shape(R)[0]):
        if(R[i][i]==0):
            LinDepList.append(i)
    
    return(LinDepList)

def FindLCL(M, X):
    LCL=[]    
    LinDepList=FindLinDep(M)
    i=0
    counter=0
    deleted_rows=[]
    states=Matrix(X[:])
    while(LinDepList!=[]):
        i=LinDepList[0]
        testM=FillwithRanNum(M)
        rowliste=list(numpy.nonzero(testM[:,i])[0])
        colliste=[i]        
        for z in range(i):
            for k in rowliste:        
                for j in range(i):
                    jliste=list(numpy.nonzero(testM[:,j])[0])
                    if(k in jliste):
                        rowliste=rowliste+jliste
                        colliste=colliste+[j]
            rowliste=list(set(rowliste))
            colliste=list(set(colliste))
        rowliste.sort()
        colliste.sort()
        colliste.pop()        
        rowlisteTry=rowliste[0:(len(colliste))]
        vec=SolveSymbLES(M[rowlisteTry,colliste],M[rowlisteTry,i])
        shufflecounter=0
        while(vec==[] and shufflecounter < 100):
            shuffle(rowliste)
            shufflecounter=shufflecounter+1
            rowlisteTry=rowliste[0:(len(colliste))]
            vec=SolveSymbLES(M[rowlisteTry,colliste],M[rowlisteTry,i])
        if(shufflecounter==100):
            print('Problems while finding conserved quantities!')
            return(0,0)
        counter=counter+1
        try:
            mat=[states[l] for l in colliste]
            test=parse_expr('0')
            for v in range(0,len(vec)):
                test=test-parse_expr(str(vec[v]))*parse_expr(str(mat[v]))
        except:
            return([],0)
        partStr=str(test)+' + '+str(states[i])
        partStr=partStr.split(' + ')
        partStr2=[]
        for index in range(len(partStr)):
            partStr2=partStr2+partStr[index].split('-')
        partStr=partStr2
        if(len(partStr) > 1):        
            CLString=LCS(str(partStr[0]),str(partStr[1]))
            for ps in range(2,len(partStr)):
                CLString=LCS(CLString,str(partStr[ps]))
        else:
            CLString=str(partStr[0])
        if(CLString==''):
            CLString=str(counter)
        LCL.append(str(test)+' + '+str(states[i])+' = '+'total'+CLString)
        M.col_del(i)
        states.row_del(i)
        deleted_rows.append(i+counter-1)
        LinDepList=FindLinDep(M)
    return(LCL, deleted_rows)

def printmatrix(M):    
    lengths=[]
    for i in range(len(M.row(0))):
        lengths.append(0)
        for j in range(len(M.col(0))):
            lengths[i]=max(lengths[i],len(str(M.col(i)[j])))          
    string=''.ljust(5)
    string2=''.ljust(5)
    for j in range(len(M.row(0))):
        string=string+(str(j)).ljust(lengths[j]+2)
        for k in range(lengths[j]+2):        
            string2=string2+('-')        
    print(string)
    print(string2)
    for i in range(len(M.col(0))):
        string=str(i).ljust(4) + '['
        for j in range(len(M.row(0))):
            if(j==len(M.row(0))-1):
                string=string+str(M.row(i)[j]).ljust(lengths[j])
            else:
                string=string+(str(M.row(i)[j])+', ').ljust(lengths[j]+2)        
        print(string+']')    
    return()
    
def printgraph(G):
    for el in G:
        print(el+': '+str(G[el]))
    return()
def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False
    
def checkNegRows(M):
    NegRows=[]
    if((M==Matrix(0,0,[])) | (M==Matrix(0,1,[])) | (M==Matrix(1,0,[]))):
        return(NegRows)
    else:        
        for i in range(len(M.col(0))):
            foundPos=False
            for j in range(len(M.row(i))):
                if(M[i,j]>0):
                    foundPos=True
            if(foundPos==False):
                NegRows.append(i)    
        return(NegRows)
    
def checkPosRows(M):
    PosRows=[]
    if((M==Matrix(0,0,[])) | (M==Matrix(0,1,[])) | (M==Matrix(1,0,[]))):
        return(PosRows)
    else: 
        for i in range(len(M.col(0))):
            foundNeg=False
            for j in range(len(M.row(i))):
                if(M[i,j]<0):
                    foundNeg=True
            if(foundNeg==False):
                PosRows.append(i)    
        return(PosRows)             

def DetermineGraphStructure(SM, F, X, neglect):
    graph={}    
    for i in range(len(SM*F)):
        liste=[]
        for j in range(len(X)):
            if((SM*F)[i]!=((SM*F)[i]).subs(X[j],1)):
                if(j==i):
                    In=((SM*F)[i]).subs(X[j],0)
                    Out=simplify(((SM*F)[i]-In)/X[j])
                    if(Out!=Out.subs(X[j],1)):
                        liste.append(str(X[j]))
                else:
                    liste.append(str(X[j]))
            else:
                if(j==i):
                    liste.append(str(X[j]))
                    
        graph[str(X[i])]=liste
        #print(graph)
    for el in neglect:
        if(parse_expr(el) in X):
            if not el in graph:
                graph[el]=[el]
            else:
                if(el not in graph[el]):
                    graph[el].append(el)
    return(graph)

def FindCycle(graph, X):
    for el in X:
        cycle=find_cycle(graph, str(el), str(el), path=[])
        if(cycle!=None):
            return(cycle)
    return(None)
    
def find_cycle(graph, start, end, path=[]):
    path = path + [start]
    if not start in graph:
        return None
    if ((start == end) & (path!=[start])):
        return path    
    for node in graph[start]:
        if node==end: return (path+[end])
        if node not in path:
            #print(node)
            newpath = find_cycle(graph, node, end, path)
            if newpath: return newpath
    return None
    
def GetBestPair(cycle, SM, fluxpars, X, LCLs, neglect):    
    for state in cycle:
        for LCL in LCLs:
            ls=parse_expr(LCL.split(' = ')[0])
            if(ls.subs(parse_expr(state),1)!=ls):
                return(0, state, None, False)
    dimList=[]
    signList=[]
    for state in cycle:
        dim, sign = GetDimension(state, X, SM, True)
        signList.append(sign)
        dimList.append(dim)
    #minOfDimList=min(dimList)
    beststate=None
    bestflux=None
    besttype=-1
    n2beat=1000
    signChanged=False
    min2beat=max(dimList)+1
    for i in range(len(dimList)):
        if(dimList[i] < min2beat):
            min2beat=dimList[i]
            sign=signList[i]
            appearList=[]
            #print(sign)
            if(sign=="minus"):
                fluxpars2use=GetNegFluxParameters(SM, fluxpars, X, cycle[i])                
            else:
                fluxpars2use=GetPosFluxParameters(SM, fluxpars, X, cycle[i])
            abort_flux=False
            for fp in fluxpars2use:
                if(str(fp) not in neglect):
                    appearList.append(GetAppearances(fp, fluxpars, SM))
                else:
                    abort_flux=True
            if(abort_flux):
                ##### Change sign
                print("Sign changed!")
                signChanged=True
                if((sign=="minus" and not signChanged) or (sign=="plus" and signChanged)):
                    fluxpars2use=GetNegFluxParameters(SM, fluxpars, X, cycle[i])                
                else:
                    fluxpars2use=GetPosFluxParameters(SM, fluxpars, X, cycle[i])
                abort_flux=False
                for fp in fluxpars2use:
                    if(str(fp) not in neglect):
                        appearList.append(GetAppearances(fp, fluxpars, SM))
                    else:
                        abort_flux=True
            if(sum(appearList) < n2beat and not abort_flux):
                n2beat=sum(appearList)
                beststate=cycle[i]
                if((sign=="minus" and not signChanged) or (sign=="plus" and signChanged)):
                    bestflux=GetNegFluxParameters(SM, fluxpars, X, cycle[i])[0]
                else:
                    bestflux=GetPosFluxParameters(SM, fluxpars, X, cycle[i])[0]
                if(min2beat==1 and max(appearList)==1):
                    besttype=1
                else:
                    if(max(appearList)==1 and min2beat>1):
                        besttype=2
                    else:
                        besttype=3
    return(besttype, beststate, bestflux, signChanged)
    
       
def GetNegFluxParameters(SM, fluxpars, X, node):
    row=list(X).index(parse_expr(node))
    liste=[]
    for i in range(len(SM.row(row))):
        if(SM.row(row)[i]<0):
            liste.append(fluxpars[i])
    return(liste)

def GetPosFluxParameters(SM, fluxpars, X, node):
    row=list(X).index(parse_expr(node))
    liste=[]    
    for i in range(len(SM.row(row))):
        if(SM.row(row)[i]>0):
            liste.append(fluxpars[i])
    return(liste)
        
def GetType(node, fp, fluxpars, LCLs):
    for LCL in LCLs:
        ls=parse_expr(LCL.split(' = ')[0])
        if(ls.subs(parse_expr(node),1)!=ls):
            return(0)
    if(GetAppearances(fp, fluxpars)==1):
        if(GetDimension(node)==1):
            return(1)
        else:
            return(2)
    else:
        return(3)
        
def GetAppearances(fp, fluxpars, SM):
    anz=0
    cols = [i for i, x in enumerate(fluxpars) if x == fp]
    #col=list(fluxpars).index(fp)
    for i in cols:
        for j in range(len(SM.col(i))):
            if(SM.col(i)[j]!=0):
                anz=anz+1
    return(anz)

def GetDimension(node, X, SM, getSign=False):
    row=list(X).index(parse_expr(node))
    anzminus=0
    anzappearminus=0
    for i in range(len(SM.row(row))):
        if(SM.row(row)[i]<0):
            anzappearminus=anzappearminus+CountNZE(SM.col(i))
            anzminus=anzminus+1
    anzplus=0
    anzappearplus=0
    for i in range(len(SM.row(row))):
        if(SM.row(row)[i]>0):
            anzappearplus=anzappearplus+CountNZE(SM.col(i))
            anzplus=anzplus+1
    if(not getSign):
        return(min(anzminus, anzplus))
    else:
        if(anzminus<anzplus or (anzminus==anzplus and anzappearminus<anzappearplus)):
            return(anzminus, "minus")
        else:
            return(anzplus, "plus")
            
def GetOutfluxes(node, X, SM, F, fluxpars):
    row=list(X).index(parse_expr(node))
    outsum=0
    out=[]
    fps=[]
    for i in range(len(SM.row(row))):
        if(SM.row(row)[i]<0):
            outsum=outsum-SM.row(row)[i]*F[i]
            out.append(-SM.row(row)[i]*F[i])
            fps.append(fluxpars[i])
    return(out, outsum, fps)

def GetInfluxes(node, X, SM, F, fluxpars):
    row=list(X).index(parse_expr(node))
    outsum=0
    out=[]
    fps=[]
    for i in range(len(SM.row(row))):
        if(SM.row(row)[i]>0):
            outsum=outsum+SM.row(row)[i]*F[i]
            out.append(SM.row(row)[i]*F[i])
            fps.append(fluxpars[i])
    return(out, outsum, fps)

def FindNodeToSolve(graph):
    for el in graph:
        if(graph[el]==[]):
            return(el)
    return(None)

def CountNZE(V):
    counter=0
    for v in V:
        if(v!=0):
            counter=counter+1
    return(counter)
    
def Sparsify(M, level, sparseIter):
    oldM=M.copy()
    if(level==3):
        ncol=len(M.row(0))
        print('0 columns of '+str(ncol) +' done')
        for i in range(ncol):            
            icol=M.col(i)
            tobeat=CountNZE(M.col(i))
            for j in range(ncol):
                if(i<j):
                    for factor_j in [1,2,-1,-2,0]:
                        for k in range(ncol):
                            if(i<k and j<k):
                                for factor_k in [1,2,-1,-2,0]:
                                    for l in range(ncol):
                                        if(i<l and j<l and k<l):
                                            for factor_l in [1,2,-1,-2,0]:
                                                test=icol+factor_j*M.col(j)+factor_k*M.col(k)+factor_l*M.col(l)
                                                if(tobeat > CountNZE(test)):
                                                    Mtest=M.copy()
                                                    Mtest.col_del(i)
                                                    Mtest=Mtest.col_insert(i,test)
                                                    if(CountNZE(test)!=0 and M.rank()==Mtest.rank()):
                                                        M=Mtest.copy()
                                                        tobeat=CountNZE(test)
                                                        #print(str(i)+'+'+str(factor_j)+'*'+str(j)+'+'+str(factor_k)+'*'+str(k)+'+'+str(factor_l)+'*'+str(l)+'    '+str(tobeat))

            print(str(i+1)+' columns of '+str(ncol) +' done')
    if(level==2):
        ncol=len(M.row(0))
        for i in range(ncol):
            icol=M.col(i)
            tobeat=CountNZE(M.col(i))
            for j in range(ncol):
                if(i<j):
                    for factor_j in [1,2,-1,-2,0]:
                        for k in range(ncol):
                            if(i<k and j<k):
                                for factor_k in [1,2,-1,-2,0]:
                                    test=icol+factor_j*M.col(j)+factor_k*M.col(k)
                                    if(tobeat > CountNZE(test)):
                                        Mtest=M.copy()
                                        Mtest.col_del(i)
                                        Mtest=Mtest.col_insert(i,test)
                                        if(CountNZE(test)!=0 and M.rank()==Mtest.rank()):
                                            M=Mtest.copy()
                                            tobeat=CountNZE(test)
                                            #print(str(i)+'+'+str(factor_j)+'*'+str(j)+'+'+str(factor_k)+'*'+str(k))
            #sys.stdout.write('\rdone %d' %i)
            #sys.stdout.flush()
            #print('\r'+str(i+1)+' columns of '+str(ncol) +' done\r')
    if(level==1):
        ncol=len(M.row(0))
        for i in range(ncol):
            icol=M.col(i)
            tobeat=CountNZE(M.col(i))
            for j in range(ncol):
                if(i<j):
                    for factor_j in [1,2,-1,-2,0]:
                        test=icol+factor_j*M.col(j)
                        if(tobeat > CountNZE(test)):
                            Mtest=M.copy()
                            Mtest.col_del(i)
                            Mtest=Mtest.col_insert(i,test)
                            if(CountNZE(test)!=0 and M.rank()==Mtest.rank()):
                                M=Mtest.copy()
                                tobeat=CountNZE(test)
    if(oldM!=M and sparseIter<10):
        oldM=M.copy() 
        print("Sparsify with level", level,", Iteration ",sparseIter, " of maximal 10")
        return(Sparsify(M,level, sparseIter=sparseIter+1))                            
    else:
        return(M)
    
def Alyssa(filename,
          injections=[],
          givenCQs=[],
          neglect=[],
          sparsifyLevel = 2,
          outputFormat='R'):
    filename=str(filename)
    file=csv.reader(open(filename), delimiter=',')
    print('Reading csv-file ...')
    L=[]
    nrrow=0
    nrcol=0
    for row in file:
        nrrow=nrrow+1
        nrcol=len(row)
        L.append(row)
        
    nrspecies=nrcol-2
    
##### Remove injections  
    counter=0
    for i in range(1,len(L)):
        if(L[i-counter][1] in injections):
            L.remove(L[i-counter])
            counter=counter+1       
    
##### Define flux vector F	
    F=[]
    
    for i in range(1,len(L)):
        F.append(L[i][1])
        #print(F)
        F[i-1]=F[i-1].replace('^','**')
        F[i-1]=parse_expr(F[i-1])
        for inj in injections:
            F[i-1]=F[i-1].subs(parse_expr(inj),0)
    F=Matrix(F)
    #print(F)
##### Define state vector X
    X=[]
    X=L[0][2:]
    for i in range(len(X)):
        X[i]=parse_expr(X[i])               
    X=Matrix(X)
    #print(X)
    Xo=X.copy()
        
##### Define stoichiometry matrix SM
    SM=[]
    for i in range(len(L)-1):
    	SM.append(L[i+1][2:])        
    for i in range(len(SM)):
    	for j in range(len(SM[0])):
    		if (SM[i][j]==''):
    			SM[i][j]='0'
    		SM[i][j]=parse_expr(SM[i][j])    
    SM=Matrix(SM)
    SM=SM.T
    SMorig=SM.copy()

    
##### Check for zero fluxes
    icounter=0
    jcounter=0
    for i in range(len(F)):
        if(F[i-icounter]==0):
            F.row_del(i-icounter)
            for j in range(len(SM.col(i-icounter))):
                if(SM[j-jcounter,i-icounter]!=0):
                    #UsedRC.append(X[j-jcounter])
                    X.row_del(j-jcounter)
                    SM.row_del(j-jcounter)
                    SMorig.row_del(j-jcounter)
                    jcounter=jcounter+1
            SM.col_del(i-icounter)
            SMorig.col_del(i-icounter)
            icounter=icounter+1
    
    print('Removed '+str(icounter)+' fluxes that are a priori zero!')
    nrspecies=nrspecies-icounter
    #printmatrix(SM)
    #print(F)
    #print(X)
    #print(UsedRC)
#####Check if some species are zero and remove them from the system
    zeroStates=[]
    NegRows=checkNegRows(SM)
    PosRows=checkPosRows(SM)
    #print(PosRows)
    #print(NegRows)
    while((NegRows!=[]) | (PosRows!=[])):
        #print(PosRows)
        #print(NegRows)
        if(NegRows!=[]):        
            row=NegRows[0]
            zeroStates.append(X[row])
            counter=0    
            for i in range(len(F)):
                if(F[i-counter].subs(X[row],1)!=F[i-counter] and F[i-counter].subs(X[row],0)==0):
                    F.row_del(i-counter)
                    SM.col_del(i-counter)                    
                    counter=counter+1
                else:
                    if(F[i-counter].subs(X[row],1)!=F[i-counter] and F[i-counter].subs(X[row],0)!=0):
                        F[i-counter]=F[i-counter].subs(X[row],0)
            X.row_del(row)
            SM.row_del(row)
        else:
            row=PosRows[0]
            zeroFluxes=[]
            for j in range(len(SM.row(row))):
                if(SM.row(row)[j]!=0):
                    zeroFluxes.append(F[j])
            for k in zeroFluxes:
                StateinFlux=[]
                for state in X:
                    if(k.subs(state,1)!=k):
                        StateinFlux.append(state)
                if(len(StateinFlux)==1):
                    zeroStates.append(StateinFlux[0])
                    row=list(X).index(StateinFlux[0])
                    counter=0            
                    for i in range(len(F)):
                        if(F[i-counter].subs(X[row],1)!=F[i-counter]):
                            if(F[i-counter].subs(X[row],0)==0):
                                F.row_del(i-counter)
                                SM.col_del(i-counter)
                            else:
                                F[i-counter]=F[i-counter].subs(X[row],0)                            
                            counter=counter+1
        #printmatrix(SM)
        NegRows=checkNegRows(SM)      
        PosRows=checkPosRows(SM)
    #printmatrix(SM)
    #print(F)
    #print(X)
    nrspecies=nrspecies-len(zeroStates)
    if(nrspecies==0):
        print('All states are zero!')
        return(0)
    else:
        if(zeroStates==[]):
            print('No states found that are a priori zero!')
        else:
            print('These states are zero:')
            for state in zeroStates:
                print('\t'+str(state))
    
    nrspecies=nrspecies+len(zeroStates)

##### Identify linearities, bilinearities and multilinearities        
    Xsquared=[]
    for i in range(len(X)):
        Xsquared.append(X[i]*X[i])        
    Xsquared=Matrix(Xsquared)
      
    BLList=[]
    MLList=[]
    for i in range(len(SM*F)):
        LHS=str(expand((SM*F)[i]))
        LHS=LHS.replace(' ','')
        LHS=LHS.replace('-','+')
        LHS=LHS.replace('**2','tothepowerof2')
        LHS=LHS.replace('**3','tothepowerof3')
        exprList=LHS.split('+')
        for expr in exprList:
            VarList=expr.split('*')
            counter=0
            factors=[]
            for j in range(len(X)):
                anz=0
                if(str(X[j]) in VarList):
                    anz=1
                    factors.append(X[j])
                if((str(X[j])+'tothepowerof2') in VarList):
                    anz=2 
                    factors.append(X[j])
                    factors.append(X[j])
                if((str(X[j])+'tothepowerof3') in VarList):
                    anz=3
                    factors.append(X[j])
                    factors.append(X[j])
                    factors.append(X[j])
                counter=counter+anz
            if(counter==2):
                string=''            
                for l in range(len(factors)):
                    if(l==len(factors)-1):
                        string=string+str(factors[l])
                    else:
                        string=string+str(factors[l])+'*'
                if(not(string in BLList)):
                    BLList.append(string)
            if(counter>2):
                string=''            
                for l in range(len(factors)):
                    if(l==len(factors)-1):
                        string=string+str(factors[l])
                    else:
                        string=string+str(factors[l])+'*'
                if(not(string in MLList)):
                    MLList.append(string)
        
    COPlusLIPlusBL=[]
    for i in range(len(SM*F)):
        COPlusLIPlusBL.append((SM*F)[i])
        for j in range(len(MLList)):
            ToSubs=expand((SM*F)[i]).coeff(MLList[j])
            COPlusLIPlusBL[i]=expand(COPlusLIPlusBL[i]-ToSubs*parse_expr(MLList[j]))
            
    COPlusLI=[]
    for i in range(len(COPlusLIPlusBL)):
        COPlusLI.append(COPlusLIPlusBL[i])
        for j in range(len(BLList)):
            ToSubs=expand((COPlusLIPlusBL)[i]).coeff(BLList[j])
            COPlusLI[i]=expand(COPlusLI[i]-ToSubs*parse_expr(BLList[j]))
    
##### C*X contains linear terms
    C=zeros(len(COPlusLI),len(X))  
    for i in range(len(COPlusLI)):
    	for j in range(len(X)):
    		C[i*len(X)+j]=expand((COPlusLI)[i]).coeff(X[j])
        
##### ML contains multilinearities
    ML=expand(Matrix(SM*F)-Matrix(COPlusLIPlusBL))
##### BL contains bilinearities
    BL=expand(Matrix(COPlusLIPlusBL)-Matrix(COPlusLI))    
#### CM is coefficient matrix of linearities
    CM=C        
#####CMBL gives coefficient matrix of bilinearities
    CMBL=[]
    if(BLList!=[]):
        for i in range(len(BLList)):
            CVBL=[]
            for k in range(len(BL)):
                CVBL.append(BL[k].coeff(BLList[i]))
            CMBL.append(CVBL)            
    else:
        CVBL=[]
        for k in range(len(BL)):
            CVBL.append(0)
        CMBL.append(CVBL)
    
    CMBL=Matrix(CMBL).T 
    
#####CMML gives coefficient matrix of multilinearities
#####Summarize multilinearities and bilinearities 
    if(MLList!=[]):
        CMML=[]
        for i in range(len(MLList)):
            CVML=[]
            for k in range(len(ML)):
                CVML.append(expand(ML[k]).coeff(MLList[i]))
            CMML.append(CVML)    
        CMML=Matrix(CMML).T  
        BLList=BLList+MLList
        CMBL=Matrix(concatenate((CMBL,CMML),axis=1))
      
    for i in range(len(BLList)):
        BLList[i]=parse_expr(BLList[i])
       
    if(BLList!=[]):    
        CMbig=Matrix(concatenate((CM,CMBL),axis=1))
    else:
        CMbig=Matrix(CM)      

#### Save ODE equations for testing solutions at the end    
    print('Rank of SM is '+str(SM.rank()) + '!')
    SMorig=SM.copy()
    ODE=SMorig*F
#### Get Flux Parameters
    fluxpars=[]
    for flux in F:
        if(flux.args!=()):
            foundFluxpar=False
            for el in flux.args:
                if(not foundFluxpar and el not in X and not is_number(str(el))):
                    if(flux.subs(el, 0)==0):
                        fluxpars.append(el)
                        foundFluxpar=True
        else:
            fluxpars.append(flux)

##### Increase Sparsity of stoichiometry matrix SM
    print('Sparsify stoichiometry matrix with sparsify-level '+str(sparsifyLevel)+'!')
    newSM=(Sparsify(SM.T, level=sparsifyLevel, sparseIter=1)).T
    if(newSM!=SM):
        print("Sparsified!")
        SM=newSM
    
#### Find conserved quantities
    
    #printmatrix(CMbig)
    #print(X)
    if(givenCQs==[]):
        print('\nFinding conserved quantities ...')
        LCLs, rowsToDel=FindLCL(CMbig.transpose(), X)
    else:
        print('\nI took the given conserved quantities!')
        LCLs=givenCQs
    if(LCLs!=[]):
        print(LCLs)
    else:
        print('System has no conserved quantities!')
#### Define graph structure
    print('\nDefine graph structure ...\n')
    
    SSgraph=DetermineGraphStructure(SM, F, X, neglect)    
    #printgraph(SSgraph)
    #print(fluxpars)
#### Check for Cycles
    cycle=FindCycle(SSgraph, X)
#### Remove cycles step by step
    gesnew=0
    eqOut=[]
    while(cycle!=None):
        print('Removing cycle '+str(counter))
        #printmatrix(SM)
        #print(F)
        minType, state2Rem, fp2Rem, signChanged = GetBestPair(cycle, SM, fluxpars, X, LCLs, neglect)
        #print(cycle)
        #print(state2Rem)
        #print(fp2Rem)
        #print(minType)
        if(minType==-1):
            print("    The cycle")
            print("       "+str(cycle))
            print("    cannot be removed. Set more parameters free or enable steady-state expressions with minus signs. The latter is not yet provided by the tool.")
            return(0)
        if(minType==0):
            for LCL in LCLs:
                ls=parse_expr(LCL.split(' = ')[0])
                if(ls.subs(parse_expr(state2Rem),1)!=ls):
                    LCL2Rem=LCL
            LCLs.remove(LCL2Rem)
            index=list(X).index(parse_expr(state2Rem))
            eqOut.append(state2Rem+' = '+state2Rem)
            print('   '+str(state2Rem)+' --> '+'Done by CQ')
        if(minType==1):
            index=list(X).index(parse_expr(state2Rem))
            eq=(SM*F)[index]
            sol=solve(eq, fp2Rem, simplify=False)[0]
            eqOut.append(str(fp2Rem)+' = '+str(sol))
            print('   '+str(state2Rem)+' --> '+str(fp2Rem))
        if(minType==2):
            anz, sign=GetDimension(state2Rem, X, SM, getSign=True)
            index=list(X).index(parse_expr(state2Rem))            
            negs, sumnegs, negfps=GetOutfluxes(state2Rem, X, SM, F, fluxpars)
            poss, sumposs, posfps=GetInfluxes(state2Rem, X, SM, F, fluxpars)
            if(anz==1):
                print("Error in Type Determination. Please report this bug!")
                return(0)
            else:
                nenner=1
                for j in range(anz):
                    if(j>0):
                        nenner=nenner+parse_expr('r_'+state2Rem+'_'+str(j))
                trafoList=[]
                if((sign=="minus" and not signChanged) or (sign=="plus" and signChanged)):
                    for j in range(len(negs)):
                        flux=negs[j]
                        fp=negfps[j]
                        prefactor=flux/fp
                        if(j==0):
                            trafoList.append(str(fp)+' = ('+str(sumposs)+')*1/('+str(nenner)+')*1/('+str(prefactor)+')')
                        else:
                            gesnew=gesnew+1
                            trafoList.append(str(fp)+' = ('+str(sumposs)+')*'+'r_'+state2Rem+'_'+str(j)+'/('+str(nenner)+')*1/('+str(prefactor)+')')                        
                    print('   '+str(state2Rem)+' --> '+str(negfps))
                    
                else:
                    for j in range(len(poss)):
                        flux=poss[j]
                        fp=posfps[j]
                        prefactor=flux/fp
                        if(j==0):
                            trafoList.append(str(fp)+' = ('+str(sumnegs)+')*1/('+str(nenner)+')*1/('+str(prefactor)+')')
                        else:
                            gesnew=gesnew+1
                            trafoList.append(str(fp)+' = ('+str(sumnegs)+')*'+'r_'+state2Rem+'_'+str(j)+'/('+str(nenner)+')*1/('+str(prefactor)+')')
                    print('   '+str(state2Rem)+' --> '+str(posfps))
                for eq in trafoList:
                    eqOut.append(eq)
        if(minType==3):
            anz, sign=GetDimension(state2Rem, X, SM, getSign=True)
            index=list(X).index(parse_expr(state2Rem))            
            negs, sumnegs, negfps=GetOutfluxes(state2Rem, X, SM, F, fluxpars)
            poss, sumposs, posfps=GetInfluxes(state2Rem, X, SM, F, fluxpars)
            if(anz==1):
                if((sign=="minus" and not signChanged) or (sign=="plus" and signChanged)):
                    fp2Rem=negfps[0]
                    flux=negs[0]
                else:
                    fp2Rem=posfps[0]
                    flux=poss[0]
                eq=(SM*F)[index]
                sol=solve(eq, fp2Rem, simplify=False)[0]
                eqOut.append(str(fp2Rem)+' = '+str(sol))
                FsearchFlux = matrix_multiply_elementwise(abs(SM[index,:]),F.T)
                colindex=list(FsearchFlux).index(flux)
                for row2repl in range(len(SM.col(0))):
                    if(SM[row2repl,colindex]!=0 and row2repl!=index):
                        SM=SM.row_insert(row2repl,SM.row(row2repl)-(SM[row2repl,colindex]/SM[index,colindex])*SM.row(index))
                        SM.row_del(row2repl+1)
            else:
                nenner=1
                for j in range(anz):
                    if(j>0):
                        nenner=nenner+parse_expr('r_'+state2Rem+'_'+str(j))
                trafoList=[]
                if((sign=="minus" and not signChanged) or (sign=="plus" and signChanged)):
                    for j in range(len(negs)):
                        flux=negs[j]
                        fp=negfps[j]
                        prefactor=flux/fp
                        if(j==0):
                            trafoList.append(str(fp)+' = ('+str(sumposs)+')*1/('+str(nenner)+')*1/('+str(prefactor)+')')
                        else:
                            gesnew=gesnew+1
                            trafoList.append(str(fp)+' = ('+str(sumposs)+')*'+'r_'+state2Rem+'_'+str(j)+'/('+str(nenner)+')*1/('+str(prefactor)+')')
                        
                        FsearchFlux = matrix_multiply_elementwise(abs(SM[index,:]),F.T)
                        colindex=list(FsearchFlux).index(flux)
                        for k in range(len(posfps)):
                            SM=SM.col_insert(len(SM.row(0)),SM.col(colindex))
                            F=F.row_insert(len(F),Matrix(1,1,[poss[k]/nenner]))
                            fluxpars.append(posfps[k])
                        SM.col_del(colindex)
                        F.row_del(colindex)
                        fluxpars.__delitem__(colindex)
                    print('   '+str(state2Rem)+' --> '+str(negfps))
                    
                else:
                    for j in range(len(poss)):
                        flux=poss[j]
                        fp=posfps[j]
                        prefactor=flux/fp
                        if(j==0):
                            trafoList.append(str(fp)+' = ('+str(sumnegs)+')*1/('+str(nenner)+')*1/('+str(prefactor)+')')
                        else:
                            gesnew=gesnew+1
                            trafoList.append(str(fp)+' = ('+str(sumnegs)+')*'+'r_'+state2Rem+'_'+str(j)+'/('+str(nenner)+')*1/('+str(prefactor)+')')
                        FsearchFlux = matrix_multiply_elementwise(abs(SM[index,:]),F.T)
                        colindex=list(FsearchFlux).index(flux)
                        for k in range(len(negfps)):
                            SM=SM.col_insert(len(SM.row(0)),SM.col(colindex))
                            F=F.row_insert(len(F),Matrix(1,1,[negs[k]/nenner]))
                            fluxpars.append(negfps[k])
                        SM.col_del(colindex)
                        F.row_del(colindex)
                        fluxpars.__delitem__(colindex)
                    print('   '+str(state2Rem)+' --> '+str(posfps))
                for eq in trafoList:
                    eqOut.append(eq)
        X.row_del(index)
        SM.row_del(index)
        SSgraph=DetermineGraphStructure(SM, F, X, neglect)
        #print(X)
        #printgraph(SSgraph)
        cycle=FindCycle(SSgraph, X)
        counter=counter+1       
    print('There is no cycle in the system!\n')              
    
#### Solve remaining equations
    eqOut.reverse()
    print('Solving remaining equations ...\n')
    while(SSgraph!={}):
        #print(SSgraph)
        node=FindNodeToSolve(SSgraph)
        #print(node)        
        index=list(X).index(parse_expr(node))
        #print((SM*F)[index])
        sol=solve((SM*F)[index],parse_expr(node), simplify=True)
        #print(sol)
        eqOut.insert(0,node+' = '+str(sol[0]))
        for f in range(len(F)):                
            F[f]=F[f].subs(parse_expr(node), sol[0])
            #print(node+' = '+str(sol[0]))
        X.row_del(index)
        SM.row_del(index)
        SSgraph=DetermineGraphStructure(SM, F, X, neglect=[])
    
#### Test Solution  
    print('Testing Steady State...\n')
    NonSteady=False
    #print(eqOut)
    #print(ODE)
    #print(SM*F)
    for i in range(len(ODE)):
        expr=parse_expr(str(ODE[i]))
        for j in range(len(zeroStates)):
            zeroState=zeroStates[j]
            expr=expr.subs(zeroState, 0)
        #print(len(eqOut))
        for j in range(len(eqOut)):
            ls, rs = eqOut[-(j+1)].split('=')
            #print(ls)
            ls=parse_expr(ls)
            #print(rs)
            rs=parse_expr(rs)
            expr=expr.subs(ls, rs)
            #print(simplify(expr))
        expr=simplify(expr)
        #print(expr)
        if(expr!=0):
            print('   Equation '+str(ODE[i]))
            print('   results:'+str(expr))
            NonSteady=True
    if(NonSteady):
        print('Solution is wrong!\n')
    else:
        print('Solution is correct!\n')
    
#### Print Equations
    print('I obtained the following equations:\n')
    if(outputFormat=='M'):
        for state in zeroStates:
            print('\tinit_'+str(state)+'  "0"'+'\n')
        eqOutReturn=[]
        for i in range(len(eqOut)):
            ls, rs = eqOut[i].split('=')
            ls=parse_expr(ls)
            rs=parse_expr(rs)
            for j in range(i,len(eqOut)):
                ls2, rs2 = eqOut[j].split('=')
                rs2=parse_expr(rs2)
                rs2=rs2.subs(ls,rs)
                eqOut[j]=str(ls2)+'='+str(rs2)
            for state in Xo:
                ls=ls.subs(state, parse_expr('init_'+str(state)))
                rs=rs.subs(state, parse_expr('init_'+str(state)))
            eqOut[i]=str(ls)+'  "'+str(rs)+'"'
                            
        for i in range(len(eqOut)):
            eqOut[i]=eqOut[i].replace('**','^')
                    
        for eq in eqOut:
            print('\t'+eq+'\n')
            eqOutReturn.append(eq)            
        
    else:
        for state in zeroStates:
            print('\t'+str(state)+' = 0'+'\n')
        eqOutReturn=[]
        for eq in eqOut:            
            ls, rs = eq.split(' = ')
            print('\t'+ls+' = "'+rs+'",'+'\n')
            eqOutReturn.append(ls+'='+rs)
    print('Number of Species:  '+str(nrspecies))
    print('Number of Equations:  '+str(len(eqOut)+len(zeroStates)))
    print('Number of new introduced variables:  '+str(gesnew))
    return(eqOutReturn)
