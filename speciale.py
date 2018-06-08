import numpy as np
import sys
import math
import sympy
import re
import pstats
import cProfile
import time
from subprocess import call
from sympy.abc  import epsilon
from fractions  import Fraction
np.set_printoptions(threshold=np.nan)
sys.setrecursionlimit(3000)

#Example LCP
M0 = np.array([[2,1], [1,2]])
b0 = np.array([[-5, -6]])
#Solution:
#  (w1,w2;z1,z2) = (0, 0; 4/3, 7/3)

#Example 2
M1 = np.array([[3,2,1], [2,3,1], [1,2,3]])
b1 = np.array([[-5,-6,-7]])

#Example 3
M2 = np.array([[1,2,3,4], [2,3,4,1], [3,4,1,2], [4,1,2,3]])
b2 = np.array([[-4,-5,-6,-7]])



class SystemInfo:

    #PivotInfo
    eRow     = 0
    eColumn  = 0
    eValue   = 0
    leftBase = 0
    
    detB     = 1
    pev      = 1
    avb      = np.array([])

    def __init__(self, M, b):
        self.B   = np.arange(len(M)).astype(np.int64)
        self.I   = np.identity(len(M)).astype(np.int64)
        self.M   = M
        self.b   = b
        self.avb = np.concatenate((self.I, -np.ones((1,len(M))).astype(np.int64).T, -M, b, self.I), axis=1)

    def getColumnsFrombAndEpsilon(self):
        return self.avb[:, -len(self.I)-len(self.b[0]):-1]
    
    def getEnteringColumn(self):
        return self.avb[:,self.eColumn]

    def prep(self):
        b1        = self.avb[:,-len(self.I)-len(self.b[0]):]
        v         = np.zeros([2*len(self.I)+1,b1.shape[1]])
        v[self.B] = b1
        return v

    def getw(self):
        v = self.prep()
        return v[0:len(self.I),:]

    def getz(self):
        v = self.prep()
        return v[len(self.I)+1:2*len(self.I)+1,:]

    def getz0(self):
        v = self.prep()
        return v[len(self.I),:]

    def getRay(self):
        ec                     = -self.getEnteringColumn()
        b1                     = self.avb[:,-len(self.I)-2:]/self.detB
        vtilde                 = np.zeros([2*len(self.I)+1,b1.shape[1]])
        vtilde[self.B,0]       = ec/self.detB
        vtilde[self.eColumn,0] = 1
        wtilde                 = vtilde[0:len(self.I),:]
        ztilde                 = vtilde[len(self.I)+1:2*len(self.I)+1,:]
        z0tilde                = vtilde[len(self.I),:]
        return wtilde,ztilde,z0tilde
    

#Give it the list of all candidates and the corresponding potential pivot values. Returns any vectors that are true
#candidates divided by their corresponding pivot element.
def findEnteringRow(SystemInfo):
    candidates    = SystemInfo.getColumnsFrombAndEpsilon().tolist()
    pivotValues   = SystemInfo.getEnteringColumn()
    newCandidates = []
    zero          = [0 for i in range(len(candidates[0]))] 
    #print('findEntering')
    for i in range(len(candidates)):
        if(pivotValues[i] != 0):
            candFrac = [Fraction(int(candidates[i][j]), int(pivotValues[i])) for j in range(len(candidates[i]))]
            if(candFrac > zero):
                candFrac.extend((pivotValues[i],i))
                newCandidates.append(candFrac)
    if(len(newCandidates) == 0):
        return ' '
    return min(newCandidates)

def LemkeMain(M,b):
    system1 = SystemInfo(M,b)

    candidates      = system1.getColumnsFrombAndEpsilon().tolist()
    system1.eRow    = candidates.index(min(candidates))
    system1.eColumn = len(system1.I)
    system1.eValue  = system1.avb[int(system1.eRow)][int(system1.eColumn)]

    while(True):
        LemkePivoting(system1)

        #d=np.ones((1,len(M)))
        #w=system1.getw()
        #z=system1.getz()
        #z0=system1.getz0()
        #Mz=M.dot(z)
        #dz0=np.outer(d,z0)

        #print((w-dz0-Mz)[:,0:2])
    
    
        #leftBase == lengthOfI hvis z0 er leaving variable.
        if(system1.leftBase == len(system1.I)):
            print('Succes, vectors:')
            return system1
        
        if(system1.leftBase > len(system1.I)):
            system1.eColumn = system1.leftBase - len(system1.I) - 1
        else:
            system1.eColumn = system1.leftBase + len(system1.I) + 1
        
        winner = findEnteringRow(system1)
        if(winner == ' '):
            with open('fejl.txt', 'w') as f4:
                f4.writelines(str(M))
                f4.writelines(str(b))
            return system1
    #if(len(candidates) == 0):
    #    wtilde,ztilde,z0tilde=system1.getRay()
    #    w1=w+wtilde
    #    z1=z+ztilde
    #    z01=z0+z0tilde
    #    Mz1=M.dot(z1)
    #    dz01=np.outer(d,z01)
    #    print('1 unit out of ray')
    #    print((w1-dz01-Mz1)[:,0:2])
    #    zz=ztilde[:,0]
    #    x=zz[0:4]
    #    y=zz[4:8]
    #    p=zz[8:11]-zz[11:14]
    #    g=zz[14:16]
    #    q=zz[16:19]-zz[19:22]
    #    h=zz[22:24]
    #    break
        system1.eRow   = winner[-1]
        system1.eValue = system1.avb[int(system1.eRow)][int(system1.eColumn)]

                
def LemkePivoting(SystemInfo):
    for i in range(len(SystemInfo.avb)):
        if(i != SystemInfo.eRow):
            #Her ganges pivot elementets værdi på række i != pivotrækken
            SystemInfo.avb[i]   = SystemInfo.eValue*SystemInfo.avb[i]
            
            #Her findes forholdet mellem, elementet i samme søjle som pivotelementet i række i, og pivotelementet.
            ratio               = SystemInfo.avb[i][int(SystemInfo.eColumn)]//SystemInfo.eValue
            
            #Her trækkes et heltals multiplum af pivotrækken fra række i.
            SystemInfo.avb[i]  -= ratio*SystemInfo.avb[int(SystemInfo.eRow)]
            
            #Her divideres med den forhenværende pivot.
            SystemInfo.avb[i] //= SystemInfo.pev
            
    #Opdater determinant og sidste pivot.
    SystemInfo.detB  *= SystemInfo.eValue
    SystemInfo.detB //= SystemInfo.pev
    SystemInfo.pev    = SystemInfo.eValue
    
    #Opdater basen efter pivot og gem den gamle base til at finde complementet.
    SystemInfo.leftBase                = SystemInfo.B[int(SystemInfo.eRow)]
    SystemInfo.B[int(SystemInfo.eRow)] = SystemInfo.eColumn
    
    return SystemInfo


#Part II: Sorting networks
def insertionSort(nrOfVar):
    sortingNetwork = []
    previousLayer  = []
    for i in range(nrOfVar-1):
        sortingNetwork.append((i,i+1))
        sortingNetwork.extend(reversed(previousLayer))
        previousLayer.append((i,i+1))
    return sortingNetwork


def Batcher(nrOfVar):
    if (nrOfVar != 0 and ((nrOfVar & (nrOfVar - 1)) == 0)):
        sortingNetwork = list(splitInTwoAndSort(0, nrOfVar-1))
    else:
        closestHigherPowerOfTwo = 2**math.ceil(math.log(nrOfVar, 2))
        sortingNetwork          = list(splitInTwoAndSort(0,closestHigherPowerOfTwo-1))
        sortingNetwork          = [(car,cdr) for car,cdr in sortingNetwork if car < nrOfVar and cdr < nrOfVar]
    return sortingNetwork


def splitInTwoAndSort(minVar, maxVar):
    if (maxVar - minVar >= 1):
        midVar = int((maxVar-minVar)/2) + minVar
        yield from splitInTwoAndSort(minVar, midVar)
        yield from splitInTwoAndSort(midVar + 1, maxVar)
        yield from mergeSortedLists(minVar, maxVar, 1)

def mergeSortedLists(minVar, maxVar, distance):
    nextDistance = distance * 2
    if (nextDistance < maxVar - minVar):
        yield from mergeSortedLists(minVar, maxVar, nextDistance)
        yield from mergeSortedLists(minVar + distance, maxVar, nextDistance)
        yield from [(i, i + distance) for i in range(minVar + distance, maxVar - distance, nextDistance)]
    else:
        yield (minVar, minVar + distance) 


#Part III: Conversion to extended formulation.
def makeVariableTable(nrOfVar):
    table = {}
    for i in range(nrOfVar):
        table[str(i)] = i
    return table
    
def setupMatrices(network, varNr):
    rowsInE = len(network) + varNr
    G       = np.zeros((len(network)*2,len(network)*2+varNr)).astype(np.int64)
    E       = np.zeros((rowsInE,len(network)*2+varNr)).astype(np.int64)
    e       = np.zeros((rowsInE, varNr)).astype(np.int64)
    return makeConstraints(varNr, G, E, e, network)

def makeConstraints(nrOfVar, matrixG, matrixE, matrixe, sortingNetwork):
    varTable    = makeVariableTable(nrOfVar)
    indexNumber = nrOfVar
    for i in range(0, 2*len(sortingNetwork), 2):
        #Check if the variables in the current pair have been previously updated and apply the update if so.
        firstVarInPair  = varTable[str(sortingNetwork[int(i/2)][0])]
        secondVarInPair = varTable[str(sortingNetwork[int(i/2)][1])]
    
        #Set the next two rows not yet modified of G, to the values given by the constraints computed over the table.
        #The constraints are always made such that the network pair numbers are compared to the smallest index number. 
        #So these entries in the matrix gets the value 1.
        matrixG[i][indexNumber]            =  1
        matrixG[i + 1][indexNumber]        =  1

        #Modifying G according to the 2nd part of the constraints. 
        matrixG[i][firstVarInPair]         = -1
        matrixG[i + 1][secondVarInPair]    = -1
    
        #Modifying E according to the equations.
        matrixE[int(i/2)][firstVarInPair]  =  1
        matrixE[int(i/2)][secondVarInPair] =  1
        matrixE[int(i/2)][indexNumber]     = -1
        matrixE[int(i/2)][indexNumber + 1] = -1
    
        #Update the indexes in the table for the next iteration.
        varTable[str(sortingNetwork[int(i/2)][0])] = indexNumber
        varTable[str(sortingNetwork[int(i/2)][1])] = indexNumber + 1
        indexNumber += 2
    
    #Update matrixE with epsilon equations. VarTable has been used to keep track of which variable the initial ones have changed to. 
    #That way we know which variable each epsilon should equal.
    for i in range(len(varTable)):
        matrixE[len(sortingNetwork) + i][varTable[str(i)]]  = 1    
        matrixe[len(sortingNetwork) + i][len(varTable)-i-1] = 1 
    return matrixG, matrixE, matrixe



#Part IV: More conversion.
def makeMAndbMatrices(A, B, E, F, G, H, e, f):
    
    #Matrices A and B are padded with zeroes to fit inside their respective cells.
    paddedB                            = np.zeros((len(F.T), len(G[0]))).astype(np.int64)
    paddedB[0:len(B.T), 0:len(B.T[0])] = -B.T
    paddedA                            = np.zeros((len(G.T), len(F[0]))).astype(np.int64)
    paddedA[0:len(A), 0:len(A[0])]     = -A
    
    #Making b
    if(len(e[0]) != len(f[0])):
        smallerOne = e
        if(len(e[0]) > len(f[0])):
            smallerOne = f
        difference       = np.absolute(len(e[0]) - len(f[0]))
        differenceZeroes = np.zeros((len(smallerOne), difference)).astype(np.int64)
        smallerOne       = np.concatenate((smallerOne, differenceZeroes), axis=1)
        if(len(e) > len(f)):
            f = smallerOne
        else:
            e = smallerOne
        
    columnsOfb = np.maximum(len(e[0]), len(f[0]))
    b1         = np.zeros((len(paddedA) + len(paddedB), columnsOfb)).astype(np.int64)
    b2         = np.zeros((len(G), columnsOfb)).astype(np.int64)
    b3         = np.zeros((len(H), columnsOfb)).astype(np.int64)
    b          = np.vstack((b1, e, -e, b2, f, -f, b3))

    
    #Making M
    #E and G have the same number of columns. The number of rows in M11 is dependant on the number of rows in E^T/G^T,
    #which is the number of columns in E/G. The number of columns in M11 are dependant on the number of columns in E/G.
    M11        = np.zeros((len(G.T), len(G[0]))).astype(np.int64)
    M61        = np.zeros((2*len(F)+len(H),len(G[0]))).astype(np.int64)
    
    #Rows and columns are created one by one to merge later.
    columnOne  = np.vstack((M11,paddedB, -E, E, -G, M61))
    
    M22To52    = np.zeros((len(paddedB)+2*len(E)+len(G), len(F[0]))).astype(np.int64)
    columnTwo  = np.vstack((paddedA, M22To52, -F, F, -H))

    M16        = np.zeros((len(G.T), 2*len(F.T[0])+len(H.T[0]))).astype(np.int64)
    rowOne     = np.concatenate((E.T, -E.T, G.T, M16), axis=1)
    
    M22To25    = np.zeros((len(F.T), 2*len(E.T[0])+len(G.T[0]))).astype(np.int64)
    rowTwo     = np.concatenate((M22To25, F.T, -F.T, H.T), axis=1)
    
    dimOfRest  = 2*len(E)+2*len(F)+len(G)+len(H)
    MRest      = np.zeros((dimOfRest, dimOfRest)).astype(np.int64)
    
    #Putting the pieces together
    lastColumn = np.vstack((rowOne, rowTwo, MRest))
    M          = np.concatenate((columnOne, columnTwo, lastColumn), axis=1)
    return M, b

def makeMatrixNeg(M):
    M -= np.amax(M) + 1
    return M

#Part V: Put the pieces together.
def findProperEquilibrium(A,B):
    print('Matrix A:\n', A, '\nMatrix B:\n', B)
    A          = makeMatrixNeg(A)
    B          = makeMatrixNeg(B)
    networkOne = Batcher(len(A))
    G, E, e    = setupMatrices(networkOne, len(A))
    H, F, f    = G, E, e
    if(len(A) != len(A[0])):
        networkTwo = Batcher(len(A[0]))
        H, F, f    = setupMatrices(networkTwo, len(A[0]))
    M, b        = makeMAndbMatrices(A,B,E,F,G,H,e,f)
    print(len(M),'\n',len(b))
    #system1     = LemkeMain(M,b)
    allEpsilons = system1.getz()[:,0:len(b[0])]
    epsilonsX   = allEpsilons[0:len(A), :]
    epsilonsY   = allEpsilons[len(G[0]): len(G[0]) + len(A[0]), :]
    xVector     = 'Proper Equilibrium\nx = ('
    xVector    += makeOutput(epsilonsX, system1) + '),'
    yVector     = 'y = ('
    yVector    += makeOutput(epsilonsY, system1) + ').'
    output      = xVector + '\n' + yVector
    return output

def zeroOrNegInMatrix(M):
    M += -np.amin(M) + 1
    return M


def findEquilibrium(A,B):
    #print('Matrix A:\n', A, '\nMatrix B:\n', B)
    A             = zeroOrNegInMatrix(A)
    B             = zeroOrNegInMatrix(B)
    b             = -np.ones([len(A)+len(B.T), 1]).astype(np.int64)
    #b[0:len(A)]  *= -1
    matrixZeroes1 = np.zeros([len(A),len(B.T[0])]).astype(np.int64)
    matrixZeroes2 = np.zeros([len(B.T), len(A[0])]).astype(np.int64)
    firstRow      = np.concatenate((matrixZeroes1,A),axis=1)
    secondRow     = np.concatenate((B.T, matrixZeroes2),axis=1)
    C             = np.vstack((firstRow,secondRow))
    C            -= np.amax(C) + 1
    #print('Matrix C:\n',C)
    M             = -C
    #print('Matrix M:\n',M,'\nMatrix b:\n', b)
    system1       = LemkeMain(M,b)
    allEpsilons   = system1.getz()[:,0:1]
    epsilonsX     = allEpsilons[0:len(A)]/system1.detB
    epsilonsY     = allEpsilons[len(A):]/system1.detB
    xVector       = 'Nash Equilibrium:\nx = ('
    for i in range(len(epsilonsX)):
        xVector += str(epsilonsX[i][0]) + ', '
    xVector = xVector[:-2] + '),'
    yVector = 'y = ('
    for i in range(len(epsilonsY)):
        yVector += str(epsilonsY[i][0]) + ', '
    yVector = yVector[:-2] + ').'
    output  = xVector + '\n' + yVector
    return output


def makeOutput(epsilonValues, system):
    vector = ''
    for i in range(0,len(epsilonValues)):
        if(i > 0):
            vector  = vector[:-3]
            vector += ', '
        if(not epsilonValues[i].any()):
            vector  = vector[:-2] 
            vector += ' + '
        else:
            for j in range(0,len(epsilonValues[0])):
                fraction = str(Fraction(int(epsilonValues[i][j]), int(system.detB)))
                if(j > 0 and Fraction(fraction) == 1):
                    fraction = ''
                if(epsilonValues[i][0] != 0):
                    if(j == 1):
                        vector  = vector[:-5] + ' + '
                if(epsilonValues[i][j] != 0):
                    if(np.sign(int(epsilonValues[i][j])/int(system.detB)) == -1):
                        vector  = vector[:-3] + ' - '
                        vector += fraction[1:] + sympy.pretty(epsilon) + str(j) + ' + '
                    else:
                        vector += fraction + sympy.pretty(epsilon) + str(j) + ' + '
    return vector[:-3]



#Tests
def runGamut(actionsp1, actionsp2, minpayoff, maxpayoff, nrOfGames):
    for i in range(0, nrOfGames):
        call(['java', '-jar', 'gamut.jar', '-g', 'RandomGame', '-players', '2', '-normalize', '-int_payoffs', '-int_mult', '1', '-min_payoff'
              , str(minpayoff), '-max_payoff', str(maxpayoff), '-f', 'Tests/BoS' + str(i) + '.game', '-actions', str(actionsp1), str(actionsp2)])


def readGamutFileReturnMatrices(fileName, actionsp1, actionsp2):
    payoffp1   = np.zeros((actionsp1, actionsp2))
    payoffp2   = np.zeros((actionsp1, actionsp2))
    fileObject = open(str(fileName), 'r')
    for line in fileObject:
        if('#' not in line):
            listOfBracketsInLine = [int(s) for s in re.findall(r'\b\d+\b', line)]
            payoffp1[int(listOfBracketsInLine[0])-1][int(listOfBracketsInLine[1])-1] = listOfBracketsInLine[2]
            payoffp2[int(listOfBracketsInLine[0])-1][int(listOfBracketsInLine[1])-1] = listOfBracketsInLine[3]
    return payoffp1, payoffp2


def run100GamesTakeAverage(actionsp1, actionsp2):
    runGamut(actionsp1, actionsp2, 0, 10, 100)
    for i in range(0, 100):
        print(i)
        name   = 'Tests/BoS' + str(i) + '.game'
        A, B   = readGamutFileReturnMatrices(name, actionsp1, actionsp2)
        cProfile.runctx('findEquilibrium(A,B)', {'findEquilibrium':findEquilibrium}, {'A': A, 'B': B}, 'out.txt')
        p 	   = pstats.Stats('out.txt')
        stream = open('statistic.txt', 'w');
        stats  = pstats.Stats('out.txt', stream=stream)
        stats.sort_stats('cumulative').print_stats(1)
        stream.close()
        result = 'Tests/' + str(actionsp1) + 'x' + str(actionsp2) + 'g' + '0to10entry' + '.txt'
        with open('statistic.txt', 'r') as file:
            for line in file:
                line = line.replace(' ', '')
                if('seconds' in line):
                    line = line.split('in',1)[1]
                    with open(result , 'a') as f1:
                        f1.writelines(line)
    average = 0.0
    with open(result, 'r') as f2:
        for line in f2:
            average += float(line.split('s',1)[0])
        average /= 100
    with open(result, 'a') as f3:
        averageString = 'Average = ' + str(average) 
        f3.writelines(averageString)


#run100GamesTakeAverage(24,24)

#print(LemkeMain(M0, b0.T))
runGamut(13,13,0,10,1)
#for i in range(0,3):
name  = 'Tests/BoS' + str(0) + '.game'
A, B  = readGamutFileReturnMatrices(name, 13, 13)
#A    *= -1
#B    *= -1
#print(findProperEquilibrium(A,B))
#A = np.array([[0,0,1,1], [0,2,0,2]])
#B = np.array([[0,0,-1,-1], [0,-2,0,-2]])
#for i in range(1,17):
   #print(len(Batcher(i)))
print(findProperEquilibrium(A,B))
#A = np.array([[3,0,10], [9,8,0], [6,7,2]])
#B = np.array([[10,9,7], [10,8,10], [4,7,6]])
#print(findEquilibrium(A,B))


# x   4  0:4
# y   4  4:8
# p'  3  8:11
# p'' 3  11:14
# g   2  14:16
# q'  3  16:19
# q'' 3  19:22
# h   2  22:24

