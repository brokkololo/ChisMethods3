from numpy import *
#import matplotlib.pyplot as plt
from scipy.linalg import norm

def initFun(tInit, x):
    """This function calculates initial function and is used to fulfill
    the first (presely, zero) time layer in the grid solution aray
    """
    res = sin(x)
    return res

def boundaryLeftFun(t, xLeft):
    """This function calculates left boundary function
    """
    res = sin(xLeft) * cos(t)
    return 0

def boundaryRightFun(t, xRight):
    """This function calculates left boundary function
    """
    res = sin(xRight) * cos(t)
    return 0

def nonLinFun(u):
    """This function calculates nonlinear part of difference operator
    """
    return exp(u)

def nonLinFunDer(u):
    """This function calculates the derivative with respect to time for
    nonlinear part of difference operator
    """
    return exp(u)

def heterFun(t, x):
    """This function calculates the heterogenious function
    """
    res = -exp(sin(x)*cos(t))*sin(x)*sin(t)+sin(x)*cos(t)
    return res

def exactSolution(t, x):
    """This function calculates the exact solution and
    is used to estimate the error
    """
    res = sin(x)*cos(t)
    return res

xLeft = 0            #xLeft - left boundary
xRight = 3.14159265  #xRight - right boundary, we solve the boundary problen on the segment [xLeft, xRight]
tInit = 0            #tInit - initial time
tFinal = 8.0         #tFinal - final time, we solve the boundary problen on the segment [tInint, tFinal]
#numberOfSegm_x = 31  #numSeg - number of subsegments in segment [x_left, x_right]
#numberOfSegm_t = 80  #numSeg - number of subsegments in segment [tInint, tFinal]
errorNorm = 0
epsNewton = 10**(-5)

#Below we calculate the time step tStep,
#construct the uniform time grid and allocate memory for solution-array
#Note, that number of dots in time grid = numSeg + 1

piVal = 3.14159265
xStep = piVal / 50 
tStep = 0.0001

numberOfSegm_x = int((xRight - xLeft) / xStep)  
numberOfSegm_t = int((tFinal - tInit) / tStep)  

spaceGrid = linspace(xLeft, xRight, numberOfSegm_x + 1)
timeGrid = linspace(tInit, tFinal, numberOfSegm_t + 1)

uSolution = zeros((numberOfSegm_t + 1, numberOfSegm_x + 1))
uError = zeros((numberOfSegm_t + 1, numberOfSegm_x + 1))
alpha = zeros(numberOfSegm_x + 1)
beta = zeros(numberOfSegm_x + 1)
Ai = tStep / (2*(xStep ** 2))
Bi = tStep / (2*(xStep ** 2))
yk = zeros(numberOfSegm_x + 1)
yk1= zeros(numberOfSegm_x + 1)

for i in range(numberOfSegm_x + 1):
    uSolution[0][i] = initFun(tInit, spaceGrid[i])


for j in range(numberOfSegm_t):
    uSolution[j+1][0] = boundaryLeftFun(timeGrid[j+1], xLeft)
    uSolution[j+1][numberOfSegm_x] = boundaryRightFun(timeGrid[j+1], xRight)
    flag = True
    yk = copy(uSolution[j])
    alpha[0] = 0
    beta[0] = boundaryLeftFun(timeGrid[j+1], xLeft)
    yk1[numberOfSegm_x] = boundaryRightFun(timeGrid[j+1], xRight)
   
    while flag:
        
        for i in range(0, numberOfSegm_x ):
            Ficorrect = - (tStep/2) * ((uSolution[j][i-1] - 2*uSolution[j][i] + uSolution[j][i+1]) / xStep**2)
            Fi = -(nonLinFun(yk[i]) - nonLinFun(uSolution[j][i]) - nonLinFunDer(yk[i])*yk[i] - tStep* heterFun(timeGrid[j], spaceGrid[i]) + Ficorrect)
            Ci = nonLinFunDer(yk[i]) + tStep/(xStep**2)
            alpha[i+1] = Bi/(Ci - alpha[i]*Ai)
            beta[i+1] = (Ai*beta[i] + Fi)/(Ci - alpha[i]*Ai)
        for i in range(numberOfSegm_x -1, 0, -1 ):
            yk1[i] = alpha[i+1]*yk1[i+1] + beta[i+1]
        flag = (norm(yk1 - yk) > epsNewton)
        yk = copy(yk1)
    
    uSolution[j+1] = copy(yk1)
    for i in range(0, numberOfSegm_x + 1):
        uError[j+1][i] = abs(uSolution[j+1][i] - exactSolution(timeGrid[j+1], spaceGrid[i]))
        if (uError[j+1][i] > errorNorm):
            errorNorm = uError[j+1][i]

#print(uError)
#print(uSolution)
print(errorNorm)
#plt.plot(timeGrid, uError[:,16])
#plt.plot(timeGrid, uSolution[:,1])

#plt.show()