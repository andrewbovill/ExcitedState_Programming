from __future__ import division
import sys
import math
import numpy as np
from numpy import linalg as la

#Part 1: Read Matrix input files and convert too matrixes.

Alpha = sys.argv[1]
Beta = sys.argv[2]
S = sys.argv[3]

Alpha_mat = []
Beta_mat = []
S_mat = []

#Defined function for reading in input file and outputting matrix
#Note: The input mat file will not work if terminated with ' '
def getmat(inputmat):
    mat = open(inputmat, 'r')
    outputmat = [[float(num) for num in line.split(' ')] for line in mat]
    return outputmat 
    mat.close()

#Re-write matrices with function defined above.
Alpha_mat= getmat(Alpha)
Alpha_mat= np.array(Alpha_mat)

print(Alpha_mat)

Beta_mat= getmat(Beta)
Beta_mat= np.array(Beta_mat)

print(Beta_mat)
S_mat= getmat(S)
S_mat= np.array(S_mat)
print(S_mat)

#Part 2: Calculate Trace of PalphaS and PbetaS: 
PalphaS = np.dot(Alpha_mat,S_mat)
print(PalphaS)
TracePalphaS = np.trace(PalphaS)
print(TracePalphaS)

PbetaS = np.dot(Beta_mat,S_mat)
print(PbetaS)
TracePbetaS = np.trace(PbetaS)
print(TracePbetaS)

#Part #: Starting from S, obtain matrix S ^(1/2) 

c = la.eig(S_mat)
