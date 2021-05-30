from __future__ import division
import sys
import math
import numpy as np

#Part 1: Read Matrix input files and convert too matrixes.

Alpha = sys.argv[1]
Beta = sys.argv[2]
Overlap = sys.argv[3]

Alpha_mat = []
Beta_mat = []
Overlap_mat = []

#Defined function for reading in input file and outputting matrix
#Note: The input mat file will not work if terminated with ' '
def getmat(inputmat):
    mat = open(inputmat, 'r')
    outputmat = [[float(num) for num in line.split(' ')] for line in mat]
    return outputmat
    mat.close()

#Re-write matrices with function defined above.
Alpha_mat= getmat(Alpha)
print(Alpha_mat)
Beta_mat= getmat(Beta)
print(Beta_mat)
Overlap_mat= getmat(Overlap)
print(Overlap_mat)

#Part 2: Calculate Trace of PalphaS and PbetaS: 
