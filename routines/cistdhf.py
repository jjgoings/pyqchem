#!/usr/bin/python

from __future__ import division
import math
import numpy as np

def cistdhf(Nelec,dim,fs,ints,printnum):
    ######################################################
    #
    #   CIS & TDHF CALCULATION 
    #
    ######################################################
    np.set_printoptions(precision=4)
    NOV = Nelec*(2*dim - Nelec) # number occupied * number virtual
    if printnum >= NOV:
        printnum = NOV
    A = np.zeros((NOV,NOV))
    B = np.zeros((NOV,NOV))
    I = -1
    for i in range(0,Nelec):
      for a in range(Nelec,dim*2):
	I = I + 1
	J = -1
	for j in range(0,Nelec):
	  for b in range(Nelec,dim*2):
	    J = J+1
	    A[I,J] = (fs[a,a] - fs[i,i]) * ( i == j ) * (a == b) + ints[a,j,i,b]
	    B[I,J] =  ints[a,b,i,j]

    #print B
    # Solve CIS matrix equation
    ECIS,CCIS = np.linalg.eig(A)
    ECIS = np.real(np.sort(ECIS))*27.2114
    print "E(CIS) = ",ECIS[:printnum],"eV"
    # Solve TDHF matrix equation
    #M = np.bmat([[A,B],[-B,-A]])
    M = np.dot((A+B),(A-B))
    ETD,CTD = np.linalg.eig(M)
    ETD = np.real(np.sort(np.sqrt(ETD)))*27.2114
    print "E(TDHF) = ",ETD[:printnum],"eV"
    del ETD
    del CTD
    del ECIS
    del CCIS
    del A
    del B
    del M

