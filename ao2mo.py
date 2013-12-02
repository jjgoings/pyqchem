#!/usr/bin/python

from __future__ import division
import math
import numpy as np

####################################
#
#   FUNCTIONS
#
####################################

# Return compund index given four indices
def eint(a,b,c,d):
    if a > b: ab = a*(a+1)/2 + b
    else: ab = b*(b+1)/2 + a
    if c > d: cd = c*(c+1)/2 + d
    else: cd = d*(d+1)/2 + c
    if ab > cd: abcd = ab*(ab+1)/2 + cd
    else: abcd = cd*(cd+1)/2 + ab
    return abcd

# Return Value of two electron integral
# Example: (12|34) = tei(1,2,3,4)
# I use chemists notation for the SCF procedure. 
def tei(a,b,c,d,twoe):
    return twoe.get(eint(a,b,c,d),0.0)

def transform_ao2mo(dim,twoe,C,E):
    ####################################################
    #
    #   2-E INTEGRAL TRANSFORM TO SPATIAL MO BASIS
    #
    ####################################################

    # put AO ERI into 4-D array for easy indexing
    INT = np.zeros((dim,dim,dim,dim))
    for i in range(0,dim):
	for j in range(0,dim):
	    for k in range(0,dim):
		for l in range(0,dim):
		    INT[i,j,k,l] = tei(i+1,j+1,k+1,l+1,twoe)
    
    del twoe
    
    # now transform integrals smartly (N^5)
    MO2 = np.zeros((dim,dim,dim,dim))
    temp = np.zeros((dim,dim,dim,dim))
    temp2 = np.zeros((dim,dim,dim,dim))
    temp3= np.zeros((dim,dim,dim,dim))
    for mu in range(0,dim):
	for i in range(0,dim):
	    temp[mu,:,:,:] += C[i,mu]*INT[i,:,:,:]
	for nu in range(0,dim):
	    for j in range(0,dim):
		temp2[mu,nu,:,:] += C[j,nu]*temp[mu,j,:,:]
	    for lam in range(0,dim):
		for k in range(0,dim):
		    temp3[mu,nu,lam,:] += C[k,lam]*temp2[mu,nu,k,:]
		for sig in range(0,dim):
		    for l in range(0,dim):
			MO2[mu,nu,lam,sig] += C[l,sig]*temp3[mu,nu,lam,l]
    
    del temp
    del temp2
    del temp3
    
    ####################################################
    #
    #  CONVERT SPATIAL TO SPIN ORBITAL MO
    #
    ####################################################

    # This makes the spin basis double bar integral (physicists' notation)

    ints=np.zeros((dim*2,dim*2,dim*2,dim*2))
    for p in range(0,dim*2):
	for q in range(0,dim*2):
	    for r in range(0,dim*2):
		for s in range(0,dim*2):
		    value1 = MO2[(p)//2,(r)//2,(q)//2,(s)//2] * (p%2 == r%2) * (q%2 == s%2)
		    value2 = MO2[(p)//2,(s)//2,(q)//2,(r)//2] * (p%2 == s%2) * (q%2 == r%2)
		    ints[p,q,r,s] = value1 - value2
    
    del MO2
   
    #####################################################
    #
    #  Spin basis fock matrix eigenvalues 
    #
    #####################################################

    fs = np.zeros((dim*2))
    for i in range(0,dim*2):
	fs[i] = E[i//2]

    fs = np.diag(fs)
    
    return ints, fs
