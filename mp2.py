#!/usr/bin/python

from __future__ import division
import math
import numpy as np

def mp2_energy(dim,Nelec,ints,E):
    #####################################################
    #
    #  MP2 ENERGY CALCULATION
    #
    #####################################################

    CC = 0.0
    for i in range(0,Nelec):
	for j in range(0,Nelec):
	    for a in range(Nelec,dim*2):
		for b in range(Nelec,dim*2):
		    CC += 0.25*(ints[i,j,a,b]*ints[i,j,a,b])/(E[i//2] + E[j//2] - E[a//2] - E[b//2])

    return CC
