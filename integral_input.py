#!/usr/bin/python

from __future__ import division
from numpy import genfromtxt

# Return compund index given four indices
def eint(a,b,c,d):
    if a > b: ab = a*(a+1)/2 + b
    else: ab = b*(b+1)/2 + a
    if c > d: cd = c*(c+1)/2 + d
    else: cd = d*(d+1)/2 + c
    if ab > cd: abcd = ab*(ab+1)/2 + cd
    else: abcd = cd*(cd+1)/2 + ab
    return abcd

def __init_integrals__(Integral_Location):

    '''
        Reads integrals and places them in 
            array format for future use in 
            quantum chemical programs (e.g.
            for an SCF procedure)
            
            Returns arrays: ENUC, Nelec, dim,
                S, T, V, Hcore, twoe

    '''    

    ####################################
    #
    #       FORM CORE HAMILTONIAN 
    #
    ####################################

    # ENUC is nuclear repulsion, S is overlap matrix, T is kinetic energy matrix,
    # V is potential energy matrix, Nelec is number electrons, dim is number of 
    # basis functions

    ENUC  = genfromtxt('./'+ Integral_Location + '/enuc.dat', dtype=float)
    S     = genfromtxt('./'+ Integral_Location + '/s.dat', dtype=None)
    T     = genfromtxt('./'+ Integral_Location + '/t.dat', dtype=None)
    V     = genfromtxt('./'+ Integral_Location + '/v.dat', dtype=None)
    Nelec = genfromtxt('./'+ Integral_Location + '/nelec.dat', dtype=int)
    dim   = genfromtxt('./'+ Integral_Location + '/nbf.dat', dtype=int)

    Hcore = T - V

    #####################################
    #
    #       TWO ELECTRON INTEGRALS
    #
    #####################################

    # Like the core hamiltonian, we need to grab the integrals from the
    # separate file, and put into ERIraw (ERI = electron repulsion integrals).
    # I chose to store the two electron integrals in a python dictionary. 
    # The function 'eint' generates a unique compund index for the unique two 
    # electron integral, and maps this index to the corresponding integral value.
    # 'twoe' is the name of the dictionary containing these.

    ERIraw = genfromtxt('./'+ Integral_Location + '/eri.dat',dtype=None)
    twoe = {eint(row[0],row[1],row[2],row[3]) : row[4] for row in ERIraw}

    return ENUC,Nelec,dim,S,T,V,Hcore,twoe

