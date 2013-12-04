#!/usr/bin/python
from __future__ import division
import sys
import numpy as np
import matplotlib.pyplot as plt
import math
from routines.integral_input import __init_integrals__ 
import routines.scf as scf
import routines.ao2mo as ao2mo
import routines.mp2 as mp2
import routines.cistdhf as cistdhf
import routines.ccsd as ccsd
import routines.eomccsd as eomccsd
import routines.eommbpt2 as eommbpt2
import routines.eommbptp2 as eommbptp2
import routines.eommbptd as eommbptd

""" Edit below to perform the calculation desired
"""

do_DIIS      = True 
do_ao2mo     = True
do_mp2       = True 
do_cistdhf   = True
do_ccsd      = True 
do_eomccsd   = True
do_eommbpt2  = False
do_eommbptp2 = False
do_eommbptd  = False
printops     = True
convergence  = 1.0e-8
printnum     = 6			# number eigenvalues to print for TDHF/CIS

""" Here the main routine runs
"""

# ask for location of integral/input files

if len(sys.argv)==1:
    print "Enter folder containing input .dat files"
    sys.exit(1)

LOCATION = sys.argv[1]

print "\n\t*** Begin quantum chemistry on:   "+sys.argv[1]

# create one and two electron integral arrays, as well as determine
# number of basis functions (dim) and number of electrons (Nelec)

print "\n\t*** Started creation of one- and two- electron integrals"

ENUC,Nelec,dim,S,T,V,Hcore,twoe = __init_integrals__(LOCATION)

print "\t*** Finished creation of one- and two- electron integrals\n"
# do SCF iteration

print "\t*** Begin SCF iteration, convergence requested: ",convergence,"a.u."
EN,orbitalE,C,P,F = scf.scf_iteration(convergence,ENUC,Nelec,dim,S,Hcore,twoe,printops,do_DIIS)

print "Total E(RHF): ", EN+ENUC,"a.u."
print "\t*** End of SCF Iteration\n"

if do_ao2mo == True:
    print "\t*** Begin AO to MO <pq||rs> transformation"
    ints,fs = ao2mo.transform_ao2mo(dim,twoe,C,orbitalE)
    print "\t*** Finished AO to MO <pq||rs> tranformation\n"


if do_mp2 == True:
    # Begin MP2 calculation
    print "\t*** Begin MP2 Energy calculation"
    mp2_corr = mp2.mp2_energy(dim,Nelec,ints,orbitalE)
    print "MP2 corr: ",mp2_corr,"a.u."
    print "Total E(MP2): ",mp2_corr+EN+ENUC,"a.u."
    print "\t*** Finished MP2 Energy calculation"

if do_cistdhf == True:
    print "\n\t*** Begin CIS/TDHF calculation"
    cistdhf.cistdhf(Nelec,dim,fs,ints,printnum)
    print "\t*** End CIS/TDHF calculation"

if do_ccsd == True:
    print "\n\t*** Begin CCSD calculation"
    ECCSD,T1,T2 = ccsd.ccsd(Nelec,dim,fs,ints,convergence,printops)
    print "E(CCSD): {0:.8f}".format(ECCSD+ENUC+EN),"a.u."
    print "\t*** End CCSD calculation"
elif do_ccsd == False:
    T1 = np.zeros((dim*2,dim*2))
    T2 = np.zeros((dim*2,dim*2,dim*2,dim*2))
    for a in range(Nelec,dim*2):
        for b in range(Nelec,dim*2):
            for i in range(0,Nelec):
                for j in range(0,Nelec):
                    T2[a,b,i,j] += ints[i,j,a,b]/(fs[i,i] + fs[j,j] - fs[a,a] - fs[b,b])

if do_eomccsd == True:
    print "\n\t*** Begin EOM-CCSD calculation"
    eomccsd.eomccsd(Nelec,dim,fs,ints,T1,T2)
    print "\t*** End EOM-CCSD calculation"

if do_eommbpt2 == True:
    print "\n\t*** Begin EOM-MBPT2 calculation"
    eommbpt2.eommbpt2(Nelec,dim,fs,ints,T1,T2)
    print "\t*** End EOM-MBPT2 calculation"

if do_eommbptd == True:
    print "\n\t*** Begin EOM-MBPT(D) calculation"
    eommbptd.eommbptd(Nelec,dim,fs,ints,T1,T2)
    print "\t*** End EOM-MBPT(D) calculation"

if do_eommbptp2 == True:
    print "\n\t*** Begin EOM-MBPT(2) calculation"
    eommbptp2.eommbptp2(Nelec,dim,fs,ints,T1,T2)
    print "\t*** End EOM-MBPT(2) calculation"



