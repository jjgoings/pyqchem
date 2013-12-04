#!/usr/bin/python

####################################
#
# QUANTUM CHEMISTRY ON H2O 
#
####################################

from __future__ import division
import sys
import math
import numpy as np


def ccsd(Nelec,dim,fs,ints,convergence,printops):
    #######################################################
    #
    #   CCSD CALCULATION
    #
    #######################################################
    dim = dim*2

    ts = np.zeros((dim,dim))
    td = np.zeros((dim,dim,dim,dim))

    # Initial guess T2

    for a in range(Nelec,dim):
      for b in range(Nelec,dim):
	for i in range(0,Nelec):
	  for j in range(0,Nelec):
	    td[a,b,i,j] += ints[i,j,a,b]/(fs[i,i] + fs[j,j] - fs[a,a] - fs[b,b])

    # Make denominator arrays
    Dai = np.zeros((dim,dim))
    for a in range(Nelec,dim):
      for i in range(0,Nelec):
	Dai[a,i] = fs[i,i] - fs[a,a]

    Dabij = np.zeros((dim,dim,dim,dim))
    for a in range(Nelec,dim):
      for b in range(Nelec,dim):
	for i in range(0,Nelec):
	  for j in range(0,Nelec):
	    Dabij[a,b,i,j] = fs[i,i] + fs[j,j] - fs[a,a] - fs[b,b]

    def taus(a,b,i,j):
      taus = td[a,b,i,j] + 0.5*(ts[a,i]*ts[b,j] - ts[b,i]*ts[a,j])
      return taus

    def tau(a,b,i,j):
      tau = td[a,b,i,j] + ts[a,i]*ts[b,j] - ts[b,i]*ts[a,j]
      return tau

    def updateintermediates(x):
      if x == True:
	Fae = np.zeros((dim,dim))
	for a in range(Nelec,dim):
	  for e in range(Nelec,dim):
	    Fae[a,e] = (1 - (a == e))*fs[a,e]
	    for m in range(0,Nelec):
	      Fae[a,e] += -0.5*fs[m,e]*ts[a,m]
	      for f in range(Nelec,dim):
		Fae[a,e] += ts[f,m]*ints[m,a,f,e] 
		for n in range(0,Nelec):
		  Fae[a,e] += -0.5*taus(a,f,m,n)*ints[m,n,e,f]

	Fmi = np.zeros((dim,dim))
	for m in range(0,Nelec):
	  for i in range(0,Nelec):
	    Fmi[m,i] = (1 - (m == i))*fs[m,i]
	    for e in range(Nelec,dim):
	      Fmi[m,i] += 0.5*ts[e,i]*fs[m,e]
	      for n in range(0,Nelec):
		Fmi[m,i] += ts[e,n]*ints[m,n,i,e] 
		for f in range(Nelec,dim):
		  Fmi[m,i] += 0.5*taus(e,f,i,n)*ints[m,n,e,f]

	Fme = np.zeros((dim,dim))
	for m in range(0,Nelec):
	  for e in range(Nelec,dim):
	    Fme[m,e] = fs[m,e]
	    for n in range(0,Nelec):
	      for f in range(Nelec,dim):
		Fme[m,e] += ts[f,n]*ints[m,n,e,f]

	Wmnij = np.zeros((dim,dim,dim,dim))
	for m in range(0,Nelec):
	  for n in range(0,Nelec):
	    for i in range(0,Nelec):
	      for j in range(0,Nelec):
		Wmnij[m,n,i,j] = ints[m,n,i,j]
		for e in range(Nelec,dim):
		  Wmnij[m,n,i,j] += ts[e,j]*ints[m,n,i,e] - ts[e,i]*ints[m,n,j,e]
		  for f in range(Nelec,dim):
		    Wmnij[m,n,i,j] += 0.25*tau(e,f,i,j)*ints[m,n,e,f]

	Wabef = np.zeros((dim,dim,dim,dim))
	for a in range(Nelec,dim):
	  for b in range(Nelec,dim):
	    for e in range(Nelec,dim):
	      for f in range(Nelec,dim):
		Wabef[a,b,e,f] = ints[a,b,e,f]
		for m in range(0,Nelec):
		  Wabef[a,b,e,f] += -ts[b,m]*ints[a,m,e,f] + ts[a,m]*ints[b,m,e,f]
		  for n in range(0,Nelec):
		    Wabef[a,b,e,f] += 0.25*tau(a,b,m,n)*ints[m,n,e,f]

	Wmbej = np.zeros((dim,dim,dim,dim))
	for m in range(0,Nelec):
	  for b in range(Nelec,dim):
	    for e in range(Nelec,dim):
	      for j in range(0,Nelec):
		Wmbej[m,b,e,j] = ints[m,b,e,j]
		for f in range(Nelec,dim):
		  Wmbej[m,b,e,j] += ts[f,j]*ints[m,b,e,f]
		for n in range(0,Nelec):
		  Wmbej[m,b,e,j] += -ts[b,n]*ints[m,n,e,j]
		  for f in range(Nelec,dim):
		    Wmbej[m,b,e,j] += -(0.5*td[f,b,j,n] + ts[f,j]*ts[b,n])*ints[m,n,e,f]

	return Fae, Fmi, Fme, Wmnij, Wabef, Wmbej

    def makeT1(x,ts,td):
      if x == True:
	tsnew = np.zeros((dim,dim))
	for a in range(Nelec,dim):
	  for i in range(0,Nelec):
	    tsnew[a,i] = fs[i,a]
	    for e in range(Nelec,dim):
	      tsnew[a,i] += ts[e,i]*Fae[a,e]
	    for m in range(0,Nelec):
	      tsnew[a,i] += -ts[a,m]*Fmi[m,i]
	      for e in range(Nelec,dim):
		tsnew[a,i] += td[a,e,i,m]*Fme[m,e]
		for f in range(Nelec,dim):
		  tsnew[a,i] += -0.5*td[e,f,i,m]*ints[m,a,e,f]
		for n in range(0,Nelec):
		  tsnew[a,i] += -0.5*td[a,e,m,n]*ints[n,m,e,i]
	    for n in range(0,Nelec):
	      for f in range(Nelec,dim): 
		tsnew[a,i] += -ts[f,n]*ints[n,a,i,f]
	    tsnew[a,i] = tsnew[a,i]/Dai[a,i]
      return tsnew

    def makeT2(x,ts,td):
      if x == True:
	tdnew = np.zeros((dim,dim,dim,dim))
	for a in range(Nelec,dim):
	  for b in range(Nelec,dim):
	    for i in range(0,Nelec):
	      for j in range(0,Nelec):
		tdnew[a,b,i,j] += ints[i,j,a,b]
		for e in range(Nelec,dim):
		  tdnew[a,b,i,j] += td[a,e,i,j]*Fae[b,e] - td[b,e,i,j]*Fae[a,e]
		  for m in range(0,Nelec):
		    tdnew[a,b,i,j] += -0.5*td[a,e,i,j]*ts[b,m]*Fme[m,e] + 0.5*td[a,e,i,j]*ts[a,m]*Fme[m,e]
		    continue
		for m in range(0,Nelec):
		  tdnew[a,b,i,j] += -td[a,b,i,m]*Fmi[m,j] + td[a,b,j,m]*Fmi[m,i]
		  for e in range(Nelec,dim):
		    tdnew[a,b,i,j] += -0.5*td[a,b,i,m]*ts[e,j]*Fme[m,e] + 0.5*td[a,b,i,m]*ts[e,i]*Fme[m,e]
		    continue
		for e in range(Nelec,dim):
		  tdnew[a,b,i,j] += ts[e,i]*ints[a,b,e,j] - ts[e,j]*ints[a,b,e,i]
		  for f in range(Nelec,dim):
		    tdnew[a,b,i,j] += 0.5*tau(e,f,i,j)*Wabef[a,b,e,f]
		    continue
		for m in range(0,Nelec):
		  tdnew[a,b,i,j] += -ts[a,m]*ints[m,b,i,j] + ts[b,m]*ints[m,a,i,j]  
		  for e in range(Nelec,dim):
		    tdnew[a,b,i,j] += td[a,e,i,m]*Wmbej[m,b,e,j] - ts[e,i]*ts[a,m]*ints[m,b,e,j]
		    tdnew[a,b,i,j] += -td[a,e,j,m]*Wmbej[m,b,e,i] + ts[e,j]*ts[a,m]*ints[m,b,e,i]
		    tdnew[a,b,i,j] += -td[b,e,i,m]*Wmbej[m,a,e,j] - ts[e,i]*ts[b,m]*ints[m,a,e,j]
		    tdnew[a,b,i,j] += td[b,e,j,m]*Wmbej[m,a,e,i] - ts[e,j]*ts[b,m]*ints[m,a,e,i]
		    continue
		  for n in range(0,Nelec):
		    tdnew[a,b,i,j] += 0.5*tau(a,b,m,n)*Wmnij[m,n,i,j]
		    continue
		tdnew[a,b,i,j] = tdnew[a,b,i,j]/Dabij[a,b,i,j] 
	return tdnew

    def ccsdenergy():
      ECCSD = 0.0
      for i in range(0,Nelec):
	for a in range(Nelec,dim):
	  ECCSD += fs[i,a]*ts[a,i]
	  for j in range(0,Nelec):
	    for b in range(Nelec,dim):
	      ECCSD += 0.25*ints[i,j,a,b]*td[a,b,i,j] + 0.5*ints[i,j,a,b]*(ts[a,i])*(ts[b,j]) 
      return ECCSD

    # CCSD iteration
    do_DIIS = True 
    ECCSD = 0.0
    num_e = 4
    Error1Set = []
    Error2Set = []
    T1Set     = []
    T2Set     = []
    for j in xrange(0,60):
        OLDCC = ECCSD
        Fae,Fmi,Fme,Wmnij,Wabef,Wmbej = updateintermediates(True)
        ts = makeT1(True,ts,td)
        td = makeT2(True,ts,td)
        
        if do_DIIS == True:
            td = td.reshape((dim*dim,dim*dim))
            if j == 1 or 2:
                T1Set.append(ts)
                T2Set.append(td)
            elif j > 2:
                if len(Error1Set) < num_e:
                    T1Set.append(ts)
                    T2Set.append(td)
                    Error1Set.append((T1Set[j-1] - T1Set[j-1]))
                    Error2Set.append((T2Set[j] - T2Set[j-1]))
                elif len(Error1Set) >= num_e:
                    del T1Set[0]
                    del T2Set[0]
                    del Error1Set[0]
                    del Error2Set[0]
                    T1Set.append(ts)
                    T2Set.append(td)
                    Error1Set.append(T1Set[-1] - T1Set[-2])
                    Error2Set.append(T2Set[-1] - T2Set[-2])
            NErr = len(Error1Set)
            if NErr > 2:
                Bmat1 =  Bmat2 = np.zeros((NErr+1,NErr+1))
                ZeroVec = np.zeros((NErr+1))
                ZeroVec[-1] = -1.0
                for a in range(0,NErr):
                    for b in range(0,a+1):
                        Bmat1[a,b] = Bmat1[b,a] = np.trace(np.dot(Error1Set[a].T,Error1Set[b]))
                        Bmat2[a,b] = Bmat2[b,a] = np.trace(np.dot(Error2Set[a].T,Error2Set[b]))
                        Bmat1[a,NErr] = Bmat1[NErr,a] = -1.0
                        Bmat2[a,NErr] = Bmat2[NErr,a] = -1.0
                try:
                    coeff1 = np.linalg.solve(Bmat1,ZeroVec)
                    coeff2 = np.linalg.solve(Bmat2,ZeroVec)
                except np.linalg.linalg.LinAlgError as err:
                    if 'Singular matrix' in err.message:
                        print '\tSingular B matrix, turning off DIIS'
                        do_DIIS = False
                else:
                    T1new = 0.0
                    T2new = 0.0
                    for i in range(0,len(coeff1)-1):
                        T1new += coeff1[i]*T1Set[i]
                        T2new += coeff2[i]*T2Set[i]
                    ts = T1new
                    td = T2new
            td = td.reshape((dim,dim,dim,dim))
        ECCSD = ccsdenergy()
        DECC = abs(ECCSD - OLDCC)
        if DECC < convergence:
            print "TOTAL ITERATIONS: ",j
            break
        if printops == True:
            print "E corr: {0:.12f}".format(ECCSD),"a.u.",'\t',"DeltaE: {0:.12f}".format(DECC)
    return ECCSD,ts,td 
        
