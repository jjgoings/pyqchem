from __future__ import division
import sys
#import re
import math
import time
import matplotlib.pyplot as plt
import numpy as np

#######################################################
#
#  EOM-MBPT(2) CALCULATION
#
#######################################################

def eommbptp2(Nelec,dim,fs,ints,ts,td):
    dim = dim*2
    NOV = Nelec*(dim-Nelec)
    print "\t\t** Constructing singles-singles block" 
    t0 = time.time()
    HSS = np.zeros((NOV,NOV))
    ia = -1
    for i in range(0,Nelec):
      for a in range(Nelec,dim):
        ia = ia + 1
        kc = -1
        for k in range(0,Nelec):
          for c in range(Nelec,dim):
            kc = kc + 1
            HSS[ia,kc] += fs[a,c]*(i==k) # 1
            HSS[ia,kc] += -fs[i,k]*(a==c) # 2 
            HSS[ia,kc] += ints[a,k,i,c] # 3
            for m in range(0,Nelec):
              for e in range(Nelec,dim):
                HSS[ia,kc] += ints[k,m,c,e]*td[e,a,m,i] # 12
                for n in range(0,Nelec):
                  HSS[ia,kc] += -0.5*(i==k)*ints[m,n,c,e]*td[a,e,m,n] # 10
                for f in range(Nelec,dim):
                  HSS[ia,kc] += -0.5*(a==c)*ints[k,m,e,f]*td[e,f,i,m] # 11
    t1 = time.time()
    print "\t\t** Finished singles-singles block. Time: ", t1-t0 


    print "\t\t** Constructing singles-doubles block" 
    t2 = time.time()
    HSD = np.zeros((NOV,NOV*NOV))
    ia = -1
    for i in range(0,Nelec):
      for a in range(Nelec,dim):
        ia = ia + 1
        kcld = -1
        for k in range(0,Nelec):
          for c in range(Nelec,dim):
            for l in range(0,Nelec):
              for d in range(Nelec,dim):
                kcld = kcld + 1
                HSD[ia,kcld] += 0.5*ints[a,l,c,d]*(i==k) # 17
                HSD[ia,kcld] += -0.5*ints[k,l,i,d]*(a==c) # 18
    t3 = time.time()
    print "\t\t** Finished singles-doubles block. Time: ", t3-t2 

    print "\t\t** Constructing doubles-singles block" 
    t4 = time.time()
    HDS = np.zeros((NOV*NOV,NOV))
    iajb = -1
    for i in range(0,Nelec):
      for a in range(Nelec,dim):
        for j in range(0,Nelec):
          for b in range(Nelec,dim):
            iajb += 1
            kc = -1
            for k in range(0,Nelec):
              for c in range(Nelec,dim):
                kc += 1
                HDS[iajb,kc] += (i==k)*ints[a,b,c,j] - \
                                (j==k)*ints[a,b,c,i] # 22
                HDS[iajb,kc] += (b==c)*ints[k,a,i,j] - \
                                (a==c)*ints[k,b,i,j] # 23
    t5 = time.time()
    print "\t\t** Finished doubles-singles block. Time: ",t5-t4 

    print "\t\t** Constructing doubles-doubles block" 
    t6 = time.time()
    HDD = np.zeros((NOV*NOV,NOV*NOV))
    iajb = -1
    for i in range(0,Nelec):
      for a in range(Nelec,dim):
        for j in range(0,Nelec):
          for b in range(Nelec,dim):
            iajb += 1
            kcld = -1
            for k in range(0,Nelec):
              for c in range(Nelec,dim):
                for l in range(0,Nelec):
                  for d in range(Nelec,dim):
                    kcld += 1
                    HDD[iajb,kcld] += (j==k)*(i==l)*(a==d)*fs[b,c] - \
                                      (j==k)*(i==l)*(b==d)*fs[a,c] # 52
                    HDD[iajb,kcld] += (j==l)*(a==d)*(b==c)*fs[k,i] - \
                                      (i==l)*(a==d)*(b==c)*fs[k,j] # 53
                    HDD[iajb,kcld] += 0.5*(i==k)*(j==l)*ints[a,b,c,d]  # 54 
                    HDD[iajb,kcld] += 0.5*(a==c)*(b==d)*ints[k,l,i,j] # 55
                    HDD[iajb,kcld] += (i==l)*(a==d)*ints[k,b,c,j] - \
                                      (j==l)*(a==d)*ints[k,b,c,i] - \
                                      (i==l)*(b==d)*ints[k,a,c,j] + \
                                      (j==l)*(b==d)*ints[k,a,c,i] # 56
    t7 = time.time()
    print "\t\t** Finished doubles-doubles block. Time: ",t7-t6 

    eomMatrix = np.bmat([[HSS,HSD],[HDS,HDD]])

    print "\t\t** Begin full diagonalization"
    print "\t\t * Matrix dinension:  ", str(len(eomMatrix)) + "x" + str(len(eomMatrix)) 
    eomEVal,eomEVec = np.linalg.eig(eomMatrix)

    print "\nExcitations (eV):"
    for excitation in (np.sort(np.real(eomEVal))*27.21138386):
      if (excitation > 1.0) and (excitation < 50.0):
        print (excitation)
        continue
    plt.imshow(eomMatrix,interpolation='nearest',cmap='jet',alpha=0.75)
    plt.colorbar()
    plt.savefig('eom-mbptP2.png',bbox_inches='tight')
    plt.close()
    del eomMatrix, eomEVal, eomEVec
