from __future__ import division
import sys
#import re
import math
import time
import matplotlib.pyplot as plt
import numpy as np

#######################################################
#
#  EOM-MBPT2 CALCULATION
#
#######################################################

def eommbpt2(Nelec,dim,fs,ints,ts,td):
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
                for e in range(Nelec,dim):
                  HDS[iajb,kc] += ints[k,a,c,e]*td[e,b,i,j] - \
                                  ints[k,b,c,e]*td[e,a,i,j] # 34 
                  for f in range(Nelec,dim):
                    HDS[iajb,kc] += 0.5*(b==c)*ints[k,a,e,f]*td[e,f,i,j] - \
                                    0.5*(a==c)*ints[k,b,e,f]*td[e,f,i,j]  # 32
                for m in range(0,Nelec):
                  HDS[iajb,kc] += ints[k,m,c,j]*td[a,b,m,i] - \
                                  ints[k,m,c,i]*td[a,b,m,j] #35
                  for n in range(0,Nelec):
                    HDS[iajb,kc] += 0.5*(i==k)*ints[m,n,c,j]*td[a,b,m,n] - \
                                    0.5*(j==k)*ints[m,n,c,i]*td[a,b,m,n]  # 33 
                  for e in range(Nelec,dim):
                    HDS[iajb,kc] += (i==k)*ints[a,m,c,e]*td[e,b,m,j] - \
                                    (i==k)*ints[b,m,c,e]*td[e,a,m,j] - \
                                    (j==k)*ints[a,m,c,e]*td[e,b,m,i] + \
                                    (j==k)*ints[b,m,c,e]*td[e,a,m,i] # 30
                    HDS[iajb,kc] += (b==c)*ints[k,m,i,e]*td[e,a,m,j] - \
                                    (b==c)*ints[k,m,j,e]*td[e,a,m,i] - \
                                    (a==c)*ints[k,m,i,e]*td[e,b,m,j] + \
                                    (a==c)*ints[k,m,j,e]*td[e,b,m,i]  # 31
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
                    for e in range(Nelec,dim):
                      HDD[iajb,kcld] += 0.5*(b==d)*ints[k,l,c,e]*td[e,a,i,j] - \
                                        0.5*(a==d)*ints[k,l,c,e]*td[e,b,i,j] # 70 
                      for f in range(Nelec,dim):
                        HDD[iajb,kcld] += 0.25*(a==c)*(b==d)*ints[k,l,e,f]*td[e,f,i,j] # 66
                    for m in range(0,Nelec):
                      HDD[iajb,kcld] += 0.5*(j==l)*ints[k,m,c,d]*td[a,b,m,i] - \
                                        0.5*(i==l)*ints[k,m,c,d]*td[a,b,m,j] # 68 
                      for n in range(0,Nelec):
                        HDD[iajb,kcld] += 0.25*(i==k)*(j==l)*ints[m,n,c,d]*td[a,b,m,n] # 65
                      for e in range(Nelec,dim):
                        HDD[iajb,kcld] += (j==l)*(b==d)*ints[k,m,c,e]*td[e,a,m,i] - \
                                          (i==l)*(b==d)*ints[k,m,c,e]*td[e,a,m,j] - \
                                          (j==l)*(a==d)*ints[k,m,c,e]*td[e,b,m,i] + \
                                          (i==l)*(a==d)*ints[k,m,c,e]*td[e,b,m,j]  # 67 
                        for n in range(0,Nelec):
                          HDD[iajb,kcld] += 0.5*(i==k)*(j==l)*(a==d)*ints[m,n,e,c]*td[b,e,n,m] - \
                                            0.5*(i==k)*(j==l)*(b==d)*ints[m,n,e,c]*td[a,e,n,m] # 71
                        for f in range(Nelec,dim):
                          HDD[iajb,kcld] += 0.5*(a==c)*(i==l)*(b==d)*ints[m,k,e,f]*td[f,e,j,m] - \
                                            0.5*(a==c)*(j==l)*(b==d)*ints[m,k,e,f]*td[f,e,i,m] # 69
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
    plt.savefig('eom-mbpt2.png',bbox_inches='tight')
    plt.close()
    del eomMatrix, eomEVal, eomEVec

