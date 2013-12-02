
from __future__ import division
import sys
#import re
import math
import time
import matplotlib.pyplot as plt
import numpy as np

#######################################################
#
#  EOM-CCSD CALCULATION
#
#######################################################

def eomccsd(Nelec,dim,fs,ints,ts,td):
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
            for e in range(Nelec,dim):
              HSS[ia,kc] += ints[a,k,e,c]*ts[e,i] #6
              HSS[ia,kc] += -(a==c)*fs[k,e]*ts[e,i] #8
            for m in range(0,Nelec):
              HSS[ia,kc] += -ints[m,k,i,c]*ts[a,m] # 7
              HSS[ia,kc] += -(i==k)*fs[m,c]*ts[a,m] # 9
              for e in range(Nelec,dim):
                HSS[ia,kc] += ints[a,m,c,e]*ts[e,m]*(i==k) # 4
                HSS[ia,kc] += -ints[k,m,i,e]*ts[e,m]*(a==c) # 5
                HSS[ia,kc] += ints[k,m,c,e]*td[e,a,m,i] # 12
                HSS[ia,kc] += -ints[m,k,e,c]*ts[e,i]*ts[a,m] # 13
                for n in range(0,Nelec):
                  HSS[ia,kc] += -0.5*(i==k)*ints[m,n,c,e]*td[a,e,m,n] # 10
                for f in range(Nelec,dim):
                  HSS[ia,kc] += -0.5*(a==c)*ints[k,m,e,f]*td[e,f,i,m] # 11
                  HSS[ia,kc] += -(a==c)*ints[k,m,e,f]*ts[e,i]*ts[f,m] # 14
              for n in range(0,Nelec):
                for f in range(Nelec,dim):
                  HSS[ia,kc] += -(i==k)*ints[m,n,c,f]*ts[a,m]*ts[f,n] # 15
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
                HSD[ia,kcld] += (i==k)*(a==c)*fs[l,d] # 16
                HSD[ia,kcld] += 0.5*ints[a,l,c,d]*(i==k) # 17
                HSD[ia,kcld] += -0.5*ints[k,l,i,d]*(a==c) # 18
                for e in range(Nelec,dim):
                  HSD[ia,kcld] += -0.5*ints[k,l,e,d]*ts[e,i]*(a==c) # 19
                for m in range(0,Nelec):
                  HSD[ia,kcld] += -0.5*ints[m,l,c,d]*ts[a,m]*(i==k) # 20
                  for e in range(Nelec,dim):
                    HSD[ia,kcld] += ints[m,k,e,c]*ts[e,m]*(i==l)*(a==d) # 21
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
                  HDS[iajb,kc] += (i==k)*ints[a,b,c,e]*ts[e,j] - \
                                  (j==k)*ints[a,b,c,e]*ts[e,i] # 24
                  HDS[iajb,kc] += (b==c)*ints[k,a,e,j]*ts[e,i] - \
                                  (b==c)*ints[k,a,e,i]*ts[e,j] - \
                                  (a==c)*ints[k,b,e,j]*ts[e,i] + \
                                  (a==c)*ints[k,b,e,i]*ts[e,j]   # 27 
                  HDS[iajb,kc] += (b==c)*fs[k,e]*td[e,a,i,j] - \
                                  (a==c)*fs[k,e]*td[e,b,i,j]  # 29
                  HDS[iajb,kc] += ints[k,a,c,e]*td[e,b,i,j] - \
                                  ints[k,b,c,e]*td[e,a,i,j] # 34 
                  for f in range(Nelec,dim):
                    HDS[iajb,kc] += 0.5*(b==c)*ints[k,a,e,f]*td[e,f,i,j] - \
                                    0.5*(a==c)*ints[k,b,e,f]*td[e,f,i,j]  # 32
                    HDS[iajb,kc] += -(a==c)*ts[e,i]*ts[f,j]*ints[k,b,e,f] + \
                                     (b==c)*ts[e,i]*ts[f,j]*ints[k,a,e,f]   # 36
                for m in range(0,Nelec):
                  HDS[iajb,kc] += (a==c)*ints[k,m,i,j]*ts[b,m] - \
                                  (b==c)*ints[k,m,i,j]*ts[a,m]  # 25
                  HDS[iajb,kc] += (j==k)*ints[m,b,c,i]*ts[a,m] - \
                                  (j==k)*ints[m,a,c,i]*ts[b,m] - \
                                  (i==k)*ints[m,b,c,j]*ts[a,m] + \
                                  (i==k)*ints[m,a,c,j]*ts[b,m] #26
                  HDS[iajb,kc] += (j==k)*fs[m,c]*td[a,b,m,i] - \
                                  (i==k)*fs[m,c]*td[a,b,m,j] # 28
                  HDS[iajb,kc] += ints[k,m,c,j]*td[a,b,m,i] - \
                                  ints[k,m,c,i]*td[a,b,m,j] #35
                  for n in range(0,Nelec):
                    HDS[iajb,kc] += 0.5*(i==k)*ints[m,n,c,j]*td[a,b,m,n] - \
                                    0.5*(j==k)*ints[m,n,c,i]*td[a,b,m,n]  # 33 
                    HDS[iajb,kc] += (i==k)*ts[a,m]*ts[b,n]*ints[m,n,c,j] - \
                                    (j==k)*ts[a,m]*ts[b,n]*ints[m,n,c,i]  # 37
                  for e in range(Nelec,dim):
                    HDS[iajb,kc] += (i==k)*ints[a,m,c,e]*td[e,b,m,j] - \
                                    (i==k)*ints[b,m,c,e]*td[e,a,m,j] - \
                                    (j==k)*ints[a,m,c,e]*td[e,b,m,i] + \
                                    (j==k)*ints[b,m,c,e]*td[e,a,m,i] # 30
                    HDS[iajb,kc] += (b==c)*ints[k,m,i,e]*td[e,a,m,j] - \
                                    (b==c)*ints[k,m,j,e]*td[e,a,m,i] - \
                                    (a==c)*ints[k,m,i,e]*td[e,b,m,j] + \
                                    (a==c)*ints[k,m,j,e]*td[e,b,m,i]  # 31
                    HDS[iajb,kc] += (i==k)*ts[a,m]*ts[e,i]*ints[m,b,c,e] - \
                                    (i==k)*ts[b,m]*ts[e,i]*ints[m,a,c,e] - \
                                    (j==k)*ts[a,m]*ts[e,j]*ints[m,b,c,e] + \
                                    (j==k)*ts[b,m]*ts[e,j]*ints[m,a,c,e] #38 
                    HDS[iajb,kc] += (a==c)*ts[e,i]*ts[b,m]*ints[k,m,e,j] - \
                                    (b==c)*ts[e,i]*ts[a,m]*ints[k,m,e,j] - \
                                    (a==c)*ts[e,j]*ts[b,m]*ints[k,m,e,i] + \
                                    (b==c)*ts[e,j]*ts[a,m]*ints[k,m,e,i] # 39 
                    #HDS[iajb,kc] +=  -(j==k)*ts[e,i] *ts[a,m]*ints[m,b,e,c] + (i==k) *ts[e,j]*ts[a,m] *ints[m,b,e,c] + (j==k)*ts[e,i]*ts[b,m]*ints[m,a,e,c] - (i==k) *ts[e,j]*ts[b,m]*ints[m,a,e,c] # 40
                    #HDS[iajb,kc] += (b==c)* ts[e,i] *ts[a,m]*ints[m,k,e,j] - (a==c)*ts[e,i]*ts[b,m]*ints[m,k,e,j] - (b==c)*ts[e,j]*ts[a,m]*ints[m,k,e,i] + (a==c)*ts[e,j]*ts[b,m]*ints[m,k,e,i] # 41 
                    HDS[iajb,kc] += ts[e,j]*td[a,b,m,i]*ints[k,m,c,e] - \
                                    ts[e,i]*td[a,b,m,j]*ints[k,m,c,e]   # 46 
                    HDS[iajb,kc] += ts[b,m]*td[e,a,i,j]*ints[k,m,c,e] - \
                                    ts[a,m]*td[e,b,i,j]*ints[k,m,c,e]  # 48 
                    for n in range(0,Nelec):
                      HDS[iajb,kc] += 0.5*(i==k)*ts[e,j]*td[a,b,m,n]*ints[m,n,c,e] - \
                                      0.5*(j==k)*ts[e,i]*td[a,b,m,n]*ints[m,n,c,e]  # 42
                      HDS[iajb,kc] += (j==k)*ts[a,m]*td[e,b,n,i]*ints[m,n,c,e] - \
                                      (i==k)*ts[a,m]*td[e,b,n,j]*ints[m,n,c,e] - \
                                      (j==k)*ts[b,m]*td[e,a,n,i]*ints[m,n,c,e] + \
                                      (i==k)*ts[b,m]*td[e,a,n,j]*ints[m,n,c,e]  # 44
                      HDS[iajb,kc] += (j==k)*ts[e,m]*td[a,b,n,i]*ints[m,n,e,c] - \
                                      (i==k)*ts[e,m]*td[a,b,n,j]*ints[m,n,e,c]  # 47
                      HDS[iajb,kc] += (j==k)*ts[e,i]*ts[a,m]*ts[b,n]*ints[m,n,e,c] - \
                                      (i==k)*ts[e,j]*ts[a,m]*ts[b,n]*ints[m,n,e,c] # 50
                    for f in range(Nelec,dim):
                      HDS[iajb,kc] += 0.5*(a==c)*ts[b,m]*td[e,f,i,j]*ints[k,m,e,f] - \
                                      0.5*(b==c)*ts[a,m]*td[e,f,i,j]*ints[k,m,e,f]  # 43
                      HDS[iajb,kc] += (b==c)*ts[e,i]*td[f,a,m,j]*ints[k,m,e,f] - \
                                      (a==c)*ts[e,i]*td[f,b,m,j]*ints[k,m,e,f] - \
                                      (b==c)*ts[e,j]*td[f,a,m,i]*ints[k,m,e,f] + \
                                      (a==c)*ts[e,j]*td[f,b,m,i]*ints[k,m,e,f]  # 45
                      HDS[iajb,kc] += (b==c)*ts[e,m]*td[f,a,i,j]*ints[m,k,e,f] - \
                                      (a==c)*ts[e,m]*td[f,b,i,j]*ints[m,k,e,f]  # 49
                      HDS[iajb,kc] += (b==c)*ts[e,i]*ts[a,m]*ts[f,j]*ints[m,k,e,f] - \
                                      (a==c)*ts[e,i]*ts[b,m]*ts[f,j]*ints[m,k,e,f]  # 51
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
                      HDD[iajb,kcld] += (a==c)*(i==l)*(b==d)*fs[k,e]*ts[e,j] - \
                                        (a==c)*(j==l)*(b==d)*fs[k,e]*ts[e,i] # 57
                      HDD[iajb,kcld] += (j==l)*(b==d)*ints[a,k,e,c]*ts[e,i] - \
                                        (i==l)*(b==d)*ints[a,k,e,c]*ts[e,j] - \
                                        (j==l)*(a==d)*ints[b,k,e,c]*ts[e,i] + \
                                        (i==l)*(a==d)*ints[b,k,e,c]*ts[e,j] # 59
                      HDD[iajb,kcld] += 0.5*(b==d)*(a==c)*ints[k,l,e,j]*ts[e,i] - \
                                        0.5*(b==d)*(a==c)*ints[k,l,e,i]*ts[e,j] # 62
                      HDD[iajb,kcld] += 0.5*(b==d)*ints[k,l,c,e]*td[e,a,i,j] - \
                                        0.5*(a==d)*ints[k,l,c,e]*td[e,b,i,j] # 70 
                      for f in range(Nelec,dim):
                        HDD[iajb,kcld] += 0.25*(a==c)*(b==d)*ints[k,l,e,f]*td[e,f,i,j] # 66
                        HDD[iajb,kcld] += 0.5*(a==c)*(b==d)*ints[k,l,e,f]*ts[e,i]*ts[f,j] # 72
                    for m in range(0,Nelec):
                      HDD[iajb,kcld] += (i==k)*(j==l)*(a==d)*fs[m,c]*ts[b,m] - \
                                        (i==k)*(j==l)*(b==d)*fs[m,c]*ts[a,m] # 58
                      HDD[iajb,kcld] += (i==l)*(b==d)*ints[m,k,j,c]*ts[a,m] - \
                                        (j==l)*(b==d)*ints[m,k,i,c]*ts[a,m] - \
                                        (i==l)*(a==d)*ints[m,k,j,c]*ts[b,m] + \
                                        (j==l)*(a==d)*ints[m,k,i,c]*ts[b,m] # 60
                      HDD[iajb,kcld] += 0.5*(j==l)*(i==k)*ints[m,a,c,d]*ts[b,m] - \
                                        0.5*(j==l)*(i==k)*ints[m,b,c,d]*ts[a,m] # 61
                      HDD[iajb,kcld] += 0.5*(j==l)*ints[k,m,c,d]*td[a,b,m,i] - \
                                        0.5*(i==l)*ints[k,m,c,d]*td[a,b,m,j] # 68 
                      for n in range(0,Nelec):
                        HDD[iajb,kcld] += 0.25*(i==k)*(j==l)*ints[m,n,c,d]*td[a,b,m,n] # 65
                        HDD[iajb,kcld] += 0.5*(i==k)*(j==l)*ints[m,n,c,d]*ts[a,m]*ts[b,n]  # 73 
                      for e in range(Nelec,dim):
                        HDD[iajb,kcld] += (i==k)*(j==l)*(b==d)*ints[m,a,e,c]*ts[e,m] - \
                                          (i==k)*(j==l)*(a==d)*ints[m,b,e,c]*ts[e,m]  # 63 
                        HDD[iajb,kcld] += (a==c)*(i==l)*(b==d)*ints[m,k,e,j]*ts[e,m] - \
                                          (a==c)*(j==l)*(b==d)*ints[m,k,e,i]*ts[e,m]   # 64 
                        HDD[iajb,kcld] += (j==l)*(b==d)*ints[k,m,c,e]*td[e,a,m,i] - \
                                          (i==l)*(b==d)*ints[k,m,c,e]*td[e,a,m,j] - \
                                          (j==l)*(a==d)*ints[k,m,c,e]*td[e,b,m,i] + \
                                          (i==l)*(a==d)*ints[k,m,c,e]*td[e,b,m,j]  # 67 
                        HDD[iajb,kcld] += (l==j)*(a==d)*ints[m,k,e,c]*ts[e,i]*ts[b,m] - \
                                          (l==j)*(b==d)*ints[m,k,e,c]*ts[e,i]*ts[a,m] - \
                                          (l==i)*(a==d)*ints[m,k,e,c]*ts[e,j]*ts[b,m] + \
                                          (l==i)*(b==d)*ints[m,k,e,c]*ts[e,j]*ts[a,m] # 74 
                        for n in range(0,Nelec):
                          HDD[iajb,kcld] += 0.5*(i==k)*(j==l)*(a==d)*ints[m,n,e,c]*td[b,e,n,m] - \
                                            0.5*(i==k)*(j==l)*(b==d)*ints[m,n,e,c]*td[a,e,n,m] # 71
                          HDD[iajb,kcld] += (i==k)*(j==l)*(a==d)*ints[m,n,e,c]*ts[e,m]*ts[b,n] - \
                                            (i==k)*(j==l)*(b==d)*ints[m,n,e,c]*ts[e,m]*ts[a,n]# 76
                        for f in range(Nelec,dim):
                          HDD[iajb,kcld] += 0.5*(a==c)*(i==l)*(b==d)*ints[m,k,e,f]*td[f,e,j,m] - \
                                            0.5*(a==c)*(j==l)*(b==d)*ints[m,k,e,f]*td[f,e,i,m] # 69
                          HDD[iajb,kcld] += (a==c)*(i==l)*(b==d)*ints[m,k,e,f]*ts[f,j]*ts[e,m] - \
                                            (a==c)*(j==l)*(b==d)*ints[m,k,e,f]*ts[f,i]*ts[e,m]# 75
    t7 = time.time()
    print "\t\t** Finished doubles-doubles block. Time: ",t7-t6 

    eomMatrix = np.bmat([[HSS,HSD],[HDS,HDD]])

    print "\t\t** Begin full diagonalization"
    print "\t\t * Matrix dimension:  ", str(len(eomMatrix)) + "x" + str(len(eomMatrix)) 
    eomEVal,eomEVec = np.linalg.eig(eomMatrix)

    print "\nExcitations (eV):"
    for excitation in (np.sort(np.real(eomEVal))*27.21138386):
      if (excitation > 1.0) and (excitation < 50.0):
        print (excitation)
        continue
    plt.imshow(eomMatrix,interpolation='nearest',cmap='jet',alpha=0.75)
    plt.colorbar()
    plt.savefig('eom-ccsd.png',bbox_inches='tight')
    plt.close()
    del eomMatrix, eomEVal, eomEVec
