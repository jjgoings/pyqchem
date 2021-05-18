import numpy as np
from routines.integral_input import __init_integrals__ 
import routines.scf as scf
import routines.ao2mo as ao2mo
import routines.eommbpt2 as eommbpt2


def do_eommbpt2(LOCATION):
    ENUC,Nelec,dim,S,T,V,Hcore,twoe = __init_integrals__(LOCATION)
    EN,orbitalE,C,P,F = scf.scf_iteration(1e-8,ENUC,Nelec,dim,S,Hcore,twoe,printops=True,do_DIIS=True)
    ints,fs = ao2mo.transform_ao2mo(dim,twoe,C,orbitalE)
    T1 = np.zeros((dim*2,dim*2))
    T2 = np.zeros((dim*2,dim*2,dim*2,dim*2))
    for a in range(Nelec,dim*2):
        for b in range(Nelec,dim*2):
            for i in range(0,Nelec):
                for j in range(0,Nelec):
                    T2[a,b,i,j] += ints[i,j,a,b]/(fs[i,i] + fs[j,j] - fs[a,a] - fs[b,b])
    EOMMBPT2 = eommbpt2.eommbpt2(Nelec,dim,fs,ints,T1,T2)
    return EOMMBPT2 

def test_eommbpt2():
    """
    G16 EOMMBPT2 output: [IOP(9/128=1), NoSymm]
    Excited State   1:      ?-?Sym   10.6572 eV  116.34 nm  f=-0.0000
    Excited State   2:      ?-?Sym   15.7087 eV   78.93 nm  f=0.6738
    Excited State   3:      ?-?Sym   26.2655 eV   47.20 nm  f=-0.0000
    Excited State   4:      ?-?Sym   30.2223 eV   41.02 nm  f=0.0000
    Excited State   5:      ?-?Sym   31.6785 eV   39.14 nm  f=0.0000
    Excited State   6:      ?-?Sym   40.2073 eV   30.84 nm  f=-0.0000
    """
    ref_eommbpt2 = np.array([10.6572, 10.6572, 10.6572, 15.7087, 26.2655, 26.2655, 26.2655, 30.2223, 31.6785, 40.2073, 40.2073, 40.2073])
    EOMMBPT2 = do_eommbpt2('h2_3-21g')
    assert np.allclose(EOMMBPT2[:12],ref_eommbpt2)


