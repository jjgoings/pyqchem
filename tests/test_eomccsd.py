import numpy as np
from routines.integral_input import __init_integrals__ 
import routines.scf as scf
import routines.ao2mo as ao2mo
import routines.ccsd as ccsd
import routines.eomccsd as eomccsd


def do_eomccsd(LOCATION):
    ENUC,Nelec,dim,S,T,V,Hcore,twoe = __init_integrals__(LOCATION)
    EN,orbitalE,C,P,F = scf.scf_iteration(1e-8,ENUC,Nelec,dim,S,Hcore,twoe,printops=True,do_DIIS=True)
    ints,fs = ao2mo.transform_ao2mo(dim,twoe,C,orbitalE)
    ECCSD,T1,T2 = ccsd.ccsd(Nelec,dim,fs,ints,1e-8,printops=True)
    EOMCCSD = eomccsd.eomccsd(Nelec,dim,fs,ints,T1,T2)
    return EOMCCSD 

def test_eomccsd():
    """
    G16 EOMCCSD output:
    Excited State   1:      ?-?Sym   10.8527 eV  114.24 nm  f=-0.0000
    Excited State   2:      ?-?Sym   15.8984 eV   77.99 nm  f=0.6531
    Excited State   3:      ?-?Sym   26.4712 eV   46.84 nm  f=-0.0000
    Excited State   4:      ?-?Sym   30.5216 eV   40.62 nm  f=-0.0000
    Excited State   5:      ?-?Sym   31.8814 eV   38.89 nm  f=0.0000
    Excited State   6:      ?-?Sym   40.4020 eV   30.69 nm  f=-0.0000
    """
    ref_eomccsd = np.array([10.8527,10.8527,10.8527,15.8984,26.4712,26.4712,26.4712,30.5216,31.8814,40.4020,40.4020,40.4020])
    EOMCCSD = do_eomccsd('h2_3-21g')
    assert np.allclose(EOMCCSD[:12],ref_eomccsd)



