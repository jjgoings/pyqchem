import numpy as np
from routines.integral_input import __init_integrals__ 
import routines.scf as scf
import routines.ao2mo as ao2mo
import routines.ccsd as ccsd


def do_ccsd(LOCATION):
    ENUC,Nelec,dim,S,T,V,Hcore,twoe = __init_integrals__(LOCATION)
    EN,orbitalE,C,P,F = scf.scf_iteration(1e-8,ENUC,Nelec,dim,S,Hcore,twoe,printops=True,do_DIIS=True)
    ints,fs = ao2mo.transform_ao2mo(dim,twoe,C,orbitalE)
    ECCSD,T1,T2 = ccsd.ccsd(Nelec,dim,fs,ints,1e-8,printops=True)
    return ENUC + EN + ECCSD 

def test_ccsd():
    ccsd_energy = do_ccsd('h2o_sto3g')
    assert np.allclose(-75.012696756,ccsd_energy)
   
    ccsd_energy = do_ccsd('h2_3-21g')
    assert np.allclose(-1.1478131290,ccsd_energy)
   



