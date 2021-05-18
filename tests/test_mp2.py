import numpy as np
from routines.integral_input import __init_integrals__ 
import routines.scf as scf
import routines.ao2mo as ao2mo
import routines.mp2 as mp2


def do_mp2(LOCATION):
    # obviously, this will also test SCF
    ENUC,Nelec,dim,S,T,V,Hcore,twoe = __init_integrals__(LOCATION)
    EN,orbitalE,C,P,F = scf.scf_iteration(1e-8,ENUC,Nelec,dim,S,Hcore,twoe,printops=True,do_DIIS=True)
    ints,fs = ao2mo.transform_ao2mo(dim,twoe,C,orbitalE)
    mp2_corr = mp2.mp2_energy(dim,Nelec,ints,orbitalE)
    return ENUC + EN + mp2_corr

def test_mp2():
    mp2_energy = do_mp2('h2o_sto3g')
    assert np.allclose(-74.99122954257679,mp2_energy)
   
    mp2_energy = do_mp2('h2_3-21g')
    assert np.allclose(-1.1402533129128027,mp2_energy)
   
    mp2_energy = do_mp2('h2_sto3g')
    assert np.allclose(-1.1298973809859785,mp2_energy)

    mp2_energy = do_mp2('h2o_3-21g')
    assert np.allclose(-75.69464050746087,mp2_energy)

    mp2_energy = do_mp2('methane_sto3g')
    assert np.allclose(-39.78306066252666,mp2_energy)


