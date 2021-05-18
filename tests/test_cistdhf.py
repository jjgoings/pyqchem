import numpy as np
from routines.integral_input import __init_integrals__ 
import routines.scf as scf
import routines.ao2mo as ao2mo
import routines.cistdhf as cistdhf


def do_cistdhf(LOCATION):
    ENUC,Nelec,dim,S,T,V,Hcore,twoe = __init_integrals__(LOCATION)
    EN,orbitalE,C,P,F = scf.scf_iteration(1e-8,ENUC,Nelec,dim,S,Hcore,twoe,printops=True,do_DIIS=True)
    ints,fs = ao2mo.transform_ao2mo(dim,twoe,C,orbitalE)
    ECIS, ETDHF = cistdhf.cistdhf(Nelec,dim,fs,ints,printnum=100000)
    
    return ECIS, ETDHF 

def test_cistdhf():
    """
    G16 CIS output:
    Excited State   1:  3.000-SGU   10.3954 eV  119.27 nm  f=0.0000  <S**2>=2.000
    Excited State   2:  1.000-SGU   15.7554 eV   78.69 nm  f=0.7983  <S**2>=0.000
    Excited State   3:  3.000-SGG   25.9084 eV   47.85 nm  f=0.0000  <S**2>=2.000
    Excited State   4:  1.000-SGG   32.1280 eV   38.59 nm  f=0.0000  <S**2>=0.000
    Excited State   5:  3.000-SGU   39.9982 eV   31.00 nm  f=0.0000  <S**2>=2.000
    Excited State   6:  1.000-SGU   46.6228 eV   26.59 nm  f=0.0714  <S**2>=0.000

    G16 TDHF output:
    Excited State   1:  3.000-SGU    9.8363 eV  126.05 nm  f=0.0000  <S**2>=2.000
    Excited State   2:  1.000-SGU   15.4997 eV   79.99 nm  f=0.6685  <S**2>=0.000
    Excited State   3:  3.000-SGG   25.7211 eV   48.20 nm  f=0.0000  <S**2>=2.000
    Excited State   4:  1.000-SGG   31.9771 eV   38.77 nm  f=0.0000  <S**2>=0.000
    Excited State   5:  3.000-SGU   39.7876 eV   31.16 nm  f=0.0000  <S**2>=2.000
    Excited State   6:  1.000-SGU   46.4065 eV   26.72 nm  f=0.0381  <S**2>=0.000
    """
    ref_ecis = np.array([10.3954, 10.3954, 10.3954, 15.7554, 25.9084, 25.9084, 25.9084, 32.1280, 39.9982, 39.9982, 39.9982, 46.6228])
    ref_etdhf = np.array([9.8363, 9.8363, 9.8363, 15.4997, 25.7211, 25.7211, 25.7211, 31.9771, 39.7876, 39.7876, 39.7876, 46.4065])
    ECIS, ETDHF = do_cistdhf('h2_3-21g')
    assert np.allclose(ECIS,ref_ecis)
    assert np.allclose(ETDHF,ref_etdhf)


    """"
    G16 CIS output:
    Excited State   1:  3.000-B1     7.8166 eV  158.62 nm  f=0.0000  <S**2>=2.000
    Excited State   2:  3.000-A1     9.3723 eV  132.29 nm  f=0.0000  <S**2>=2.000
    Excited State   3:  1.000-B1     9.6998 eV  127.82 nm  f=0.0023  <S**2>=0.000
    Excited State   4:  3.000-A2     9.9591 eV  124.49 nm  f=0.0000  <S**2>=2.000
    Excited State   5:  3.000-B2    10.7353 eV  115.49 nm  f=0.0000  <S**2>=2.000
    Excited State   6:  1.000-A2    11.3219 eV  109.51 nm  f=0.0000  <S**2>=0.000
    Excited State   7:  1.000-A1    13.7588 eV   90.11 nm  f=0.0649  <S**2>=0.000
    Excited State   8:  3.000-B2    13.9945 eV   88.59 nm  f=0.0000  <S**2>=2.000
    Excited State   9:  1.000-B2    15.1075 eV   82.07 nm  f=0.0155  <S**2>=0.000
    Excited State  10:  3.000-A1    15.3215 eV   80.92 nm  f=0.0000  <S**2>=2.000

    G16 TDHF output:
    Excited State   1:  3.000-B1     7.7597 eV  159.78 nm  f=0.0000  <S**2>=2.000
    Excited State   2:  3.000-A1     8.1564 eV  152.01 nm  f=0.0000  <S**2>=2.000
    Excited State   3:  3.000-B2     9.5955 eV  129.21 nm  f=0.0000  <S**2>=2.000
    Excited State   4:  1.000-B1     9.6540 eV  128.43 nm  f=0.0021  <S**2>=0.000
    Excited State   5:  3.000-A2     9.9357 eV  124.79 nm  f=0.0000  <S**2>=2.000
    Excited State   6:  1.000-A2    11.3014 eV  109.71 nm  f=0.0000  <S**2>=0.000
    Excited State   7:  1.000-A1    13.6084 eV   91.11 nm  f=0.0548  <S**2>=0.000
    Excited State   8:  3.000-B2    13.8958 eV   89.22 nm  f=0.0000  <S**2>=2.000
    Excited State   9:  3.000-A1    14.8594 eV   83.44 nm  f=0.0000  <S**2>=2.000
    Excited State  10:  1.000-B2    15.0036 eV   82.64 nm  f=0.0140  <S**2>=0.000
    """

    ref_ecis = np.array([7.8166,7.8166,7.8166,9.3723,9.3723,9.3723,9.6998,9.9591,9.9591,9.9591,10.7353,10.7353,10.7353,11.3219,13.7588,13.9945,13.9945,13.9945,15.1075,15.3215,15.3215,15.3215])
    ref_etdhf = np.array([7.7597,7.7597,7.7597,8.1564,8.1564,8.1564,9.5955,9.5955,9.5955,9.6540,9.9357,9.9357,9.9357,11.3014,13.6084,13.8958,13.8958,13.8958,14.8594,14.8594,14.8594,15.0036])
    ECIS, ETDHF = do_cistdhf('h2o_sto3g')
    assert np.allclose(ECIS[:22],ref_ecis)
    assert np.allclose(ETDHF[:22],ref_etdhf)
  



