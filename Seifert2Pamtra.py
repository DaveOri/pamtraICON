import numpy as np

def transform(a=np.nan,b=np.nan,nu=np.nan,mu=np.nan,xmin=np.nan,xmax=np.nan,A=np.nan,lam=np.nan):
    """
    This is a converter from the Axel Seifert MAP (2006) ice parametrization scheme to pamtra
    Anything that is not set or not computable goes to nan
    """
    bm  = 1.0/b
    am  = 1.0/a**(bm)
    N0  = am**2.0*bm*A
    mum = nu*bm+bm-1.0
    LAM = lam*am**mu
    gam = bm*mu
    Dmin= a*xmin**b
    Dmax= a*xmax**b

    print('am = ',  am)
    print('bm = ',  bm)
    print('N0 = ',  N0)
    print('mum= ', mum)
    print('LAM= ', LAM)
    print('gam= ', gam)
    print('Dmn= ',Dmin)
    print('Dmx= ',Dmax)
