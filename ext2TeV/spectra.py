import scipy.interpolate as sinterp
import numpy as np
from astropy.table import Table
from astropy import units as u

import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.optimize as sop
import scipy.integrate as sint

class Spectrum(object):
    def __init__(self):
        self.redshift = 0
        self.E0 = 0.1
    
    def set_model_redshift(self,z):
        self.redshift = z
    
    def set_model_spectralclass(self,sptype):
        if sptype in ['EHSP','EHBL']:
            self.Ecutoff = 10
        elif sptype in ['HSP','HBL']:
            self.Ecutoff = 11
        else:
            self.Ecutoff = 0.1
        
        self.Ecutoff *= 1./(1+self.redshift)
        
    def powerlaw(self, E, theta):
        N0, Gamma = theta
        sed = 10**N0 * (E/self.E0)**Gamma
        return(sed)
    def logparabola(self, E, theta):
        beta = theta[-1]
        sed = self.powerlaw(E,theta[:-1])*(E/self.E0)**(-beta*np.log(E/self.E0))
        return(sed)
    def pwl_expcutoff(self, E, theta):
        Ecut = 10**theta[-1]
        sed = self.powerlaw(E,theta[:-1])*np.exp(-(E/Ecut))
        return(sed)
    def pwl_fermiexpcutoff(self, E, theta):
        Ecut = 10**theta[-1]
        sed = self.powerlaw(E,theta[:-1])*np.exp(-(E/Ecut)**(2/3.))
        return(sed)
    def pwl_supexpcutoff(self, E, theta):
        Ecut,Gamma2 = 10**theta[-2],theta[-1]
        sed = self.powerlaw(E,theta[:-2])*np.exp(-(E/Ecut)**Gamma2)
        return(sed)
    def logparabola_cutoff(self,E,theta):
        Ecut = 10**theta[-1]
        sed = self.logparabola(E,theta[:-1])*np.exp(-(E/Ecut)**Gamma2)
        return(sed)
    def logparabola_ctacut(self,E,theta):        
        sed = self.logparabola(E,theta)*np.exp(-(E/self.Ecutoff))
        return(sed)
    
    @staticmethod
    def lat_model_interpreter(E,model):
        ### lat energies are in MeV by default ... accept TeV?
        E = E*u.Unit("TeV")
        
        PivotE = model['PivotE']
        Fnorm = model['Fnorm']
        EFnorm = model['EFnorm']
        Gamma = model['Gamma']
        EGamma = model['EGamma']
        if 'Beta' in model:
            Beta = model['Beta']
            EBeta = model['EBeta']
        else:
            Beta = 0
            EBeta = 0
        if 'ExpF' in model:
            ExpF = model['ExpF']
            EExpF = model['EExpF']
            ExpI = model['ExpI']
            EExpI = model['EExpI']
        else:
            ExpF = 0
            EExpF = 0
            ExpI = 0
            EExpI = 0
    
        # (1.60218 * 1e-6) = TeV2/cm2/s/MeV -> erg/cm2/s 
        # attention, it is the natural log
        sed = Fnorm * (E ** 2) * \
              ((E / PivotE) ** (-Gamma - Beta * np.log(E / PivotE)))*\
              np.exp(ExpF * (PivotE.to("MeV").value ** ExpI - (E.to("MeV").value) ** ExpI))
        
        return(sed.to("erg/(cm2*s)").value)