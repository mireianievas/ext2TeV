import scipy.interpolate as sinterp
import numpy as np
from astropy.table import Table
from astropy import units as u

import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.optimize as sop
import scipy.integrate as sint

from ext2TeV.spectra import Spectrum

class Source(Spectrum):
    def __init__(self):
        super().__init__()
        self.fermi  = None
        self.fermi2 = None
        self.vhe    = None
    
    def add_fermi(self,spec):
        self.fermi = spec

    def add_fermi2(self,spec):
        self.fermi2 = spec
    
    def add_vhe(self,spec):
        self.vhe = spec

    def set_ebl(self,ebl):
        self.ebl = ebl
        
    def build_broadband_spectra(self):
        self.gamma = np.asarray([{
            'name_vhe': self.vhe[0]['name'],
            'name_fermi': self.fermi['name'],
            'redshift': self.fermi['redshift'],
            'class': self.fermi['class'],
            'E':  np.append(self.fermi['E'],self.vhe[k]['E']),
            'Eerr':  np.append(self.fermi['Eerr'],self.vhe[k]['Eerr'],axis=1),
            'EF': np.append(self.fermi['EF'],self.vhe[k]['EF']),
            'EFerr': np.append(self.fermi['EFerr'],self.vhe[k]['EFerr'],axis=1),
            'is_ul': np.append(self.fermi['EFerr'],self.vhe[k]['EFerr'],axis=1)[1]==0,
        } for k in range(len(self.vhe))])

    def get_state(self):
        '''
        Fit broadband gamma data using an EBL-absorbed LP to obtain the lowest and highest states
        '''
        z = self.fermi['redshift']
        logLP = lambda x, *args: np.log10((10**x)*(10**x)*self.logparabola(10**x,args))
        logLPabs = lambda x, *args: np.log10((10**x)*(10**x)*self.logparabola(10**x,args)) +\
            np.log10(self.ebl.ebl_absorption(z, 10**x))
        
        p0f = [-11., -1.5, 0.1]
        
        self.integrated_fluxes = [None for k in self.vhe]
        for k,gamma in enumerate(self.gamma):
            is_ul = gamma['is_ul']
            _xfit = np.log10(gamma['E'].to("TeV").value)[~is_ul]
            _yfit = np.log10(gamma['EF'].to("erg/(cm2*s)").value)[~is_ul]
            _yerrp = np.log10(gamma['EFerr'].to("erg/(cm2*s)").value)[1][~is_ul]
            
            # Sort by energy, otherwise bispev err may happen if there is a LAT point after the first VHE point.
            _ksort = np.argsort(_xfit)
            _xfit = _xfit[_ksort]
            _yfit = _yfit[_ksort]
            
            # do not use weights ... 
            _yerrfake = np.ones(_yfit.shape) * 1e-13  # np.log10(vhe['eflux'].value+vhe['efluxerr'][1].value)-_yfit
            popt, pcov = sop.curve_fit(logLPabs, _xfit, _yfit, p0f, sigma=_yerrfake)
            
            # Integrate
            logLPdFdE = lambda x: self.logparabola(10**x,popt)
            self.integrated_fluxes[k] = sint.quad(logLPdFdE, -1, 1)[0]
        
        self.index_lowstate  = np.argmin(self.integrated_fluxes)
        self.index_highstate = np.argmax(self.integrated_fluxes)
        #print(self.integrated_fluxes)
    