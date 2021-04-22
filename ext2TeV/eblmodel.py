import scipy.interpolate as sinterp
import numpy as np
from astropy.table import Table
from astropy import units as u

import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.optimize as sop
import scipy.integrate as sint

class Eblmodel(object):
    def load_ebl_dominguez(self):
        # Load EBL table from Dominguez et al. 2011:
        eblmodeltab = Table.read(
            "https://raw.githubusercontent.com/me-manu/ebltable/"+\
            "master/ebltable/data/tau_dominguez11_cta.out",
            format='ascii')
        eblstruc = eblmodeltab.as_array()
        eblarray = eblstruc.view(np.float64).reshape(eblstruc.shape + (-1,))

        self.eblmodel = dict()
        self.eblmodel['redshifts'] = eblarray[0][1:]
        self.eblmodel['energies'] = np.log10(eblarray[:, 0][1:])
        self.eblmodel['tau2d'] = eblarray[1:, 1:]
        self.eblmodel['E0']  = 0.1 # 100 GeV
    
    def set_interpolator(self):
        self.ebl_tau_interp = sinterp.RectBivariateSpline(
            self.eblmodel['redshifts'], 
            self.eblmodel['energies'], 
            self.eblmodel['tau2d'].transpose()
        )
    
    def ebl_absorption(self,z,E): 
        try:
            self.ebl_tau_interp
        except AttributeError:
            self.load_ebl_dominguez()
            self.set_interpolator()
        
        return np.exp(-self.ebl_tau_interp(z, np.log10(E)))[0]

