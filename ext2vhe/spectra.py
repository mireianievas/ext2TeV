from astropy.table import Table, Column, vstack
import scipy.interpolate as sinterp

import numpy as np


class Spectra(object):

    def __init__(self):


        # Load EBL table from Dominguez et al. 2011:
        eblmodel = Table.read(
            "https://raw.githubusercontent.com/me-manu/ebltable/master/ebltable/data/tau_dominguez11_cta.out",
            format='ascii')
        eblstruc = eblmodel.as_array()
        eblarray = eblstruc.view(np.float64).reshape(eblstruc.shape + (-1,))

        ebl_absorption = dict()
        ebl_absorption['redshifts'] = eblarray[0][1:]
        ebl_absorption['energies'] = np.log10(eblarray[:, 0][1:])
        ebl_absorption['tau2d'] = eblarray[1:, 1:]
        z2d, e2d = np.meshgrid(ebl_absorption['redshifts'], ebl_absorption['energies'])
        ebl_absorption['z2d'] = z2d
        ebl_absorption['e2d'] = e2d
        # Normalization energy:
        self.E0 = 0.1  # 100 GeV
        self.ebl_absorption = ebl_absorption

    def ebl_absorption(self, zz, ee):
        # eblatten = np.exp(-sinterp.griddata((z2d.ravel(),e2d.ravel()),
        #                                tau2d.ravel(),
        #                                (zz,ee),
        #                                method='cubic',
        #                                fill_value=0))

        try:
            self.ebl_eval
        except:
            self.ebl_eval = sinterp.RectBivariateSpline(self.ebl_absorption['redshifts'],
                                                        self.ebl_absorption['energies'],
                                                        self.ebl_absorption['tau2d'], kx=1, ky=1)

        eblatten = self.ebl_eval(zz, ee)
        # eblatten[ee>np.max(energies)] = 0
        return np.exp(-eblatten[0])

    def power_law(self, E, theta):
        N0, Gamma = theta
        dNdE = 10 ** N0 * (E / 10 ** self.E0) ** Gamma
        sed = np.power(E, 2) * dNdE
        return (sed)

    def logparabola(self, E, theta):
        N0, Gamma, beta = theta
        beta = beta ** 2
        dNdE = np.power(10, N0) * np.power(np.power(E / 10, self.E0),
                                           Gamma - beta * np.log10(E / np.power(10, self.E0)))
        sed = np.power(E, 2) * dNdE
        return sed

    def pwl_exp_supcutoff(self, E, theta):
        N0, Gamma1, logEcutoff, Gamma2 = theta
        dNdE = 10 ** N0 * (E / 10 ** self.E0) ** (Gamma1) * np.exp(-(E / 10 ** logEcutoff) ** Gamma2)
        sed = np.power(E, 2) * dNdE
        return sed

    @staticmethod
    def flexible_model(e, model):
        # Not prepared for fitting.
        # This model, depending on the Beta and ExpI/F can represent a pure powerlaw
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

        # SED = Fnorm * E * E * (E/PivotE)**(-Gamma-Beta*np.log(E/PivotE)) * \
        #    np.exp(-ExpF*(PivotE**ExpI - E**Beta)) * 1.60218*1e-6 #erg/cm2/s

        # logSED = np.log10(Fnorm) + 2*np.log10(E) -\
        #         np.log10(E/PivotE)*(Gamma + Beta*np.log(E/PivotE)) -\
        #         ExpF*(PivotE**ExpI - E**Beta)/np.log(10) +\
        #         np.log10(1.60218*1e-6)

        sed = (1.60218 * 1e-6) * Fnorm * (e ** 2) * \
              ((e / PivotE) ** (-Gamma - Beta * np.log(e / PivotE))) * \
              np.exp(ExpF * (PivotE ** ExpI - e ** ExpI))

        return sed
