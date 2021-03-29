import os
from urllib.request import urlopen
#import reproject
from io import StringIO,BytesIO
from astropy.io import ascii
from astropy.io import fits as pyfits
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.modeling import models, fitting
from astropy.visualization.wcsaxes.frame import EllipticalFrame
from astropy.utils.data import get_pkg_data_filename
from astropy.table import Table,Column, vstack
from astropy.wcs import WCS
import matplotlib.cm as cm
#import astropy_healpix as ahp
#import healpy as hp
import numpy as np
import scipy.interpolate as sinterp
import scipy.integrate as sint
import scipy.optimize as sop
import emcee
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import FuncFormatter

from ext2vhe.residuals import Residuals
from ext2vhe.spectra import Spectra
from ext2vhe.utils import logparabola


class Source(object):

    def __init__(self):
        self.xval = []
        self.xerrn = []
        self.xerrp = []
        self.yval = []
        self.yerrn = []
        self.yerrp = []
        self.residuals = Residuals()
        self.spectra = Spectra()
        # Read Fermi-LAT catalogs:
        self.fhl3 = pyfits.open("https://fermi.gsfc.nasa.gov/ssc/data/access/lat/3FHL/gll_psch_v13.fit")
        self.fgl4 = pyfits.open("https://fermi.gsfc.nasa.gov/ssc/data/access/lat/10yr_catalog/gll_psc_v26.fit")
        self.lac4 = pyfits.open("https://fermi.gsfc.nasa.gov/ssc/data/access/lat/4LACDR2/table-4LAC-DR2-h.fits")
        self.lac4_low = pyfits.open("https://fermi.gsfc.nasa.gov/ssc/data/access/lat/4LACDR2/table-4LAC-DR2-l.fits")

    def set_residuals(self, residuals):
        self.residuals = residuals

    @staticmethod
    def match_with_catalog(catalog, src, tol=0.5 * u.Unit("deg")):
        coord_catalog = SkyCoord(
            ra=catalog['RAJ2000'] * u.Unit("deg"),
            dec=catalog['DEJ2000'] * u.Unit("deg"))

        c = SkyCoord(
            ra=src['ra'] * u.Unit("deg"),
            dec=src['dec'] * u.Unit("deg")
        )
        idx, d2d, d3d = c.match_to_catalog_3d(coord_catalog)
        if d2d > tol:
            return None
        return catalog[idx]

    def set_ebl_interpolator(self):
        eblmodel = Table.read(
            "https://raw.githubusercontent.com/me-manu/ebltable/master/ebltable/data/tau_dominguez11_cta.out",
            format='ascii')
        eblstruc = eblmodel.as_array()
        eblarray = eblstruc.view(np.float64).reshape(eblstruc.shape + (-1,))
        redshifts = eblarray[0][1:]
        energies = np.log10(eblarray[:, 0][1:])
        tau2d = np.transpose(eblarray[1:, 1:])
        self.ebl_tau = sinterp.RectBivariateSpline(redshifts, energies, tau2d)
        self.ebl_absorption = lambda z, e: np.exp(-self.ebl_tau(z, np.log10(e)))

    def get_fermi_points(self, fermisource, cat='3FHL'):
        if cat == '3FHL':
            raw_energies = self.fhl3[4].data
        elif cat == '4FGL':
            raw_energies = np.asarray([50, 100, 300, 1e3, 3e3, 1e4, 3e4, 3e5]) * 1e-3
            raw_energies = np.transpose([raw_energies[:-1], raw_energies[1:], 0 * raw_energies[1:]])

        print(fermisource[0])

        fermi_energies = [10 ** np.mean(np.log10(k[:-1])) for k in raw_energies]
        fermi_energyerrn = np.asarray([fermi_energies[k] - fe[0] for k, fe in enumerate(raw_energies)])
        fermi_energyerrp = np.asarray([fe[1] - fermi_energies[k] for k, fe in enumerate(raw_energies)])
        fermi_fluxes = fermisource['Flux_Band'] * u.Unit("cm**-2 * s**-1")
        fermi_fluxerr = np.transpose(np.abs(fermisource['Unc_Flux_Band'])) * u.Unit("cm**-2 * s**-1")
        fermi_energies = fermi_energies * u.Unit("GeV")
        fermi_efluxes = (fermi_fluxes * fermi_energies).to("erg/(cm2*s)")
        fermi_efluxerr = (fermi_fluxerr * fermi_energies).to("erg/(cm2*s)")

        isnotnan = (~np.isnan(fermi_fluxes)) * (~np.isnan(fermi_fluxerr[0])) * (fermi_fluxes > 0)

        try:
            self.result
        except:
            self.result = dict()

        self.result['name'] = fermisource['Source_Name']
        self.result['energy'] = (fermi_energies[isnotnan]).to("TeV")
        self.result['energyerr'] = (
                    np.asarray([fermi_energyerrn[isnotnan], fermi_energyerrp[isnotnan]]) * u.Unit("GeV")).to("TeV")

        self.result['eflux'] = fermi_efluxes[isnotnan]
        self.result['efluxerr'] = np.asarray([fermi_efluxerr[0][isnotnan],
                                              fermi_efluxerr[1][isnotnan]]) * u.Unit("erg/(cm2*s)")

        '''
        for k,E in enumerate(self.result['energy']):
            self.xval.append(np.log10(self.result['energy'][k].value))
            self.xerrp.append(np.log10(self.result['energy'][k].value+self.result['energyerr'][1][k].value)-self.xval[-1])
            self.xerrn.append(self.xerrp[-1])
            self.yval.append(np.log10(self.result['eflux'][k].value))
            self.yerrp.append(np.log10(self.result['eflux'][k].value+self.result['efluxerr'][1][k].value)-self.yval[-1])
            self.yerrn.append(self.yval[-1]-np.log10(self.result['eflux'][k].value-self.result['efluxerr'][0][k].value))
        '''

        # Get parameters + redshift
        self.result['redshift'] = 0
        if cat == '3FHL' and 'REDSHIFT' in fermisource:
            self.result['redshift'] = fermisource['REDSHIFT']
        elif '4FGL' in self.result['name']:
            srcnamefilt = self.lac4[1].data['Source_Name'] == self.result['name']
            if np.sum(srcnamefilt) == 1:
                self.redshift = 0
                if srcnamefilt in self.lac4[1].data:
                    self.redshift = self.lac4[1].data[srcnamefilt]['Redshift']
                elif srcnamefilt in self.lac4_low[1].data:
                    self.redshift = self.lac4_low[1].data[srcnamefilt]['Redshift']
                if self.redshift >= 0:
                    self.result['redshift'] = self.redshift

            self.get_model_parameters(fermisource)

        return self.result

    def get_model_parameters(self, fermisource):
        # Get model parameters from 4FGL/4LAC:
        self.result['pwl'] = dict()
        self.result['lp'] = dict()
        self.result['plec'] = dict()
        # Powerlaw
        self.result['pwl']['PivotE'] = fermisource['Pivot_Energy']
        self.result['pwl']['Fnorm'] = fermisource['PL_Flux_Density']
        self.result['pwl']['EFnorm'] = fermisource['Unc_PL_Flux_Density']
        self.result['pwl']['Gamma'] = fermisource['PL_Index']
        self.result['pwl']['EGamma'] = fermisource['Unc_PL_Index']
        # LogParabola
        self.result['lp']['PivotE'] = fermisource['Pivot_Energy']
        self.result['lp']['Fnorm'] = fermisource['LP_Flux_Density']
        self.result['lp']['EFnorm'] = fermisource['Unc_LP_Flux_Density']
        self.result['lp']['Gamma'] = fermisource['LP_Index']
        self.result['lp']['EGamma'] = fermisource['Unc_LP_Index']
        self.result['lp']['Beta'] = fermisource['LP_beta']
        self.result['lp']['EBeta'] = fermisource['Unc_LP_beta']
        # Powerlaw with exponential cutoff
        self.result['plec']['PivotE'] = fermisource['Pivot_Energy']
        self.result['plec']['Fnorm'] = fermisource['PLEC_Flux_Density']
        self.result['plec']['EFnorm'] = fermisource['Unc_PLEC_Flux_Density']
        self.result['plec']['Gamma'] = fermisource['PLEC_Index']
        self.result['plec']['EGamma'] = fermisource['Unc_PLEC_Index']
        self.result['plec']['ExpF'] = fermisource['PLEC_Expfactor']
        self.result['plec']['EExpF'] = fermisource['Unc_PLEC_Expfactor']
        self.result['plec']['ExpI'] = fermisource['PLEC_Exp_Index']
        self.result['plec']['EExpI'] = fermisource['Unc_PLEC_Exp_Index']

    def get_vhe_points(self, vhesource):
        results = []
        for k, sed in enumerate(vhesource['sed']):
            energyerr = None

            # replace ECSV dump by its Table meaning
            sed = Table.read(sed, format='ascii.ecsv')

            if 'e_ref' in sed.keys():
                energies = sed['e_ref'].to("TeV")
            if 'dnde' in sed.keys():
                flux = sed['dnde']
            if 'dnde_err' in sed.keys():
                fluxerr = sed['dnde_err']
                fluxerr = [flux - (10 ** (2 * np.log10(flux) - np.log10(flux + fluxerr))), fluxerr]
            if 'dnde_errp' in sed.keys():
                fluxerr = [sed['dnde_errn'], sed['dnde_errp']]

            isnotnan = ~np.isnan(flux) * (flux > 0)
            energies = energies[isnotnan]
            flux = flux[isnotnan]
            fluxerr = np.asarray([fluxerr[0][isnotnan], fluxerr[1][isnotnan]])

            if 'e_min' in sed.keys() and 'e_max' in sed.keys():
                energyerrn = (energies - sed['e_min'][isnotnan]).to(energies.unit)
                energyerrp = (sed['e_max'][isnotnan] - energies).to(energies.unit)

            if energyerr is None:
                energylog = np.log10(energies.value)
                energylogsep = np.min(energylog[1:] - energylog[:-1])
                energyerrn = energies.value - (10 ** (energylog - energylogsep / 2.))
                energyerrp = 10 ** (energylog + energylogsep / 2.) - energies.value

            e2flux = energies ** 2 * flux
            e2fluxerr = energies ** 2 * fluxerr

            result = {}
            result['name'] = vhesource['common_name']
            result['energy'] = energies
            result['energyerr'] = np.asarray([energyerrn, energyerrp]) * energies.unit
            result['eflux'] = e2flux
            result['efluxerr'] = e2fluxerr
            results.append(result)

            '''
            for k,E in enumerate(result['energy']):
                self.xval.append(np.log10(result['energy'][k].value))
                self.xerrp.append(np.log10(result['energy'][k].value+result['energyerr'][1][k].value)-self.xval[-1])
                self.xerrn.append(self.xerrp[-1])
                self.yval.append(np.log10(result['eflux'][k].value))
                self.yerrp.append(np.log10(result['eflux'][k].value+result['efluxerr'][1][k].value)-self.yval[-1])
                self.yerrn.append(self.yval[-1]-np.log10(result['eflux'][k].value-result['efluxerr'][0][k].value))
            '''
        return results

    def broadband_sed(self, dataset, model='pwl', figure_path=''):
        redshift = None
        ### extract the sed
        if 'sed' not in dataset.keys():
            print('no sed')
            return None

        fig = plt.figure(figsize=(6, 3))
        # sedplot = fig.add_subplot(111)
        sedplot = fig.add_axes((.1, .3, .8, .6))
        resplot = fig.add_axes((.1, .1, .8, .2))

        self.list_of_spectra = dict()
        self.list_of_spectra['xlog'] = []
        self.list_of_spectra['ylog'] = []
        self.list_of_spectra['yerrlog'] = []

        vhe_array = self.get_vhe_points(dataset)

        ### Match with 4FGL
        matched = self.match_with_catalog(self.fgl4[1].data, dataset)
        if matched is not None:
            fermi = self.get_fermi_points(matched, '4FGL')
            # print(fermi)
            gammamodel = dict(fermi)  # make a copy to preserve the 4FGL spectral model.
            if (redshift == None): redshift = fermi['redshift']
            for j in range(len(fermi['energy'].value)):
                self.list_of_spectra['xlog'].append(fermi['energy'][j].value)
                self.list_of_spectra['ylog'].append(fermi['eflux'][j].value)
                self.list_of_spectra['yerrlog'].append(fermi['efluxerr'][1][j].value)

            sedplot.errorbar(x=fermi['energy'].value,
                             xerr=fermi['energyerr'].value,
                             y=fermi['eflux'].value,
                             yerr=fermi['efluxerr'].value,
                             ls='', marker='o', color='gray', ms=4, mfc='white',
                             label='{0} (4FGL)'.format(fermi['name']))

            xfit = fermi['energy'].value
            fitatten = self.ebl_absorption(redshift, xfit)[0]
            reasonable = fitatten > 1e-3
            xfit = xfit[reasonable]
            fitatten = fitatten[reasonable]
            yfit = self.spectra.flexible_model(xfit * 1e6, model=gammamodel[model]) * fitatten
            res = (fermi['eflux'].value[reasonable] - yfit) / yfit
            reserr = fermi['efluxerr'][1].value[reasonable] / yfit
            resplot.errorbar(x=xfit,
                             y=res,
                             yerr=reserr,
                             ls='', marker='o', color='gray', ms=4, mfc='white')


        else:
            print('no match in 4FGL')

        ### Match with 3FHL
        matched = self.match_with_catalog(self.fhl3[1].data, dataset)
        if matched is not None:
            fermi = self.get_fermi_points(matched, '3FHL')
            if (redshift == None): redshift = fermi['redshift']

            for j in range(len(fermi['energy'].value)):
                self.list_of_spectra['xlog'].append(fermi['energy'][j].value)
                self.list_of_spectra['ylog'].append(fermi['eflux'][j].value)
                self.list_of_spectra['yerrlog'].append(fermi['efluxerr'][1][j].value)

            sedplot.errorbar(x=fermi['energy'].value,
                             xerr=fermi['energyerr'].value,
                             y=fermi['eflux'].value,
                             yerr=fermi['efluxerr'].value,
                             ls='', marker='o', color='black', ms=4, mfc='white',
                             label='{0} (3FHL)'.format(fermi['name']))

            xfit = fermi['energy'].value
            fitatten = self.ebl_absorption(redshift, xfit)[0]
            reasonable = fitatten > 1e-3
            xfit = xfit[reasonable]
            fitatten = fitatten[reasonable]
            yfit = self.spectra.flexible_model(xfit * 1e6, model=gammamodel[model]) * fitatten
            res = (fermi['eflux'].value[reasonable] - yfit) / yfit
            reserr = fermi['efluxerr'][1].value[reasonable] / yfit
            resplot.errorbar(x=xfit,
                             y=res,
                             yerr=reserr,
                             ls='', marker='o', color='black', ms=4, mfc='white')

        else:
            print('no match in 3FHL')

        # Fit to a log-parabola to check which is the lowest state
        logLP = lambda x, *par: np.log10(x * x * logparabola(x, par))
        FlLogabove100_array = []

        for k, vhe in enumerate(vhe_array):
            logLPabs = lambda x, *args: logLP(10 ** x, *args) + np.log10(self.ebl_absorption(redshift, 10 ** x)[0])
            p0f = [-11., -2.]

            _xfit = np.log10(vhe['energy'].value)
            _yfit = np.log10(vhe['eflux'].value)
            _yerr = np.ones(_yfit.shape) * 1e-13  # np.log10(vhe['eflux'].value+vhe['efluxerr'][1].value)-_yfit

            popt, pcov = sop.curve_fit(logLPabs, _xfit, _yfit,
                                       p0f,
                                       sigma=_yerr)

            logLPfit = lambda x: logparabola(10 ** x, popt)
            FlLogabove100 = sint.quad(logLPfit, -0.5, 0.8)
            FlLogabove100_array.append(FlLogabove100[0])
            vhe['k'] = k

        print(FlLogabove100_array)
        k_lowstate = np.argmin(FlLogabove100_array)
        k_higstate = np.argmax(FlLogabove100_array)

        for k, vhe in enumerate(vhe_array):
            # plot every vhe spectrum / dataset
            for j in range(len(vhe['energy'].value)):
                self.list_of_spectra['xlog'].append(vhe['energy'][j].value)
                self.list_of_spectra['ylog'].append(vhe['eflux'][j].value)
                self.list_of_spectra['yerrlog'].append(vhe['efluxerr'][1][j].value)

            sedplot.errorbar(x=vhe['energy'].value,
                             xerr=vhe['energyerr'].value,
                             y=vhe['eflux'].value,
                             yerr=vhe['efluxerr'].value,
                             ls='', marker='D', ms=4, mfc='white' if vhe['k'] != k_lowstate else 'black',
                             label='{0}, {1} (VHE)'.format(vhe['name'], k + 1))

            if vhe['k'] == k_lowstate:
                for ki, E in enumerate(vhe['energy']):
                    self.xval.append(np.log10(vhe['energy'][ki].value))
                    self.xerrp.append(np.log10(vhe['energy'][ki].value + vhe['energyerr'][1][ki].value) - self.xval[-1])
                    self.xerrn.append(self.xerrp[-1])
                    self.yval.append(np.log10(vhe['eflux'][ki].value))
                    self.yerrp.append(np.log10(vhe['eflux'][ki].value + vhe['efluxerr'][1][ki].value) - self.yval[-1])
                    self.yerrn.append(self.yval[-1] - np.log10(vhe['eflux'][ki].value - vhe['efluxerr'][0][ki].value))

            xfit = vhe['energy'].value
            fitatten = self.ebl_absorption(redshift, xfit)[0]
            reasonable = fitatten > 1e-3
            xfit = xfit[reasonable]
            fitatten = fitatten[reasonable]
            yfit = self.spectra.flexible_model(xfit * 1e6, model=gammamodel[model]) * fitatten
            res = (vhe['eflux'].value[reasonable] - yfit) / yfit
            reserr = vhe['efluxerr'][1].value[reasonable] / yfit
            resplot.errorbar(x=xfit,
                             y=res,
                             yerr=reserr,
                             ls='', marker='D', ms=4, mfc='white')

        for model in ["pwl", "lp", "plec"]:
            # print(fermi[model])
            xfit = 10 ** np.linspace(-4, 2, 500)  # TeV
            fitatten = self.ebl_absorption(redshift, xfit)[0]
            xfit = xfit[fitatten > 1e-3]
            fitatten = fitatten[fitatten > 1e-3]
            yfit = self.spectra.flexible_model(xfit * 1e6, model=gammamodel[model]) * fitatten
            sedplot.plot(xfit, yfit, label=model)

            for k, xval in enumerate(self.xval):
                fitatten_p = self.ebl_absorption(redshift, 10 ** xval).flatten()
                ymodel_p = (self.spectra.flexible_model((10 ** xval) * 1e6, model=gammamodel[model]) * fitatten_p)[0]
                point = {
                    'srccls': dataset['classes'],
                    'model': model,
                    'xlogval': self.xval[k],
                    'ylogval': self.yval[k],
                    'ylogerr': self.yerrp[k] if self.yval[k] < ymodel_p else self.yerrm[k],
                    'ylogmod': np.log10(ymodel_p),
                }
                self.residuals.add_value(point)

        # box = sedplot.get_position()
        # sedplot.set_position([box.x0, box.y0,
        #                      box.width * 0.8, box.height])

        sedplot.set_xlim(45 * 1e-6, 100)
        resplot.set_xlim(45 * 1e-6, 100)
        sedplot.set_yscale('log')
        sedplot.set_xscale('log')
        sedplot.set_ylim([1e-15, 1e-10])
        sedplot.axes.get_xaxis().set_visible(False)
        # resplot.set_yscale('log')
        resplot.set_xscale('log')
        sedplot.set_ylabel('$\mathrm{SED}\ [erg\ cm^{-2} s^{-1}]$')
        sedplot.set_xlabel('$\log_{10}(\mathrm{Energy})\ [TeV]$')
        sedplot.legend(fontsize='xx-small',
                       loc='center left',
                       bbox_to_anchor=(1, 0.5))

        plt.tight_layout()
        # Make sure the output folder exists:
        if not os.path.exists(figure_path):
            os.mkdir(figure_path)
        plt.savefig("{}/FermiAndVHE{}.png".format(figure_path, vhe['name'].replace(" ", "")),
                    dpi=200, bbox_inches='tight')

        return fig
