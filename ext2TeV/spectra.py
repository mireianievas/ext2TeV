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
    
class Fitter(Source):
    def __init__(self):
        super().__init__()

    def create_figure(self):
        self.fig = plt.figure(figsize=(4,3), dpi=150)
        self.plo = self.fig.add_subplot(111)
        self.plo.set_yscale('log')
        self.plo.set_xscale('log')
        mpl.rcParams['text.usetex'] = False
    
    def finish_figure(self):
        self.plo.set_ylim(
            np.max([self.plo.get_ylim()[0],1e-15]),
            np.min([self.plo.get_ylim()[1],1e-8])
        )
        
        self.plo.set_ylabel('$\mathrm{E^2 dF/dE\, [erg/cm2/s]}$')
        self.plo.set_xlabel('$\mathrm{Energy\, [TeV]}$')
        self.plo.set_title('{0} ({1})'.format(self.gamma[0]['name_vhe'],self.gamma[0]['name_fermi']))
        self.plo.legend(fontsize='5',ncol=2,loc=3)
    
    def plot_fermi_data(self):
        for fermi in [self.fermi,self.fermi2]:
            if fermi is not None:
                self.plo.errorbar(
                    x=fermi['E'].to("TeV").value,
                    xerr=fermi['Eerr'].to("TeV").value,
                    y=fermi['EF'].to("erg/(cm2*s)").value,
                    yerr=fermi['EFerr'].to("erg/(cm2*s)").value,
                    ls='None',
                    marker='D',
                    ms=4,
                    mfc='white',
                    label='{0:.10s} | {1}'.format('Fermi',fermi['name'])
                )
        
    
    def plot_vhe_data(self):
        lowstatealreadylabel = None
        for k,spec in enumerate(self.vhe):
            try:
                self.index_lowstate
            except:
                mfc = 'white'
                marker = '+'
            else:
                if (k == self.index_lowstate):
                    mfc = 'white'
                    marker='o'
                    lw=1.5
                    xerr = spec['Eerr'].to("TeV").value
                    yerr = spec['EFerr'].to("erg/(cm2*s)").value
                    mfc  = 'white'
                    color  = 'C2'
                    alpha  = 1
                    label  = '{0:.10s} | lowest | {1}'.format('GammaCat',spec['name'])
                else:
                    color  = 'C2'
                    alpha  = 0.5
                    mfc    = 'C2'
                    marker = '.'
                    xerr = None
                    yerr = None
                    lw=0.8
                    if lowstatealreadylabel == 'has_label':
                        label = None
                    else:
                        label = '{0:.10s} | elevated | {1}'.format('GammaCat',spec['name'])
                        lowstatealreadylabel = 'has_label'
            
            self.plo.errorbar(
                x=spec['E'].to("TeV").value,
                xerr=xerr,
                y=spec['EF'].to("erg/(cm2*s)").value,
                yerr=yerr,
                ls='None',
                marker=marker,
                mfc=mfc,
                color=color,
                ms=4,
                lw=lw,
                label=label,
                alpha=alpha
            )
        
    
        
    def predefined_models(self,residuals=None,without_ebl=False,with_ebl=True):
        '''
        Fit the broadband data using available models
        '''
        # Predefined LAT models with fixed pars (do not fit, just calculate residuals)
        
        #for k,spec in  enumerate(self.gamma):
            #if k != self.index_lowstate: continue
        
        
        k = self.index_lowstate
        Econt = np.logspace(-4,1,100)
        is_ul = self.vhe[k]['EF'].value<=0
        
        k_sort = np.argsort(self.vhe[k]['E'].value[~is_ul])
        E = self.vhe[k]['E'].value[~is_ul][k_sort]
        EF = self.vhe[k]['EF'].to("erg/(cm2*s)").value[~is_ul][k_sort]
        EFerrp = self.vhe[k]['EFerr'].to("erg/(cm2*s)").value[1][~is_ul][k_sort]
        
        E_log = np.log10(E)
        EF_log = np.log10(EF)
        EFerrp_log = np.log10(EF+EFerrp)-EF_log

        if without_ebl:
            for model in ['pwl','lp','plec']:
                EF_model = np.log10(self.lat_model_interpreter(E,self.fermi[model]))            
                if residuals is not None:
                    for kp, P in enumerate(EF_model):
                        Point = {
                            'srccls': self.gamma[0]['class'],
                            'model': model,
                            'xlogval': E_log[kp],
                            'ylogval': EF_log[kp],
                            'ylogerr': EFerrp_log[kp],
                            'ylogmod': EF_model[kp],
                        }
                        residuals.add_value(Point)

            self.plo.plot(
                Econt,
                self.lat_model_interpreter(Econt,self.fermi['pwl']),
                label='Fermi/PWL',
                ls='dotted',
                color='hotpink',
            )

            self.plo.plot(
                Econt,
                self.lat_model_interpreter(Econt,self.fermi['lp']),
                label='Fermi/LP',
                ls='dashed',
                color='hotpink',
            )

            self.plo.plot(
                Econt,
                self.lat_model_interpreter(Econt,self.fermi['plec']),
                label='Fermi/PLEC',
                ls='dashdot',
                color='hotpink',
            )
            
            # Same, with EBL  
        
        if with_ebl:
            ebl_abs = self.ebl.ebl_absorption(self.redshift, E)
            for model in ['pwl','lp','plec']:
                EF_model = np.log10(self.lat_model_interpreter(E,self.fermi[model])*ebl_abs)
                if residuals is not None:
                    for kp, P in enumerate(EF_model):
                        Point = {
                            'srccls': self.gamma[0]['class'],
                            'model': model,
                            'xlogval': E_log[kp],
                            'ylogval': EF_log[kp],
                            'ylogerr': EFerrp_log[kp],
                            'ylogmod': EF_model[kp],
                        }
                        residuals.add_value(Point)
            
            ebl_abs = self.ebl.ebl_absorption(self.redshift, Econt)
            self.plo.plot(
                Econt,
                self.lat_model_interpreter(Econt,self.fermi['pwl'])*ebl_abs,
                label='Fermi/PWL+EBL',
                ls='dotted',
                color='crimson',
            )

            self.plo.plot(
                Econt,
                self.lat_model_interpreter(Econt,self.fermi['lp'])*ebl_abs,
                label='Fermi/LP+EBL',
                ls='dashed',
                color='crimson',
            )

            self.plo.plot(
                Econt,
                self.lat_model_interpreter(Econt,self.fermi['plec'])*ebl_abs,
                label='Fermi/PLEC+EBL',
                ls='dashdot',
                color='crimson',
            )
            
        # Following CTA's paper
        logELPcta = lambda x, *args: np.log10((10**x)*(10**x)*self.logparabola_ctacut(10**x,args)) +\
            np.log10(self.ebl.ebl_absorption(self.redshift, 10**x))
        self.E0 = self.fermi['lp']['PivotE'].to("TeV").value
        pars = [
            np.log10(self.fermi['lp']['Fnorm'].to("erg/(cm2*s*TeV2)").value),
            -self.fermi['lp']['Gamma'],
            self.fermi['lp']['Beta'],
        ]
        
        if residuals is not None:
            EF_model = logELPcta(E_log,*pars)
            for kp, P in enumerate(EF_model):
                Point = {
                    'srccls': self.gamma[0]['class'],
                    'model': 'LP+EBL+Biteau-Cutoff',
                    'xlogval': E_log[kp],
                    'ylogval': EF_log[kp],
                    'ylogerr': EFerrp_log[kp],
                    'ylogmod': EF_model[kp],
                }
                residuals.add_value(Point)
        
        self.plo.plot(
            Econt,
            10**logELPcta(np.log10(Econt),*pars),
            label='Fermi/LP+EBL+Biteau-Cutoff',
            ls='dashed',
            color='indigo',
        )
        
    
    def fit_models(self,residuals=None):
        # fit with common models
        logPWLabs = lambda x, *args: np.log10((10**x)*(10**x)*self.powerlaw(10**x,args)) +\
            np.log10(self.ebl.ebl_absorption(self.redshift, 10**x)[0])
        logLPabs = lambda x, *args: np.log10((10**x)*(10**x)*self.logparabola(10**x,args)) +\
            np.log10(self.ebl.ebl_absorption(self.redshift, 10**x)[0])
        logEPWLabs = lambda x, *args: np.log10((10**x)*(10**x)*self.pwl_expcutoff(10**x,args)) +\
            np.log10(self.ebl.ebl_absorption(self.redshift, 10**x)[0])
        logELPcta = lambda x, *args: np.log10((10**x)*(10**x)*self.logparabola_ctacut(10**x,args)) +\
            np.log10(self.ebl.ebl_absorption(self.redshift, 10**x)[0])        
        
        ### PWL
        model_name = 'PWL'
        p0f = [-11., -2.0]
        popt, pcov = sop.curve_fit(logPWLabs, E_log, EF_log, p0f, sigma=EFerrp_log/EFerrp_log)
        self.plo.plot(
            Econt,
            10**logPWLabs(np.log10(Econt),*popt),
            label='Fit/PWL',
            ls='dashdot',
            color='indigo',
        )
            
        ### LP
        model_name = 'LP'
        p0f = [-11., -1.5, 0.1]
        popt, pcov = sop.curve_fit(logLPabs, E_log, EF_log, p0f, sigma=EFerrp_log/EFerrp_log)
        self.plo.plot(
            Econt,
            10**logLPabs(np.log10(Econt),*popt),
            label='Fit/LP',
            ls='dashdot',
            color='indigo',
        )
        ### EPWL
        model_name = 'EPWL'
        p0f = [-11., -1.5, 0]
        popt, pcov = sop.curve_fit(logEPWLabs, E_log, EF_log, p0f, sigma=EFerrp_log/EFerrp_log)
        self.plo.plot(
            Econt,
            10**logEPWLabs(np.log10(Econt),*popt),
            label='Fit/EPWL',
            ls='dashdot',
            color='indigo',
        )

        ### ELPCTA
        model_name = 'ELPCTA'
        p0f = [-11., -1.5, 0.01]
        popt, pcov = sop.curve_fit(logELPcta, E_log, EF_log, p0f, sigma=EFerrp_log/EFerrp_log)
        self.plo.plot(
            Econt,
            10**logELPcta(np.log10(Econt),*popt),
            label='Fit/ELPCTA',
            ls='dashdot',
            color='indigo',
        )
            
    def savefig(self,filename=None):
        if filename == None:
            name_vhe = self.gamma[0]['name_vhe'].replace(" ","_")
            filename = "results/Spectra_{}_ext2TeV.pdf".format(name_vhe)

        self.fig.savefig(filename,bbox_inches='tight')
