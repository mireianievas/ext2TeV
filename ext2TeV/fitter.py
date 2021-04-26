import scipy.interpolate as sinterp
import numpy as np
from astropy.table import Table
from astropy import units as u

import matplotlib as mpl
import matplotlib.pyplot as plt
import scipy.optimize as sop
import scipy.integrate as sint

from ext2TeV.source import Source

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
                    label  = '{0:.10s} | lowest | {1}'.format('VHE',spec['name'])
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
                        label = '{0:.10s} | elevated | {1}'.format('VHE',spec['name'])
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
                            'ids': self.gamma[0]['name_vhe'].replace(" ","_"),
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
                            'ids': self.gamma[0]['name_vhe'].replace(" ","_"),
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
                    'model': 'LP+EBL+CTAGammaProp-Cutoff',
                    'xlogval': E_log[kp],
                    'ylogval': EF_log[kp],
                    'ylogerr': EFerrp_log[kp],
                    'ylogmod': EF_model[kp],
                    'ids': self.gamma[0]['name_vhe'].replace(" ","_"),
                }
                residuals.add_value(Point)
        
        self.plo.plot(
            Econt,
            10**logELPcta(np.log10(Econt),*pars),
            label='Fermi/LP+EBL+CTAGammaProp-Cutoff',
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
