import json
import os
import numpy as np

from astropy.io import ascii
from astropy.io import fits as pyfits
from astropy.io.misc import yaml
from astropy.io import fits as pyfits
from astropy.coordinates import SkyCoord
from astropy.table import Table
from astropy import units as u

from ext2TeV import PACKAGE_DATA
from ext2TeV.utils import read_yaml

class Catalog(object):
    def reset_spectral_points(self):
        self.xval = []
        self.xerrn = []
        self.xerrp = []
        self.yval = []
        self.yerrn = []
        self.yerrp = []
        #self.residuals = Residuals()
        #self.spectra   = Spectra()
    
    #@staticmethod
    def match_src_catalog(self, catalog, src, tol=0.5*u.Unit("deg")):
        coord_catalog = SkyCoord(
            ra=catalog['ra'] * u.Unit("deg"),
            dec=catalog['dec'] * u.Unit("deg")
        )

        c = SkyCoord(
            ra=src['ra'] * u.Unit("deg"),
            dec=src['dec'] * u.Unit("deg")
        )

        idx, d2d, d3d = c.match_to_catalog_3d(coord_catalog)
        #print(idx, d2d, d3d)
        if d2d > tol:
            return None

        return catalog['id'][idx]
    
    def reset_results(self):
        self.result = dict()
  
class GammaCAT(Catalog):
    def read(self,f):
        with open(f) as temp:
            self.contents = yaml.load(temp.read())
            #return(self.contents)
    
    def build(self,f,ind,outf="gammacat_export.yaml"):
        
        os.chdir(ind)
        with open("gammacat.yaml") as gammacat:
            temp = yaml.load(gammacat.read())
            gammacat = dict()
            for k,source_id in enumerate(temp):
                gammacat[source_id[0][1]] = dict()
                for pair in temp[k]:
                    gammacat[source_id[0][1]][pair[0]] = pair[1]
        
        with open("gammacat-datasets.json") as datasets:
            temp = json.load(datasets)
            datasets = dict()
            for item in temp:
                if item['source_id'] not in datasets:
                    datasets[item['source_id']] = dict()
                
                # copy identifiers from gammacat
                if item['source_id'] in gammacat:
                    print('Appending data from gammacat[{0}]'.format(item['source_id']))
                    for tocopy in gammacat[item['source_id']]:
                        datasets[item['source_id']][tocopy] = gammacat[item['source_id']][tocopy]
                else:
                    print('Fallback, getting the info from the sources dir')
                    info = read_yaml("sources/tev-{0:06d}.yaml".format(item['source_id']))
                    for tocopy in info:
                        datasets[item['source_id']][tocopy] = info[tocopy]

                if item['type'] == 'sed':
                    try:
                        content = open(item['location']).read() #Table.read(item['location'])
                    except:
                        print('Skipping {0}'.format(item['location']))

                    if 'sed' in datasets[item['source_id']]:
                        datasets[item['source_id']]['sed'].append(content)
                    else:
                        datasets[item['source_id']]['sed'] = [content]
                    
                    
        with open("gammacat_export_20210215.yaml","w+") as gcf:
            yaml.dump(datasets,gcf)
            
        self.contents = datasets
        
    def match_src(self,src,tol=0.5*u.Unit("deg")):
        list_items = self.contents.keys()
        catalog = {}
        catalog['dec'] = [self.contents[item]['dec'] \
                          for item in self.contents \
                          if 'dec' in self.contents[item]]
        catalog['ra']  = [self.contents[item]['ra'] \
                          for item in self.contents \
                          if 'dec' in self.contents[item]]
        catalog['id']  = [item for item in self.contents \
                          if 'dec' in self.contents[item]]
        
        try:
            self.matched = self.contents[self.match_src_catalog(catalog,src,tol)]
        except KeyError:
            self.matched = None
            return
        
    def get_spectralpoints(self):
        self.results = []
        for k, sed in enumerate(self.matched['sed']):
            energyerr = None
            # replace ECSV dump by its Table version
            sed = Table.read(sed, format='ascii.ecsv')
            #print(sed)
            energies = sed['e_ref'].to("TeV")
            flux = sed['dnde']
            if 'dnde_err' in sed.keys():
                fluxerr = sed['dnde_err']
                logFlux = 2 * np.log10(flux) - np.log10(flux + fluxerr)
                fluxerr = [flux - (10**(logFlux)), fluxerr]
            elif 'dnde_errp' in sed.keys():
                fluxerr = [sed['dnde_errn'], sed['dnde_errp']]
            else:
                print('No dnde_err*, check' )
                print(sed)
            
    
            isnotnan = ~np.isnan(flux) * (flux > 0)
            energies = energies[isnotnan]
            flux = flux[isnotnan]
            fluxerr = np.asarray([fluxerr[0][isnotnan], fluxerr[1][isnotnan]])
            
            # add units
            fluxerr = fluxerr*flux.unit
            
            
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
            result['name'] = self.matched['common_name']
            result['set'] = k
            result['E'] = energies.to("TeV")
            result['Eerr'] = np.asarray([energyerrn, energyerrp]) * energies.unit
            result['Eerr'] = result['Eerr'].to("TeV")
            result['EF'] = e2flux.to("erg/(cm2*s)")
            result['EFerr'] = e2fluxerr.to("erg/(cm2*s)")
            self.results.append(result)

class VTSCAT(Catalog):
    def read(self,f):
        with open(f) as temp:
            self.contents = yaml.load(temp.read())  
            #return(self.contents)
    
    def build(self,ind,outf="vtscat_export.yaml"):
        #os.chdir(ind)
        vts_sources_yamls = sorted(glob.glob(ind+"/sources_vtscat/tev-*.yaml"))
        VTS_srcs = {}
        for fn in vts_sources_yamls:
            with open(fn) as f:
                YAMLcontents = yaml.load(f.read())
                VTS_srcs[YAMLcontents['source_id']]=YAMLcontents
            
        for indx in VTS_srcs:
            VTS_srcs[indx]['seds'] = []
            VTS_srcs[indx]['lcs'] = []
            
        list_of_papers = sorted(glob.glob(ind+"/20*/*"))
        list_of_seds = sorted(glob.glob(ind+"/20*/*/*-sed*.ecsv"))
        list_of_lcs = sorted(glob.glob(ind+"/20*/*/*-lc*.ecsv"))
        
        for sedfile in list_of_seds:
            sed = Table.read(sedfile, format='ascii.ecsv')
            sed_raw = open(sedfile).read()
            if sed.meta['source_id'] not in VTS_egal_srcs:
                continue

            src_dict = VTS_srcs[sed.meta['source_id']]
            src_dict['seds'].append(sed_raw)

        for lcfile in list_of_lcs:
            lc = Table.read(lcfile, format='ascii.ecsv')
            lc_raw = open(lcfile).read()
            if lc.meta['source_id'] not in VTS_srcs:
                continue

            src_dict = VTS_srcs[sed.meta['source_id']]
            src_dict['lcs'].append(lc_raw)

        with open("vtscat_export_20210215.yaml","w+") as vcf:
            yaml.dump(VTS_srcs,vcf)
                
        self.contents = VTS_srcs
        
    def reset_results(self):
        self.result = dict()
    
    def match_src(self,src,tol=0.5*u.Unit("deg")):
        list_items = self.contents.keys()
        catalog = {}
        
        try:
            src['ra']
        except KeyError:
            src['ra'] = src['pos']['ra']
            src['dec'] = src['pos']['dec']
        
        catalog['dec'] = [self.contents[item]['pos']['dec'] \
                          for item in self.contents \
                          if 'pos' in self.contents[item]]
        catalog['ra']  = [self.contents[item]['pos']['ra'] \
                          for item in self.contents \
                          if 'pos' in self.contents[item]]
        catalog['id']  = [self.contents[item]['source_id'] \
                          for item in self.contents \
                          if 'pos' in self.contents[item]]
        
        try:
            self.matched = self.contents[self.match_src_catalog(catalog,src,tol)]
        except KeyError:
            self.matched = None
            return
        
        if 'seds' not in self.matched:   
            self.matched['seds']=[]
        if 'lcs' not in self.matched:
            self.matched['lcs']=[]
        
    def get_spectralpoints(self):
        self.results = []

        for k, sed in enumerate(self.matched['seds']):
            energyerr = None
            # replace ECSV dump by its Table version
            sed = Table.read(sed, format='ascii.ecsv')
            
            if 'dnde' in sed.colnames:
                flag_vts = np.isfinite(sed['dnde'])
            elif 'e2dnde' in sed.colnames:
                flag_vts = np.isfinite(sed['e2dnde'])
            
            if 'observatory' in sed.colnames:
                flag_vts *= sed['observatory']=='VERITAS'
            
            sed = sed[flag_vts]
            
            # Check if spectrum is valid, otherwise skip
            if len(sed)<=0:
                continue
            
            energies = sed['e_ref'].to("TeV")
            
            if 'dnde' in sed.colnames:
                flux = sed['dnde']
                if 'dnde_err' in sed.keys():
                    fluxerr = sed['dnde_err']
                    logFlux = 2*np.log10(flux/fluxerr.unit)-np.log10((flux/fluxerr.unit + fluxerr/fluxerr.unit))
                    fluxerr = [flux - (10**(logFlux))*fluxerr.unit, fluxerr]
                elif 'dnde_errp' in sed.keys():
                    fluxerr = [sed['dnde_errn'], sed['dnde_errp']]
                else:
                    print('No dnde_err*, check' )            
            elif 'e2dnde' in sed.colnames:
                flux = sed['e2dnde']/(energies**2)
                if 'e2dnde_err' in sed.keys():
                    fluxerr = sed['e2dnde_err']/(energies.to("TeV")**2)
                    logFlux = 2*np.log10(flux/fluxerr.unit) - np.log10((flux/fluxerr.unit + fluxerr/fluxerr.unit))
                    fluxerr = [flux - (10**(logFlux))*fluxerr.unit, fluxerr]
                elif 'e2dnde_errp' in sed.keys():
                    fluxerr = [sed['e2dnde_errn']/(energies.to("TeV")**2), 
                               sed['e2dnde_errp']/(energies.to("TeV")**2)]
                else:
                    print('No e2dnde_err*, check' )
            else:
                print('unexpected error happened, SED does not seem to have a dnde or e2dnde column ... ')
                continue
            
            isnotnan = ~np.isnan(flux) * (flux > 0)
            energies = energies[isnotnan]
            flux = flux[isnotnan]
            fluxerr = np.asarray([fluxerr[0][isnotnan], fluxerr[1][isnotnan]])
            
            # add units
            fluxerr = fluxerr*flux.unit
            
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
            result['name'] = self.matched['common_name']
            result['set'] = k
            result['E'] = energies.to("TeV")
            result['Eerr'] = np.asarray([energyerrn, energyerrp]) * energies.unit
            result['Eerr'] = result['Eerr'].to("TeV")
            result['EF'] = e2flux.to("erg/(cm2*s)")
            result['EFerr'] = e2fluxerr.to("erg/(cm2*s)")
            self.results.append(result)

class FermiCAT(Catalog):
    def __init__(self):
        pass
        
    def set_catalog_urls(self):
        baseurl = "https://fermi.gsfc.nasa.gov/ssc/data/access/lat"
        self.fhl3 = pyfits.open(baseurl+"/3FHL/gll_psch_v13.fit")
        self.fgl4 = pyfits.open(baseurl+"/10yr_catalog/gll_psc_v26.fit")
        self.lac4 = pyfits.open(baseurl+"/4LACDR2/table-4LAC-DR2-h.fits")
        self.lac4_low = pyfits.open(baseurl+"/4LACDR2/table-4LAC-DR2-l.fits")
    
    def set_active_catalog(self,cat):
        self.cat_name = cat
        if cat=='4LAC':
            self.get_spectralpoints_4fgl()
            import numpy.lib.recfunctions
            concat = numpy.lib.recfunctions.stack_arrays(\
                (lac4[1].data,lac4_low[1].data), 
                autoconvert=True, usemask=False)

            self.cat = concat
        elif cat=='4FGL':
            self.get_spectralpoints_4fgl()
            self.cat = self.fgl4[1].data
        elif cat=='3FHL':
            self.get_spectralpoints_3fhl()
            self.cat = self.fhl3[1].data
    
    def get_spectralpoints_3fhl(self):
        self.raw_energies = self.fhl3[4].data
    
    def get_spectralpoints_4fgl(self):
        self.raw_energies = np.asarray([50, 100, 300, 1e3, 3e3, 1e4, 3e4, 3e5])*1e-3
        self.raw_energies = np.transpose([self.raw_energies[:-1], 
                                          self.raw_energies[1:], 
                                          0*self.raw_energies[1:]])
    
    def match_src(self,src,tol=0.5*u.Unit("deg")):
        catalog = dict()
        catalog['dec'] = self.cat['DEJ2000']
        catalog['ra']  = self.cat['RAJ2000']
        catalog['id']  = np.arange(len(catalog['ra']))
        self.matched = self.cat[self.match_src_catalog(catalog,src,tol)]
        if '4FGL' in self.matched['Source_Name']:
            self.lac4_from_fgl4()
        
    def extract_model_parameters(self):
        # Get model parameters from 4FGL/4LAC:
        self.result['pwl'] = dict()
        self.result['lp'] = dict()
        self.result['plec'] = dict()
        # Powerlaw
        self.result['pwl']['PivotE'] = self.matched['Pivot_Energy']*u.Unit("MeV")
        self.result['pwl']['Fnorm'] = self.matched['PL_Flux_Density']*u.Unit("MeV^-1/(s*cm2)")
        self.result['pwl']['EFnorm'] = self.matched['Unc_PL_Flux_Density']*u.Unit("MeV^-1/(s*cm2)")
        self.result['pwl']['Gamma'] = self.matched['PL_Index']
        self.result['pwl']['EGamma'] = self.matched['Unc_PL_Index']
        # LogParabola
        self.result['lp']['PivotE'] = self.matched['Pivot_Energy']*u.Unit("MeV")
        self.result['lp']['Fnorm'] = self.matched['LP_Flux_Density']*u.Unit("MeV^-1/(s*cm2)")
        self.result['lp']['EFnorm'] = self.matched['Unc_LP_Flux_Density']*u.Unit("MeV^-1/(s*cm2)")
        self.result['lp']['Gamma'] = self.matched['LP_Index']
        self.result['lp']['EGamma'] = self.matched['Unc_LP_Index']
        self.result['lp']['Beta'] = self.matched['LP_beta']
        self.result['lp']['EBeta'] = self.matched['Unc_LP_beta']
        # Powerlaw with exponential cutoff
        self.result['plec']['PivotE'] = self.matched['Pivot_Energy']*u.Unit("MeV")
        self.result['plec']['Fnorm'] = self.matched['PLEC_Flux_Density']*u.Unit("MeV^-1/(s*cm2)")
        self.result['plec']['EFnorm'] = self.matched['Unc_PLEC_Flux_Density']*u.Unit("MeV^-1/(s*cm2)")
        self.result['plec']['Gamma'] = self.matched['PLEC_Index']
        self.result['plec']['EGamma'] = self.matched['Unc_PLEC_Index']
        self.result['plec']['ExpF'] = self.matched['PLEC_Expfactor']
        self.result['plec']['EExpF'] = self.matched['Unc_PLEC_Expfactor']
        self.result['plec']['ExpI'] = self.matched['PLEC_Exp_Index']
        self.result['plec']['EExpI'] = self.matched['Unc_PLEC_Exp_Index']
        # Preferred model (4FGL)
        if self.matched['SpectrumType']=="PowerLaw":
            self.result['spectrumtype'] = "pwl"
        elif self.matched['SpectrumType']=="LogParabola":
            self.result['spectrumtype'] = "lp"
        elif self.matched['SpectrumType']=="PLExpCutoff":
            self.result['spectrumtype'] = "plec"

    def fill_spectrum(self):
        Eraw = self.raw_energies
        E = np.asarray([10 ** np.mean(np.log10(k[:-1])) for k in Eraw])
        Eerrn = np.asarray([E[k] - re[0] for k, re in enumerate(Eraw)])
        Eerrp = np.asarray([re[1] - E[k] for k, re in enumerate(Eraw)])
        
        # Add units
        F      = self.matched['Flux_Band']* u.Unit("cm**-2 * s**-1")
        Ferr   = np.transpose(np.abs(self.matched['Unc_Flux_Band']))*u.Unit("cm**-2 * s**-1")
        E      = E * u.Unit("GeV")
        SED    = (F * E).to("erg/(cm2*s)")
        SEDerr = (Ferr * E).to("erg/(cm2*s)")
        
        isnotnan = (~np.isnan(F))*\
                   (~np.isnan(Ferr[0])) * (F > 0)
        
        
        
        self.result['E'] = E[isnotnan].to("TeV")
        Eerr = np.asarray([Eerrn[isnotnan],Eerrp[isnotnan]])*u.Unit("GeV")
        self.result['Eerr'] = Eerr.to("TeV")
        self.result['EF'] = SED[isnotnan]
        self.result['EFerr'] = np.asarray(
            [SEDerr[0][isnotnan],SEDerr[1][isnotnan]])
        self.result['EFerr'] = self.result['EFerr']*SED.unit
        #return(self.result)
        
    def lac4_from_fgl4(self):
        srcnamefilt = self.lac4[1].data['Source_Name'] == self.matched['Source_Name']
        if np.sum(srcnamefilt) == 0:
            srcnamefilt = self.lac4_low[1].data['Source_Name'] == self.matched['Source_Name']
            if np.sum(srcnamefilt) == 1:
                self.matched_lac4 = self.lac4_low[1].data[srcnamefilt]
            else:
                self.matched_lac4 = None
        else:
            self.matched_lac4 = self.lac4[1].data[srcnamefilt]
        
        return(self.matched_lac4[0])
        
    def get_properties(self):
        self.result['name']  = self.matched['Source_Name']
        if '3FHL' in self.result['name']:
            self.result['redshift'] = self.matched['REDSHIFT']
            if self.matched['NuPeak_obs']<1e14:
                self.result['class'] = 'LSP'
            elif self.matched['NuPeak_obs']<1e15:
                self.result['class'] = 'ISP'
            elif self.matched['NuPeak_obs']<1e17:
                self.result['class'] = 'HSP'
            elif self.matched['NuPeak_obs']>1e18:
                self.result['class'] = 'EHSP'
            else:
                self.result['class'] = ''
        elif '4FGL' in self.result['name']: # 
            self.result['redshift'] = self.matched_lac4['Redshift'][0]
            if self.matched_lac4['nu_syn']<1e14:
                self.result['class'] = 'LSP'
            elif self.matched_lac4['nu_syn']<1e15:
                self.result['class'] = 'ISP'
            elif self.matched_lac4['nu_syn']<1e17:
                self.result['class'] = 'HSP'
            elif self.matched_lac4['nu_syn']>1e18:
                self.result['class'] = 'EHSP'
            else:
                self.result['class'] = ''