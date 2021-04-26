import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt 

from ext2TeV.catalog import GammaCAT, VTSCAT, FermiCAT
from ext2TeV.residuals import Residuals
from ext2TeV.fitter import Fitter
from ext2TeV.eblmodel import Eblmodel

from astropy.io import fits as pyfits

# Load LAT and GammaCAT catalogs 

# GammaCAT (VHE collection of figures)
gammacat = GammaCAT()
gammacat.read("ext2TeV/data/gammacat_export_20210215.yaml")

# Fermi catalogs (4FGL/4LAC and 3FHL)
fermicat = FermiCAT()
fermicat.set_catalog_urls()
fermicat3fhl = FermiCAT()
fermicat3fhl.set_catalog_urls()

# Create a residuals object (we will it with content later)
residuals = Residuals()

# Load EBL
ebl = Eblmodel()
ebl.load_ebl_dominguez()
ebl.set_interpolator()

# Look forextragalactic sources
extragalactic = []
for k in gammacat.contents.keys():
    src = gammacat.contents[k]
    if np.asarray([src['classes'],]).flatten()[0] in ['hbl','lbl', 'ibl', 'fsrq', 'agn']:
        extragalactic.append(k)

vtscat = VTSCAT()
vtscat.read("ext2TeV/data/vtscat_reduced_v0.1a_20210423.yaml")

extragalactic_in_vts = []
for indx,element in enumerate(vtscat.contents):
    try:    assert(vtscat.contents[element]['where']=='egal')
    except: continue
    else:   extragalactic_in_vts.append(element)

# Loop over each extragalactic object
for indx,element in enumerate(extragalactic+extragalactic_in_vts):
    # Load GammaCAT spectra
    myfitter = Fitter()

    if indx < len(extragalactic):
        # srcs in gammacat
        src = gammacat.contents[element]
        gammacat.reset_results()
        gammacat.match_src(src)
        gammacat.get_spectralpoints()
        myfitter.add_vhe(gammacat.results)
        vtscat.match_src(src)
        if vtscat.matched != None:
            # srcs in gammacat and vtscat
            vtscat.get_spectralpoints()
            myfitter.add_vhe(vtscat.results)
    else:
        # srcs in vtscat
        src = vtscat.contents[element]
        vtscat.reset_results()
        vtscat.match_src(src)
        vtscat.get_spectralpoints()
        myfitter.add_vhe(vtscat.results)
        
        gammacat.match_src(src)
        if gammacat.matched != None:
            # if it is already in gammacat skip it, already processed.
            continue
    
    
    if myfitter.vhe == []: continue
    
    # Locate 3FHL counterpart (if any) and its spectrum 
    fermicat3fhl.reset_results()
    fermicat3fhl.reset_spectral_points()
    fermicat3fhl.set_active_catalog("3FHL")
    fermicat3fhl.match_src(src)
    fermicat3fhl.fill_spectrum()
    fermicat3fhl.get_properties()
    
    # Locate 4LAC counterpart (if any) and its spectrum / Likelihood models 
    fermicat.reset_results()
    fermicat.reset_spectral_points()
    fermicat.set_active_catalog("4FGL")
    fermicat.match_src(src)
    fermicat.fill_spectrum()
    fermicat.extract_model_parameters()
    fermicat.get_properties()
    
    # Create a Fitter object, add the previous spectral points and compute residuals
    mpl.rcParams['text.usetex'] = False
    #myfitter = Fitter()

    # Add LAT and VHE data points
    myfitter.add_fermi(fermicat.result)
    myfitter.add_fermi2(fermicat3fhl.result)
    #myfitter.add_vhe(gammacat.results)
    
    # Combine spectra
    myfitter.build_broadband_spectra()
    myfitter.set_ebl(ebl)

    # Set redshift and spectral class (to test CTA's proposed extension to TeV)
    myfitter.set_model_redshift(fermicat.result['redshift'])
    myfitter.set_model_spectralclass(fermicat.result['class'])

    # Create an empty figurure
    myfitter.create_figure()

    # Get lowest and elevated states
    myfitter.get_state()

    # Plot data
    myfitter.plot_fermi_data()
    myfitter.plot_vhe_data()

    # Compute LAT-predefined models and CTA's adopted solution, then fill residuals.
    myfitter.predefined_models(residuals)

    # Styling and save figure
    myfitter.finish_figure()
    myfitter.savefig()

# Plot residuals
residuals.create_figure()
residuals.plot_panels()
residuals.savefig()