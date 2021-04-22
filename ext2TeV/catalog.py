#import gammacat
#import ruamel.yaml as yaml
import json
from astropy.io.misc import yaml
import os
from datetime import datetime
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
from astropy.table import Table
from astropy.io import fits as pyfits

from ext2TeV import PACKAGE_DATA
from ext2TeV.utils import read_yaml


def build_gammacat(gammacat_path=""):
    """
    Export gammacat into a yaml file, to be used by the etx2vhe library

    Parameters
    ----------

    gammacat_path: str
        Path containing the "gammacat.yaml" and "gammacat-datasets.json" files.
    """
    dt = datetime.now()

    os.chdir(gammacat_path)
    # p = pyfits.open("gammacat.fits.gz")

    with open("gammacat.yaml") as gammacat:
        temp = yaml.load(gammacat.read())
        gammacat = dict()
        for k, source_id in enumerate(temp):
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
                    content = open(item['location']).read()  # Table.read(item['location'])
                except:
                    print('Skipping {0}'.format(item['location']))

                if 'sed' in datasets[item['source_id']]:
                    datasets[item['source_id']]['sed'].append(content)
                else:
                    datasets[item['source_id']]['sed'] = [content]

    output_file = os.path.join(PACKAGE_DATA, "data/gammacat_export_{}{}{}.yaml".format(dt.year, dt.month, dt.day))
    with open(output_file, "w+") as gcf:
        yaml.dump(datasets, gcf)

    return datasets
