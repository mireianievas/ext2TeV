#!/bin/env python

import yaml
import glob
import os
import matplotlib.pyplot as plt

from ext2vhe.residuals import Residuals
from ext2vhe import PACKAGE_DATA
from ext2vhe.catalog import build_gammacat
from ext2vhe.source import Source

res = Residuals()
source = Source()
source.set_residuals(res)
source.set_ebl_interpolator()

background_dcs_file = os.path.join(PACKAGE_DATA, 'background_dcs/background_dcs.txt')


dataset = None
try:
    print('Trying to read latest available gammacat catalog:')
    gammacats = glob.glob("{}/gammacat_export_*.yaml".format(PACKAGE_DATA))
    gammacats.sort()
    print(gammacats)
    # Open latest gammacat export:
    with open(gammacats[-1]) as gcf:
        dataset = yaml.load(gcf)
except:
    #raise
    print('building catalog')
    dataset = build_gammacat(gammacat_path="../Catalogs/gamma-cat/output")

if dataset is None:
    print("Not able to generate the dataset. Exitting...")
    exit(-1)

agn_classes = ['fsrq', 'hbl', 'ibl', 'lbl', 'fri']

# Go through all sources belonging to any of the agn classes listed above, and compare their extrapolation
# to VHE from LAT catalogs with actual VHE observations.
for k, src in enumerate(dataset):
    # if 'PKS 1424+240' not  in datasets[src]['common_name']: continue

    if dataset[src]['classes'] not in agn_classes:
        continue

    ### extract the sed
    if 'sed' not in dataset[src].keys():
        print('no sed')
        continue

    src_data = dataset[src]
    source.broadband_sed(dataset=src_data, model='plec', figure_path='results/individual_fits')


def compare_class_extrapolations(agn_class, source):
    fig, spls = plt.subplots(nrows=3, ncols=1, sharex=True, sharey=True,
                             squeeze=True, subplot_kw=None, gridspec_kw=None, figsize=(6, 4), dpi=150)
    for k, mod in enumerate(source.residuals.xlogval[agn_class]):
        spls[k].scatter(
            source.residuals.xlogval[agn_class][mod],
            source.residuals.ylogres[agn_class][mod],
            label=mod,
            s=2,
            marker=['o', 's', 'D'][k],
            facecolor='None',
            color='C{}'.format(k),
            alpha=0.3,
        )

        spls[k].axhline(0, lw=0.2, ls='dashed', color='black', zorder=-10)
        spls[k].set_xlim([-1, 1.5])
        spls[k].set_ylim([-40, 40])
        spls[k].legend()
        spls[k].grid(ls='dotted', lw=0.33)

    spls[0].set_title(agn_class)
    spls[1].set_ylabel('Residues (log-fluxes)')
    spls[-1].set_xlabel('log E [TeV]')
    return fig


for agn_class in agn_classes:
    fig = compare_class_extrapolations(agn_class, source)
    plt.savefig("results/compare_extrapolations_class_{}.png".format(agn_class),
                dpi=200, bbox_inches='tight')
