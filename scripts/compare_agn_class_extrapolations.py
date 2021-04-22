#!/bin/env python

import yaml
import glob
import os
import matplotlib.pyplot as plt
import numpy as np

from ext2TeV import PACKAGE_DATA
from ext2TeV.residuals import Residuals
from ext2TeV.catalog import build_gammacat
from ext2TeV.source import Source

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
	raise
	print('building catalog')
	dataset = build_gammacat(gammacat_path="../Catalogs/gamma-cat/output")

if dataset is None:
	print("Not able to generate the dataset. Exitting...")
	exit(-1)

agn_classes = ['fsrq', 'hbl', 'ibl', 'lbl', 'fri']


def compare_class_extrapolations(agn_class, source, lbl_tail = ''):
	fig, spls = plt.subplots(nrows=3, ncols=1, sharex=True, sharey=True,
							 squeeze=True, subplot_kw=None, gridspec_kw=None, figsize=(6, 4), dpi=150)

	with open('results/data_{}_{}.yml'.format(agn_class,lbl_tail), 'w') as outfile:
		yaml.dump(source.residuals, outfile, default_flow_style=True)
	
	for k, mod in enumerate(source.residuals.xlogval[agn_class]):
		
		spls[k].scatter(
			source.residuals.xlogval[agn_class][mod],
			source.residuals.ylogres[agn_class][mod],
			label=mod if mod not in lbl_tail else lbl_tail,
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


# Go through all sources belonging to any of the agn classes listed above, and compare their extrapolation
# to VHE from LAT catalogs with actual VHE observations.


model_pool = ['pwl','lp']
model_pool = []

for ExpIndex in np.linspace(0.5,1.0,13):
	model_pool.append('plec,ExpI={0:.3f}'.format(ExpIndex))

for agn_class in agn_classes:
	for model_sel in model_pool:
		source.residuals.reset_cls_model(agn_class,model_sel.split(",")[0])
		
		print("-----> Using model {} now ~~~~~~~".format(model_sel))
		# build the broadband sed
		for k, src in enumerate(dataset):
			### extract the sed
			if 'sed' not in dataset[src].keys():
				print('no sed')
				continue
			
			source.broadband_sed(dataset=dataset[src], 
				model=model_sel, 
				figure_path='results/individual_fits')

		model_lbl = model_sel.replace(",","_")
		if model_lbl in ['pwl','lp']:
			continue

		fig = compare_class_extrapolations(agn_class, source, lbl_tail='_{}'.format(model_lbl))
		plt.savefig("results/compare_extrapolations_class_{}_{}.png".format(agn_class,model_lbl),
			 dpi=200, bbox_inches='tight')

