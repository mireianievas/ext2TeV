import yaml
import glob
import os

from ext2TeV.residuals import Residuals
from ext2TeV import PACKAGE_DATA
from ext2TeV.catalog import build_gammacat
from ext2TeV.source import Source

res  = Residuals()
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
