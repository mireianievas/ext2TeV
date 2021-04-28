import numpy as np
import matplotlib.pyplot as plt

class Residuals(object):
    def __init__(self):
        # dict with log of values of X and Y-residues for different type of sources.
        self.ids = {}
        self.xlogval = {}
        self.ylogval = {}
        self.ylogmod = {}
        self.ylogres = {}
        self.ylogerr = {}

    def reset_cls(self,cls):
        self.ids[cls] = {}
        self.xlogval[cls] = {}
        self.ylogval[cls] = {}
        self.ylogmod[cls] = {}
        self.ylogres[cls] = {}
        self.ylogerr[cls] = {}

    def reset_cls_model(self,cls,model):
        if cls not in self.xlogval:
            self.reset_cls(cls)

        self.ids[cls][model] = []
        self.xlogval[cls][model] = []
        self.ylogval[cls][model] = []
        self.ylogmod[cls][model] = []
        self.ylogres[cls][model] = []
        self.ylogerr[cls][model] = []

    def add_value(self, point):
        ids = point['ids']
        cls = point['srccls']
        model = point['model']
        xlogval = point['xlogval']
        ylogval = point['ylogval']
        ylogerr = point['ylogerr']
        ylogmod = point['ylogmod']
        # ylogres = point['ylogres']

        if cls not in self.xlogval:
            self.reset_cls(cls)

        if model not in self.xlogval[cls]:
            self.reset_cls_model(cls,model)

        self.ids[cls][model].append(ids)
        self.xlogval[cls][model].append(xlogval)
        self.ylogval[cls][model].append(ylogval)
        self.ylogmod[cls][model].append(ylogmod)
        self.ylogres[cls][model].append((ylogval -ylogmod ) /ylogerr)
        
    

    def create_figure(self):
        self.fig = plt.figure(figsize=(8,9),dpi=120)

    def plot_panels(self):
        _num_ = 0
        energy_bins = np.arange(-1.2,1.0,0.3)
        
        for i,src_class in enumerate(list(self.xlogval.keys())):

            for j,model in enumerate(list(self.xlogval[src_class].keys())):

                # extract residual values
                ids    = self.ids[src_class][model]
                E_log  = self.xlogval[src_class][model]
                EP_log = self.ylogval[src_class][model]
                EPerr_log = self.ylogval[src_class][model]
                EPmodel_log = self.ylogmod[src_class][model]
                EPresid_log = self.ylogres[src_class][model]
                
                # add subplot on its position
                _num_ += 1 
                plo = self.fig.add_subplot(len(list(self.xlogval.keys())),
                                           len(list(self.xlogval[src_class].keys())),
                                           _num_)

                
                # iterate over energy bins
                for ibin,_ in enumerate(energy_bins[:-1]):
                    emin = energy_bins[ibin]
                    emax = energy_bins[ibin+1]

                    filt = (E_log>=emin)*(E_log<emax)*np.isfinite(EPresid_log)

                    if np.sum(filt)>0:
                        ave = np.mean(np.asarray(EPresid_log)[filt])
                        std = np.std(np.asarray(EPresid_log)[filt])

                        plo.bar(
                            x=(emin+emax)/2., 
                            height=ave,
                            yerr=std,
                            width=(emax-emin)*0.95,
                            align='center', 
                            alpha=0.5, 
                            ecolor=(0,0,0,0.5), 
                            capsize=4,
                            ls='None',
                            #marker='o',
                            color='C{}'.format(j),
                            #ms=1.5,
                            #mew=1.5,
                            #lw=0.5,
                            #mfc='white',
                            zorder=5,
                        )

                plo.errorbar(
                    x=E_log,
                    y=EPresid_log,
                    ls='None',
                    marker='o',
                    ms=1.5,
                    mew=0.1,
                    mfc='C{}'.format(j),
                    color='black',
                    alpha=0.5,
                    zorder=-10,
                )


                # Add number of spectra
                props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
                plo.text(0.8,0.1,len(np.unique(sorted(ids))), 
                         transform=plo.transAxes, fontsize='large', bbox=props)

                plo.axhline(0,color='black',ls='solid',lw=0.33)
                plo.set_xlim([-1.3,1])

                if j==0:
                    plo.set_ylabel('$log-Residuals \ (\sigma)$',fontsize='small')

                if i==len(list(self.xlogval.keys()))-1:
                    plo.set_xlabel('$\mathrm{\log_{10} (Energy/TeV)}$')

                    
                title = src_class if src_class != '' else 'other'
                if 'EBL' not in title.upper():
                    title = title.upper()+'+EBL'
                
                plo.set_title('{0}/{1}'.format(src_class if src_class != '' else 'other',
                                               model),fontsize='medium')

                plo.set_ylim(-np.max(np.abs(plo.get_ylim())),np.max(np.abs(plo.get_ylim())))
                plo.grid(lw=0.5,ls='dotted',alpha=0.5,color='black')

        self.fig.tight_layout(pad=0.0)
        
    def savefig(self,filename=None):
        if filename==None:
            filename = "results/Summary_residuals_ext2TeV.pdf"
        
        self.fig.savefig(filename, bbox_inches='tight')
