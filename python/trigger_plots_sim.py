from __future__ import division



import hists


from MODPlot import *


import trigger_parse_sim





from subprocess import call

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator
from matplotlib.ticker import FixedLocator
from matplotlib.ticker import LogLocator
from matplotlib.ticker import FormatStrFormatter

from sets import Set

import time as time
import copy


import sys
import math
from collections import defaultdict

# matplotlib
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
from matplotlib.legend_handler import HandlerLine2D
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm


# RootPy
from rootpy.plotting import Hist, HistStack, Legend
import rootpy.plotting.root2matplotlib as rplt
from rootpy.plotting import Hist2D


# Stuff for calculating areas.
from scipy.integrate import simps
from scipy import interpolate
from scipy import optimize

from numpy import trapz


from matplotlib import gridspec

import matplotlib.ticker as mtick

from matplotlib.cbook import get_sample_data
from matplotlib._png import read_png
from matplotlib.backends.backend_pdf import PdfPages

from matplotlib.offsetbox import TextArea, DrawingArea, OffsetImage, AnnotationBbox, AnchoredOffsetbox, HPacker

from mpl_toolkits.axes_grid.anchored_artists import AnchoredDrawingArea

from scipy.stats import norm
from scipy.stats import gamma
from scipy import arange, array, exp

from scipy.stats import binned_statistic

import rootpy.plotting.views




mpl.rcParams['axes.linewidth'] = 5.0 #set the value globally
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
mpl.rcParams['text.latex.preamble'] = [r'\boldmath']

plt.rc('font', family='serif', size=43)




default_dir = "plots2/Version1/"


logo_location = "/home/preksha/Documents/mengproject/MODAnalyzer/mod_logo.png"
logo_text = "Preliminary (16\%)"




parsed_linear = trigger_parse_sim.load_root_files_to_hist()
print(parsed_linear)


def logo_box(x, y):
	
	logo_offset_image = OffsetImage(read_png(get_sample_data(logo_location, asfileobj=False)), zoom=0.31, resample=1, dpi_cor=1)
	text_box = TextArea(logo_text, textprops=dict(color='#444444', fontsize=70, weight='bold'))

	logo_and_text_box = HPacker(children=[logo_offset_image, text_box], align="center", pad=0, sep=25)

	anchored_box = AnchoredOffsetbox(loc=2, child=logo_and_text_box, pad=0.8, frameon=False, borderpad=0., bbox_to_anchor=[x, y], bbox_transform = plt.gcf().transFigure)

	return anchored_box



def normalize_hist(hist):
	bin_width = (hist.upperbound() - hist.lowerbound()) / hist.nbins()
	
	if hist.GetSumOfWeights() != 0.0:
		#hist.Scale(1.0 / ( hist.GetSumOfWeights() * bin_width ))
                hist.Scale(1.0 / ( bin_width ))

	return hist



def trigger_turn_on_curves_half2():

	mod_hists = parsed_linear[0]

	print(mod_hists.keys())

        colors = ["#999999", 'red', "#f781bf",  "#a65628", "#99cc00",
                  "#ff7f00", "#984ea3", "#4daf4a", "#377eb8",
                  "#e41a1c",'#8da0cb',  '#66c2a5', '#fc8d62'][7:]
	labels = [ '$ p_T^{\mathrm{Gen}} \in$ [15,30] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [30,50] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [50,80] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [80,120] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [120,170] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [170,300] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [300,470] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [470,600] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [600,800] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [800,1000] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [1000,1400] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [1400,1800] GeV',
                   '$ p_T^{\mathrm{Gen}} \ge$ 1800 GeV'][7:]
	hist_labels = ['QCD_Pt-15to30_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-30to50_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-50to80_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-80to120_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-120to170_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-170to300_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-300to470_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-470to600_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-600to800_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-800to1000_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-1800_TuneZ2_7TeV_pythia6'][7:]
	colors.reverse()
	labels.reverse()
	hist_labels.reverse()


        print(len(colors), len(labels), len(hist_labels))

	for i in range( len(hist_labels)):
		
		hist = normalize_hist( mod_hists[hist_labels[i]].hist() )
		print(hist_labels[i])
		#hist = mod_hists[hist_labels[i]].hist()

		n_bins = hist.nbins()

		#print(n_bins)

		# Loop through those bins.
		for j in range(n_bins):
			current_bin_x = hist.GetBinCenter(j)
			current_bin_y = hist.GetBinContent(j),
			#print(current_bin_x, current_bin_y)

			"""
			if current_bin_x < lower_pTs[i]:
				hist.SetBinContent(j, 0.)
                        """

		hist.SetColor(colors[i])
		hist.SetTitle(labels[i])

		rplt.errorbar(hist, zorder=range(len(hist_labels))[len(hist_labels) - i - 1],
                              emptybins=False, xerr=1, yerr=1, ls='None', marker='o',
                              markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=1.0)
		print(i)


	extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
	handles = [extra]
	labels = ["Pythia 6 Tune Z2 (Simulated) \n" + "AK5; $\left| \eta \\right| < 2.4$"]

	
	info_legend = plt.gca().legend(handles, labels, loc=7, frameon=0, borderpad=0.1, fontsize=59, bbox_to_anchor=[0.52, 0.9])

	plt.gca().add_artist(info_legend)


	plt.gca().set_xlabel("Hardest Jet $p_T$ [GeV]", fontsize=80, labelpad=50)
	plt.gca().set_ylabel("Cross Section [pb] ", fontsize=80, labelpad=25)

	plt.gca().add_artist(logo_box(0.13,0.985))

	plt.gca().xaxis.set_minor_locator(MultipleLocator(50))
	plt.gca().set_yscale('log')

	handles, labels = plt.gca().get_legend_handles_labels()

	legend = plt.legend(handles[::-1], labels[::-1],  frameon=0, fontsize=55, bbox_to_anchor=[0.99, 0.99])

	ax = plt.gca().add_artist(legend)

	outside_text = plt.gca().legend( [extra], ["CMS 2011 Open Data"], frameon=0, borderpad=0, fontsize=55, bbox_to_anchor=(1.0, 1.005), loc='lower right')
	plt.gca().add_artist(outside_text)


	plt.autoscale()
	plt.gca().set_ylim(1e-10, 1e2)
	plt.gca().set_xlim(0, 3000)
	
	plt.tick_params(which='major', width=5, length=25, labelsize=70)
	plt.tick_params(which='minor', width=3, length=15)

	plt.gcf().set_size_inches(31, 31, forward=1)

	plt.tight_layout()

	plt.savefig(default_dir + "qcd_all_datasets_hardest_pT_half2_mod.pdf")

	plt.clf()


def trigger_turn_on_curves_half1():

	mod_hists = parsed_linear[0]

	print(mod_hists.keys())

        colors = ["#999999", 'red', "#f781bf",  "#a65628", "#99cc00",
                  "#ff7f00", "#984ea3", "#4daf4a", "#377eb8",
                  "#e41a1c",'#8da0cb',  '#66c2a5', '#fc8d62'][:7]
	labels = [ '$ p_T^{\mathrm{Gen}} \in$ [15,30] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [30,50] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [50,80] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [80,120] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [120,170] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [170,300] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [300,470] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [470,600] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [600,800] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [800,1000] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [1000,1400] GeV',
                   '$ p_T^{\mathrm{Gen}} \in$ [1400,1800] GeV',
                   '$ p_T^{\mathrm{Gen}} \ge$ 1800 GeV'][:7]
	hist_labels = ['QCD_Pt-15to30_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-30to50_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-50to80_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-80to120_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-120to170_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-170to300_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-300to470_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-470to600_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-600to800_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-800to1000_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6',
                       'QCD_Pt-1800_TuneZ2_7TeV_pythia6'][:7]
	colors.reverse()
	labels.reverse()
	hist_labels.reverse()


        print(len(colors), len(labels), len(hist_labels))

	for i in range( len(hist_labels)):
		
		hist = normalize_hist( mod_hists[hist_labels[i]].hist() )
		print(hist_labels[i])
		#hist = mod_hists[hist_labels[i]].hist()

		n_bins = hist.nbins()

		#print(n_bins)

		# Loop through those bins.
		for j in range(n_bins):
			current_bin_x = hist.GetBinCenter(j)
			current_bin_y = hist.GetBinContent(j),
			#print(current_bin_x, current_bin_y)

			"""
			if current_bin_x < lower_pTs[i]:
				hist.SetBinContent(j, 0.)
                        """

		hist.SetColor(colors[i])
		hist.SetTitle(labels[i])

		rplt.errorbar(hist, zorder=range(len(hist_labels))[len(hist_labels) - i - 1],
                              emptybins=False, xerr=1, yerr=1, ls='None', marker='o',
                              markersize=10, pickradius=8, capthick=5, capsize=8, elinewidth=5, alpha=1.0)
		print(i)


	extra = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
	handles = [extra]
	labels = ["Pythia 6 Tune Z2 (Simulated) \n" + "AK5; $\left| \eta \\right| < 2.4$"]

	
	info_legend = plt.gca().legend(handles, labels, loc=7, frameon=0, borderpad=0.1, fontsize=59, bbox_to_anchor=[0.52, 0.9])

	plt.gca().add_artist(info_legend)


	plt.gca().set_xlabel("Hardest Jet $p_T$ [GeV]", fontsize=80, labelpad=50)
	plt.gca().set_ylabel("Cross Section [pb] ", fontsize=80, labelpad=25)

	plt.gca().add_artist(logo_box(0.12, 0.985))

	plt.gca().xaxis.set_minor_locator(MultipleLocator(10))
	plt.gca().set_yscale('log')

	
	handles, labels = plt.gca().get_legend_handles_labels()

	legend = plt.legend(handles[::-1], labels[::-1],  frameon=0, fontsize=55, bbox_to_anchor=[0.99, 0.99])

	ax = plt.gca().add_artist(legend)

	outside_text = plt.gca().legend( [extra], ["CMS 2011 Open Data"], frameon=0, borderpad=0, fontsize=55, bbox_to_anchor=(1.0, 1.005), loc='lower right')
	plt.gca().add_artist(outside_text)


	plt.autoscale()
	plt.gca().set_ylim(1e-3, 1e9)
	plt.gca().set_xlim(0, 700)


	
	plt.tick_params(which='major', width=5, length=25, labelsize=70)
	plt.tick_params(which='minor', width=3, length=15)

	plt.gcf().set_size_inches(30, 30, forward=1)

	plt.tight_layout()

	plt.savefig(default_dir + "qcd_all_datasets_hardest_pT_half1.pdf")

	plt.clf()

start = time.time()

trigger_turn_on_curves_half1()
#trigger_turn_on_curves_half2()

end = time.time()

print "Finished all plotting in {} seconds.".format(end - start)
