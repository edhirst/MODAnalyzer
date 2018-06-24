from __future__ import division
import numpy as np
import copy
import math

from collections import defaultdict
# RootPy
from rootpy.plotting import Hist, Hist2D


upper_left = (-0.075, 0.97)
upper_right = (0.80, 0.97)
# mid_right = (0.73, 0.6)  # => Use for JEC corrected pT.
mid_right = (0.76, 0.6)  # => Regular pT spectrum.
# mid_right = (0.80, 0.6)


mid_plus_epsilon_right = (0.80, 0.5)

upper_left_for_log_x = (-0.08, 0.87)
big_5_upper_right = {'85': (0.87, 0.97), '150': (
    0.87, 0.97), '250': (0.94, 0.97)}


big_5_mid_right = (0.97, 0.53)


class MODHist:

    def __init__(self, hist, conditions=[], use_prescale=True, x_label="", y_label="", x_scale='linear', y_scale='linear', mark_regions=[], x_range=(0, -1), y_range=(0, -1), additional_text=[], legend_location=('upper right', (1., 1.)), axes_label_pi=False):
        # hist should be a RootPy Hist object (TH1).
        # conditions should be a list of lambda functions that need to be
        # satisfied by 'value' for this hist to be filled.

        self._hist = hist
        self._conditions = conditions
        self._use_prescale = use_prescale

        self._x_label = x_label
        self._y_label = y_label

        self._x_scale = x_scale
        self._y_scale = y_scale

        self._mark_regions = mark_regions

        self._x_range = x_range
        self._y_range = y_range

        self._additional_text = additional_text

        self._axes_label_pi = axes_label_pi

        self._legend_location = legend_location

    def replace_hist(self, hist):
        self._hist = hist

    def hist(self):
        return self._hist

    def conditions(self):
        return self._conditions

    def use_prescale(self):
        return self._use_prescale

    def x_label(self):
        return self._x_label

    def y_label(self):
        return self._y_label

    def x_scale(self):
        return self._x_scale

    def y_scale(self):
        return self._y_scale

    def x_range(self):
        return self._x_range

    def y_range(self):
        return self._y_range

    def additional_text(self):
        return self._additional_text

    def axes_label_pi(self):
        return self._axes_label_pi

    def mark_regions(self):
        return self._mark_regions

    def legend_location(self):
        return self._legend_location


def multi_page_log_plot_hist_templates():
    all_hists = {}

    hardest_pT_hist = Hist(np.logspace(
        math.log(float(5), math.e), math.log(1505, math.e), 150 + 1, base=np.e))
    pT_D_hist = Hist(np.logspace(math.log(float(0.1), math.e),
                                 math.log(1, math.e), 75 + 1, base=np.e))

    frac_pT_loss_hist = Hist(np.logspace(
        math.log(float(0.001), math.e), math.log(1.0, math.e), 50 + 1, base=np.e))
    softkill_pT_loss_hist = Hist(np.logspace(
        math.log(float(0.0001), math.e), math.log(0.1, math.e), 50 + 1, base=np.e))

    lha_hist = Hist(np.logspace(math.log(float(0.1), math.e),
                                math.log(1.0, math.e), 50 + 1, base=np.e))
    width_hist = Hist(np.logspace(math.log(float(0.01), math.e),
                                  math.log(1.0, math.e), 50 + 1, base=np.e))
    thrust_hist = Hist(np.logspace(math.log(float(0.001), math.e),
                                   math.log(1.0, math.e), 50 + 1, base=np.e))

    zg_10_hist = Hist(np.logspace(math.log(float((0.1**2) / 0.5), math.e),
                                  math.log((0.5**2) / 0.1, math.e), (50 * 3 + 1), base=np.e))
    rg_10_hist = Hist(np.logspace(math.log(float(0.01), math.e),
                                  math.log(1.0, math.e), 50 + 1, base=np.e))
    e1_10_hist = Hist(np.logspace(math.log(float(0.001), math.e),
                                  math.log(1.0, math.e), 50 + 1, base=np.e))
    e05_10_hist = Hist(np.logspace(math.log(float(0.01), math.e),
                                   math.log(1.0, math.e), 50 + 1, base=np.e))
    e2_10_hist = Hist(np.logspace(math.log(float(0.0001), math.e),
                                  math.log(1.0, math.e), 50 + 1, base=np.e))

    zg_05_hist = Hist(np.logspace(math.log(float((0.05**2) / 0.5), math.e),
                                  math.log((0.5**2) / 0.05, math.e), (50 * 3 + 1), base=np.e))
    rg_05_hist = Hist(np.logspace(math.log(float(0.005), math.e),
                                  math.log(1.0, math.e), 50 + 1, base=np.e))
    e1_05_hist = Hist(np.logspace(math.log(float(0.0005), math.e),
                                  math.log(1.0, math.e), 50 + 1, base=np.e))
    e05_05_hist = Hist(np.logspace(math.log(float(0.005), math.e),
                                   math.log(1.0, math.e), 50 + 1, base=np.e))
    e2_05_hist = Hist(np.logspace(math.log(float(0.00005), math.e),
                                  math.log(1.0, math.e), 50 + 1, base=np.e))

    zg_20_hist = Hist(np.logspace(math.log(float((0.2**2) / 0.5), math.e),
                                  math.log((0.5**2) / 0.2, math.e), (50 * 3 + 1), base=np.e))
    rg_20_hist = Hist(np.logspace(math.log(float(0.01), math.e),
                                  math.log(1.0, math.e), 50 + 1, base=np.e))
    e1_20_hist = Hist(np.logspace(math.log(float(0.001), math.e),
                                  math.log(1.0, math.e), 50 + 1, base=np.e))
    e05_20_hist = Hist(np.logspace(math.log(float(0.01), math.e),
                                   math.log(1.0, math.e), 50 + 1, base=np.e))
    e2_20_hist = Hist(np.logspace(math.log(float(0.0001), math.e),
                                  math.log(1.0, math.e), 50 + 1, base=np.e))

    pT_boundaries = [85, 115, 150, 200, 250]
    # eta_boundaries = [(0.0, 1.9), (0.0, 2.4), (1.9, 2.4)]
    eta_boundaries = [(-2.4, 2.4)]

    all_hists['frac_pT_loss'] = []
    all_hists['softkill_pT_loss'] = []

    all_hists['LHA_pre_SD'], all_hists['track_LHA_pre_SD'] = [], []
    all_hists['LHA_post_SD'], all_hists['track_LHA_post_SD'] = [], []
    all_hists['width_pre_SD'], all_hists['track_width_pre_SD'] = [], []
    all_hists['width_post_SD'], all_hists['track_width_post_SD'] = [], []
    all_hists['thrust_pre_SD'], all_hists['track_thrust_pre_SD'] = [], []
    all_hists['thrust_post_SD'], all_hists['track_thrust_post_SD'] = [], []

    all_hists['zg_10'], all_hists['track_zg_10'] = [], []
    all_hists['rg_10'], all_hists['track_rg_10'] = [], []
    all_hists['e1_10'], all_hists['track_e1_10'] = [], []
    all_hists['e2_10'], all_hists['track_e2_10'] = [], []
    all_hists['e05_10'], all_hists['track_e05_10'] = [], []

    all_hists['zg_05'], all_hists['track_zg_05'] = [], []
    all_hists['rg_05'], all_hists['track_rg_05'] = [], []
    all_hists['e1_05'], all_hists['track_e1_05'] = [], []
    all_hists['e2_05'], all_hists['track_e2_05'] = [], []
    all_hists['e05_05'], all_hists['track_e05_05'] = [], []

    all_hists['zg_20'], all_hists['track_zg_20'] = [], []
    all_hists['rg_20'], all_hists['track_rg_20'] = [], []
    all_hists['e1_20'], all_hists['track_e1_20'] = [], []
    all_hists['e2_20'], all_hists['track_e2_20'] = [], []
    all_hists['e05_20'], all_hists['track_e05_20'] = [], []

    all_hists['hardest_pT'] = []
    all_hists['uncor_hardest_pT'] = []
    all_hists['track_pT_D_pre_SD'] = []
    all_hists['track_pT_D_post_SD'] = []
    all_hists['pT_D_pre_SD'] = []
    all_hists['pT_D_post_SD'] = []

    for j in xrange(len(eta_boundaries)):

        default_85_conditions = [('hardest_eta', eta_boundaries[j]),
                                 ('hardest_pT', (85, None))]
        default_150_conditions = [('hardest_eta', eta_boundaries[j]),
                                  ('hardest_pT', (150, None))]
        default_250_conditions = [('hardest_eta', eta_boundaries[j]),
                                  ('hardest_pT', (250, None))]

        default_270_conditions = [('hardest_eta', eta_boundaries[j]),
                                  ('hardest_pT', (270, None))]

        eta_boundary_label = ""

        if (abs(eta_boundaries[j][0] - 0.0) < 1e-5):
             eta_boundary_label = r"\left| \eta \right| < " + \
                 str(eta_boundaries[j][1])
        else:
             eta_boundary_label = r"\left| \eta \right| < " + \
                 str(eta_boundaries[j][1])


        additional_text = [((0.9, 1.05), 'upper right',
                           "\n AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} > 270~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$")]
        all_hists['rg_10'].append(MODHist(copy.deepcopy(rg_10_hist), conditions=default_270_conditions, x_scale='log', use_prescale=True, x_label="$\\theta_g$",
                                          y_label="$ \\frac{\\theta_g}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} \\theta_g}$", x_range=(0.01, 1.0), y_range=(0., 0.6), legend_location=('upper left', (0., 1.0)), additional_text=additional_text))
        
    return all_hists


def multi_page_plot_hist_templates():

    all_hists = {}

    hardest_pT_hist = Hist(230, 0, 2300, title="pT")
    hardest_phi_hist = Hist(50, 0.0, 2 * np.pi, title="phi")
    hardest_eta_hist = Hist(50, -5, 5, title="eta")
    constituent_mul_hist, track_constituent_mul_hist = Hist(
        101, -0.5, 100.5, title="mul"), Hist(101, -0.5, 100.5, title="mul")
    pT_D_hist = Hist(50, 0, 1, title="pT_D")
    mass_hist, track_mass_hist = Hist(
        100, 0, 100, title="mass"), Hist(100, 0, 100, title="mass")

    lha_hist = Hist(50, 0, 1)
    width_hist = Hist(60, 0, 0.6)
    thrust_hist = Hist(50, 0, 0.5)

    jec_hist = Hist(100, 0.95, 1.20)
    hardest_area_hist = Hist(20, 0.7, 0.9)

    zg_hist = Hist(50, 0.0, 0.5)
    rg_hist = Hist(60, 0.0, 1.2)
    e1_hist = Hist(60, 0.0, 0.6)
    e05_hist = Hist(60, 0.0, 0.6)
    e2_hist = Hist(60, 0.0, 0.6)

    pT_boundaries = [85, 115, 150, 200, 250]
    # eta_boundaries = [(0.0, 1.9), (0.0, 2.4), (1.9, 2.4)]
    eta_boundaries = [(-2.4, 2.4)]

    all_hists['hardest_pT'] = []
    all_hists['uncor_hardest_pT'] = []
    all_hists['frac_pT_loss'] = []
    all_hists['jec'] = []
    all_hists['hardest_area'] = []
    all_hists['softkill_pT_loss'] = []
    all_hists['pT_after_SD'] = []
    all_hists['hardest_eta'] = []
    all_hists['hardest_phi'] = []
    all_hists['mul_pre_SD'], all_hists['track_mul_pre_SD'] = [], []
    all_hists['mul_post_SD'], all_hists['track_mul_post_SD'] = [], []
    all_hists['pT_D_pre_SD'], all_hists['track_pT_D_pre_SD'] = [], []
    all_hists['pT_D_post_SD'], all_hists['track_pT_D_post_SD'] = [], []
    all_hists['mass_pre_SD'], all_hists['track_mass_pre_SD'] = [], []
    all_hists['mass_post_SD'], all_hists['track_mass_post_SD'] = [], []

    all_hists['LHA_pre_SD'], all_hists['track_LHA_pre_SD'] = [], []
    all_hists['LHA_post_SD'], all_hists['track_LHA_post_SD'] = [], []
    all_hists['width_pre_SD'], all_hists['track_width_pre_SD'] = [], []
    all_hists['width_post_SD'], all_hists['track_width_post_SD'] = [], []
    all_hists['thrust_pre_SD'], all_hists['track_thrust_pre_SD'] = [], []
    all_hists['thrust_post_SD'], all_hists['track_thrust_post_SD'] = [], []

    all_hists['zg_10'], all_hists['track_zg_10'] = [], []
    all_hists['rg_10'], all_hists['track_rg_10'] = [], []
    all_hists['e1_10'], all_hists['track_e1_10'] = [], []
    all_hists['e2_10'], all_hists['track_e2_10'] = [], []
    all_hists['e05_10'], all_hists['track_e05_10'] = [], []

    all_hists['zg_05'], all_hists['track_zg_05'] = [], []
    all_hists['rg_05'], all_hists['track_rg_05'] = [], []
    all_hists['e1_05'], all_hists['track_e1_05'] = [], []
    all_hists['e2_05'], all_hists['track_e2_05'] = [], []
    all_hists['e05_05'], all_hists['track_e05_05'] = [], []

    all_hists['zg_20'], all_hists['track_zg_20'] = [], []
    all_hists['rg_20'], all_hists['track_rg_20'] = [], []
    all_hists['e1_20'], all_hists['track_e1_20'] = [], []
    all_hists['e2_20'], all_hists['track_e2_20'] = [], []
    all_hists['e05_20'], all_hists['track_e05_20'] = [], []

    for j in xrange(len(eta_boundaries)):

        default_90_conditions = [('hardest_eta', eta_boundaries[j]),
                                 ('hardest_pT', (90, None))]
        default_270_conditions = [('hardest_eta', eta_boundaries[j]),
                                 ('hardest_pT', (270, None))]
        default_150_conditions = [('hardest_eta', eta_boundaries[j]),
                                  ('hardest_pT', (150, None))]
        default_250_conditions = [('hardest_eta', eta_boundaries[j]),
                                  ('hardest_pT', (250, None))]

        eta_boundary_label = ""

        if (abs(eta_boundaries[j][0] - 0.0) < 1e-5):
             eta_boundary_label = r"\left| \eta \right| < " + \
                 str(eta_boundaries[j][1])
        else:
             eta_boundary_label = r"\left| \eta \right| < " + \
                 str(eta_boundaries[j][1])

        additional_text = [((0.70, 0.70), 'upper right',
                            "AK5; $" + eta_boundary_label  + "$ \n $p_T^{\mathrm{jet}} > 90~\mathrm{GeV}$")]
        all_hists['hardest_pT'].append(MODHist(copy.deepcopy(hardest_pT_hist), conditions=[('hardest_eta', eta_boundaries[j]),
                                 ('hardest_pT', (90, None))], use_prescale=True, x_label="$p_T$ [GeV]",
                                               y_label="Diff. Cross Section [pb/GeV]", y_scale='log', x_range=(0, 1800), y_range=(1e-8, 1e6), mark_regions=[(270, 1e0, 'right', 500000.0, 40)], additional_text=additional_text))
    
        additional_text = [((0.7, 0.95), 'upper right',
                            "AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} > 270~\mathrm{GeV}$")]
        all_hists['hardest_phi'].append(MODHist(copy.deepcopy(hardest_phi_hist), conditions=default_270_conditions, use_prescale=True, x_label="$\phi$", x_range=(
            0, 2 * np.pi), y_label="Diff. Cross Section [nb] ", y_range=(0, 1.0),  mark_regions=[(150, 1e4, 'right', 0.04, 25)], additional_text=additional_text, axes_label_pi=True, legend_location=('upper left', (0., 1.0))))


        additional_text = [((0.70, 0.95), 'upper right',
                            "AK5 \n $p_T^{\mathrm{jet}} > 270~\mathrm{GeV}$")]
        all_hists['hardest_eta'].append(MODHist(copy.deepcopy(hardest_eta_hist), conditions=[('hardest_pT', (270, None))], use_prescale=True, x_label="$\eta$", x_range=(
            -4., 4.), y_label="Diff. Cross Section [nb] ", y_scale='linear',  y_range=(0, 1.8), mark_regions=[(-2.4, 1.0, 'right', 0.8, 0.3), (2.4, 1.0, 'left', 0.8, -0.3)], additional_text=additional_text, legend_location=('upper left', (0., 1.0))))


        additional_text = [((0.70, 0.7), 'upper right',
                             "AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} > 270~\mathrm{GeV}$")]
        all_hists['mass_pre_SD'].append(MODHist(copy.deepcopy(mass_hist), conditions=default_270_conditions, use_prescale=True, x_label="Mass [GeV]",
                                                y_label="Diff. Cross Section [pb/GeV]", x_range=(0.0, 90.0), y_scale='linear',  y_range=(0.0, 160.0), additional_text=additional_text))


        additional_text = [((0.70, 0.7), 'upper right', "AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} > 270~\mathrm{GeV}$")]
        all_hists['mul_pre_SD'].append(MODHist(copy.deepcopy(constituent_mul_hist), conditions=default_270_conditions, use_prescale=True, x_label="Constituent Multiplicity", x_range=(
            0, 90.0), y_scale='linear', y_range=(0.0, 140.0), y_label="Diff. Cross Section [pb]", additional_text=additional_text, legend_location=('upper right', (1.0, 1.0))))

        additional_text = [((0.70, 0.70), 'upper right', "AK5; $" + eta_boundary_label + "$\n $p_T^{\mathrm{jet}} > 270~\mathrm{GeV}$")]
        all_hists['track_mass_pre_SD'].append(MODHist(copy.deepcopy(track_mass_hist), conditions=default_270_conditions, use_prescale=True, x_label="Track Mass [GeV]",
                                                y_label="Diff. Cross Section [pb/GeV]", x_range=(0.0, 60.0), y_scale='linear',  y_range=(0.0, 200.0), additional_text=additional_text))


        additional_text = [((0.7, 0.7), 'upper right',"AK5; $" + eta_boundary_label + "$\n $p_T^{\mathrm{jet}} > 270~\mathrm{GeV}$")]
        all_hists['track_mul_pre_SD'].append(MODHist(copy.deepcopy(track_constituent_mul_hist), conditions=default_270_conditions, use_prescale=True, x_label="Track Constituent Multiplicity", x_range=(
            0.0, 60.0), y_scale='linear', y_range=(0.0, 200.0), y_label="Diff. Cross Section [pb]", additional_text=additional_text, legend_location=('upper right', (1.0, 1.0))))


        additional_text = [((0.85, 0.75), 'upper right',
                            "\n AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} > 270~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$")]
        all_hists['zg_10'].append(MODHist(copy.deepcopy(zg_hist), conditions=default_270_conditions, use_prescale=True, x_label="$z_g$",
                                          y_label="$ \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), y_range=(0, 9), additional_text=additional_text))

    
        additional_text = [((0.85, 0.75), 'upper right',
                            "\n AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} > 270 ~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$")]
        all_hists['rg_10'].append(MODHist(copy.deepcopy(rg_hist), conditions=default_270_conditions, use_prescale=True, x_label="$r_g$",
                                          y_label="$ \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} r_g}$", x_range=(0.0, 1.0), y_range=(0, 7), additional_text=additional_text))


    return all_hists


def trigger_hists():
    all_hists = {}

    pT_hist = Hist(50, 0, 300)

    eta_boundary_label = "\left| \eta \\right| < 2.4"

    # all_hists['corr_hardest_pT'] = []

    jet_trigger_names = ["HLT_Jet370", "HLT_Jet300", "HLT_Jet240", "HLT_Jet190", "HLT_Jet150", "HLT_Jet110", "HLT_Jet80", "HLT_Jet60", "HLT_Jet30" ]

    qcd_names = ['QCD_Pt-15to30_TuneZ2_7TeV_pythia6',
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
                       'QCD_Pt-1800_TuneZ2_7TeV_pythia6']

    trigger_names = qcd_names + jet_trigger_names

    for trigger_name in qcd_names[:7]:
        additional_text = [
            (upper_left, 'upper left', "AK5; " + eta_boundary_label)]
        all_hists[trigger_name] = MODHist(Hist(500, 0, 3000), conditions=[],
                                          use_prescale=True, x_scale='log', additional_text=additional_text)


    for trigger_name in qcd_names[7:]:
        additional_text = [
            (upper_left, 'upper left', "AK5; " + eta_boundary_label)]
        all_hists[trigger_name] = MODHist(Hist(500, 0, 3000), conditions=[],
                                          use_prescale=True, x_scale='log', additional_text=additional_text)
        # all_hists[trigger_name] =  MODHist(Hist(50, 0, 300),
        # conditions=[(['trig_jet_matched', 1], lambda x, y: y == x)],
        # use_prescale=True, x_scale='log', x_label="Fractional $p_T$ Loss",
        # y_label="A.U.", y_range=(0., 1.2), additional_text=additional_text )
        

    for trigger_name in jet_trigger_names:

        additional_text = [
            (upper_left, 'upper left', "AK5; " + eta_boundary_label)]
        all_hists[trigger_name] = MODHist(Hist(200, 0, 2000),  conditions=[(['jet_quality', 1], lambda x, y: y >= x), ([
                                                        'trig_jet_matched', 1], lambda x, y: y == x)],
                                          use_prescale=True, x_scale='log', x_label="Fractional $p_T$ Loss", y_label="A.U.", y_range=(0., 1.2), additional_text=additional_text)


        #  all_hists[trigger_name + "_prescale"] =
        # MODHist(Hist(np.logspace(math.log(float(1e0), math.e), math.log(1e4,
        # math.e), (100 + 1), base=np.e)), conditions=[(['trig_jet_matched',
        # 1], lambda x, y: y == x)], use_prescale=False, x_scale='log',
        # y_scale='log', x_label="Trigger Prescale", y_label="A.U.",
        # y_range=(0., 1.2), additional_text=additional_text )
        

    return all_hists

def get_pfc_hists():
    pfc_hists = {}

    pfc_hists['pfc_pT'] = []
    pfc_hists['pfc_eta'] = []

    pT_boundaries = [85, 115, 150, 200, 250]

    eta_boundary_label = "\left| \eta \\right| < 2.4"

    # Charged PFCs


    additional_text = [((0.78, 0.6), 'upper right',
                        "Charged PFCs \n AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} > 220~\mathrm{GeV}$")]
    pfc_hists['pfc_pT'].append(MODHist(Hist(50, 0, 10), conditions=[('pfc_pdgId', ("in", 11, -11, 13,
                                                                                  -13, 15, -15, 211, -211, 321, -321,
                                                                                  2212, -2212, 3222, -3222, 3312, -3312, 3334, -3334)),
                                                                   ('jet_eta', (-2.4, 2.4)), ('hardest_pT', (220, None))],
                                       mark_regions=[(1.0, 0.45, 'right', 0.1, 0.5)], use_prescale=True,
                                       x_label="Charged PFC $p_T$ [GeV]", y_label="A.U.",
                                       x_range=(0, 10),  y_scale='linear', y_range=(0, 120000), additional_text=additional_text))


    additional_text = [((0.78, 0.6), 'upper right',
                        "Neutral PFCs \n AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} > 220~\mathrm{GeV}$")]
    pfc_hists['pfc_pT'].append(MODHist(Hist(50, 0, 10), conditions=[('pfc_pdgId', ("not in", 11, -11, 13,
                                                                                  -13, 15, -15, 211, -211, 321, -321,
                                                                                  2212, -2212, 3222, -3222, 3312, -3312, 3334, -3334)),
                                                                     ('jet_eta', (-2.4, 2.4)),
                                                                     ('hardest_pT', (220, None))],
                                       use_prescale=True, y_scale='linear', x_label="Neutral PFC $p_T$ [GeV]",
                                       y_label="A.U.", x_range=(0, 10), y_range=(0, 120000), additional_text=additional_text))

    return pfc_hists


def two_dim_hists():
    all_hists = {}

    zg_hist = Hist2D(25, 0.0, 0.5, 25, 0.0, 1.0)

    pT_boundaries = [85, 115, 150, 200, 250]

    eta_boundary_label = "\left| \eta \\right| < 2.4"

    all_hists[('zg_10', 'rg_10')] = []
    all_hists[('track_zg_10', 'track_rg_10')] = []

    default_85_conditions = [('hardest_eta', (-2.4, 2.4)),
                             ('hardest_pT', (85, None))]
    default_150_conditions = [('hardest_eta', (-2.4, 2.4)),
                              ('hardest_pT', (150, None))]
    default_250_conditions = [('hardest_eta', (-2.4, 2.4)),
                              ('hardest_pT', (250, None))]

    for i in xrange(len(pT_boundaries) - 1):

        default_conditions = [('hardest_eta', (-2.4, 2.4)),
                              ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))]

        additional_text = [((1.0, 0.62), 'upper right', "$p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$ \n AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} \in [" + str(
            pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$")]
        all_hists[('zg_10', 'rg_10')].append(MODHist(copy.deepcopy(zg_hist), conditions=default_conditions, use_prescale=False, x_label="$z_g$",
                                                     y_label="$ \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text))

        additional_text = [((1.0, 0.62), 'upper right', "$p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$ \n AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} \in [" + str(
            pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$")]
        all_hists[('track_zg_10', 'track_rg_10')].append(MODHist(copy.deepcopy(zg_hist), conditions=default_conditions, use_prescale=False, x_label="$z_g$",
                                                                 y_label="$ \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text))

    additional_text = [((1.0, 0.62), 'upper right',
                        "$p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$ \n AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} > 85~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$")]
    all_hists[('zg_10', 'rg_10')].append(MODHist(copy.deepcopy(zg_hist), conditions=default_85_conditions, use_prescale=True, x_label="$z_g$",
                                                 y_label="$ \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text))

    additional_text = [((1.0, 0.62), 'upper right',
                        "$p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$ \n AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} > 150~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$")]
    all_hists[('zg_10', 'rg_10')].append(MODHist(copy.deepcopy(zg_hist), conditions=default_150_conditions, use_prescale=True, x_label="$z_g$",
                                                 y_label="$ \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text))

    additional_text = [((1.0, 0.62), 'upper right',
                        "$p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$ \n AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} > 250~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$")]
    all_hists[('zg_10', 'rg_10')].append(MODHist(copy.deepcopy(zg_hist), conditions=default_250_conditions, use_prescale=True, x_label="$z_g$",
                                                 y_label="$ \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text))

    additional_text = [((1.0, 0.62), 'upper right',
                        "$p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$ \n AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} > 85~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$")]
    all_hists[('track_zg_10', 'track_rg_10')].append(MODHist(copy.deepcopy(zg_hist), conditions=default_85_conditions, use_prescale=True,
                                                             x_label="$z_g$", y_label="$ \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text))

    additional_text = [((1.0, 0.62), 'upper right',
                        "$p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$ \n AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} > 150~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$")]
    all_hists[('track_zg_10', 'track_rg_10')].append(MODHist(copy.deepcopy(zg_hist), conditions=default_150_conditions, use_prescale=True,
                                                             x_label="$z_g$", y_label="$ \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text))

    additional_text = [((1.0, 0.62), 'upper right',
                        "$p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$ \n AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} > 250~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$")]
    all_hists[('track_zg_10', 'track_rg_10')].append(MODHist(copy.deepcopy(zg_hist), conditions=default_250_conditions, use_prescale=True,
                                                             x_label="$z_g$", y_label="$ \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text))

    return all_hists


def two_dim_log_hists():
    all_hists = {}

    # zg_hist = Hist2D(50, 0.0, 0.5, 50, 0, 1.0)
    # zg_10_hist = Hist(np.logspace(math.log(float(0.01), math.e),
    # math.log(0.5, math.e), 50, base=np.e))
    zg_hist = Hist2D(np.logspace(math.log(float(0.1), math.e), math.log(0.5, math.e), (25 + 1), base=np.e),
                     np.logspace(math.log(float(0.01), math.e), math.log(1.0, math.e), (25 + 1), base=np.e))

    pT_boundaries = [70, 110, 130, 190, 220, 275, 350, 420]

    default_70_conditions = [('hardest_eta', (-2.4, 2.4)),
                             ('hardest_pT', (70, None))]
    default_150_conditions = [('hardest_eta', (-2.4, 2.4)),
                              ('hardest_pT', (150, None))]
    default_250_conditions = [('hardest_eta', (-2.4, 2.4)),
                              ('hardest_pT', (250, None))]

    all_hists[('zg_10', 'rg_10')] = []
    all_hists[('track_zg_10', 'track_rg_10')] = []

    eta_boundary_label = "\left| \eta \\right| < 2.4"

    for i in xrange(len(pT_boundaries) - 1):

        default_conditions = [('hardest_eta', (-2.4, 2.4)),
                              ('hardest_pT', (pT_boundaries[i], pT_boundaries[i + 1]))]

        additional_text = [((1.0, 0.62), 'upper right', "$p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$ \n AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} \in [" + str(
            pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$")]
        all_hists[('zg_10', 'rg_10')].append(MODHist(copy.deepcopy(zg_hist), conditions=default_conditions, use_prescale=False, x_label="$z_g$",
                                                     y_label="$ \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text))

        additional_text = [((1.0, 0.62), 'upper right', "$p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$ \n AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} \in [" + str(
            pT_boundaries[i]) + ", " + str(pT_boundaries[i + 1]) + "]~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$")]
        all_hists[('track_zg_10', 'track_rg_10')].append(MODHist(copy.deepcopy(zg_hist), conditions=default_conditions, use_prescale=False, x_label="$z_g$",
                                                                 y_label="$ \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text))

    additional_text = [((1.0, 0.62), 'upper right',
                        "$p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$ \n AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} > 70~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$")]
    all_hists[('zg_10', 'rg_10')].append(MODHist(copy.deepcopy(zg_hist), conditions=default_70_conditions, use_prescale=True, x_label="$z_g$",
                                                 y_label="$ \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text))

    additional_text = [((1.0, 0.62), 'upper right',
                        "$p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$ \n AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} > 150~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$")]
    all_hists[('zg_10', 'rg_10')].append(MODHist(copy.deepcopy(zg_hist), conditions=default_150_conditions, use_prescale=True, x_label="$z_g$",
                                                 y_label="$ \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text))

    additional_text = [((1.0, 0.62), 'upper right',
                        "$p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$ \n AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} > 250~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$")]
    all_hists[('zg_10', 'rg_10')].append(MODHist(copy.deepcopy(zg_hist), conditions=default_250_conditions, use_prescale=True, x_label="$z_g$",
                                                 y_label="$ \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text))

    additional_text = [((1.0, 0.62), 'upper right',
                        "$p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$ \n AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} > 70~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$")]
    all_hists[('track_zg_10', 'track_rg_10')].append(MODHist(copy.deepcopy(zg_hist), conditions=default_70_conditions, use_prescale=True,
                                                             x_label="$z_g$", y_label="$ \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text))

    additional_text = [((1.0, 0.62), 'upper right',
                        "$p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$ \n AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} > 150~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$")]
    all_hists[('track_zg_10', 'track_rg_10')].append(MODHist(copy.deepcopy(zg_hist), conditions=default_150_conditions, use_prescale=True,
                                                             x_label="$z_g$", y_label="$ \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text))

    additional_text = [((1.0, 0.62), 'upper right',
                        "$p_T^{\mathrm{PFC}} > 1.0~\mathrm{GeV}$ \n AK5; $" + eta_boundary_label + "$ \n $p_T^{\mathrm{jet}} > 250~\mathrm{GeV}$ \n mMDT / SD_{\\beta = 0}: z_{\mathrm{cut}} = 0.1$")]
    all_hists[('track_zg_10', 'track_rg_10')].append(MODHist(copy.deepcopy(zg_hist), conditions=default_250_conditions, use_prescale=True,
                                                             x_label="$z_g$", y_label="$ \\frac{1}{\sigma} \\frac{ \mathrm{d} \sigma}{ \mathrm{d} z_g}$", x_range=(0.0, 0.6), additional_text=additional_text))

    return all_hists


def get_hist_template(var):
    all_hists = all_hist_templates()
    return all_hists[var]


if __name__ == "__main__":

    hist = Hist2D(1, 0.1, 0.2, 1, 0, 0.2)

    hist.fill_array([[0.1, 0.3]])
