from __future__ import division

import math
import time
import sys
import hists


from MODPlot import *


import parse
import pfc_parse


def compile_sources(parsed_hists):

    data_hists = parsed_hists[0]

    return [data_hists]


def compile_data_and_pythia(all_hists, variables):

    data_hists, pythia_hists = all_hists[0], all_hists[1]

    data, pythia = [], []
    for var in variables:
        data.append(data_hists[var])
        pythia.append(pythia_hists[var])

    compilation = []
    for i in range(len(data[0])):
        temp = [data[0][i], data[1][i], pythia[0][i], pythia[1][i]]
        compilation.append(temp)

    return compilation


def compile_hists(var, parsed_hists, x_scale='linear'):

    compilation = []

    data_hists, pythia_hists, herwig_hists, sherpa_hists = parsed_hists[
        0], parsed_hists[1], parsed_hists[2], parsed_hists[3]

    # print data_hists, data_hists, var

    max_index = len(data_hists[var])

    for i in range(max_index):

        sub_list = [data_hists[var][i], pythia_hists[var][
            i], herwig_hists[var][i], sherpa_hists[var][i]]
        # sub_list = [ data_hists[var][i] ]

        compilation.append(sub_list)

    return compilation


def compile_hists_with_theory(var, parsed_hists, x_scale='linear'):

    compilation = []

    data_hists, theory_hists, pythia_hists, herwig_hists, sherpa_hists = parsed_hists[
        0], parsed_hists[1], parsed_hists[2], parsed_hists[3], parsed_hists[4]
    max_index = len(data_hists[var])

    for i in range(max_index):

        theory_var = var

        if "track" in theory_var:
            theory_var = theory_var.replace("track_", "")

        # Get the correct variable name to use for theory.
        theory_var = theory_var.split("_")[0]

        if data_hists[var][i].conditions()[1][1][1] == None:
            # Open pT boundaries. Use long.
            # Don't hardcode position of pT condition. Find the pT condition
            # using "" in.
            theory_var += "_" + \
                str(data_hists[var][i].conditions()[1][1][0]) + "long"
        else:
            # Closed pT boundaries.
            # Don't hardcode position of pT condition. Find the pT condition
            # using "" in.
            theory_var += "_" + str(data_hists[var][i].conditions()[1][1][0])

        if x_scale == "log":

            theory_var = theory_var.split(
                "_")[0] + "_log" + theory_var.split("_")[1]

            # theory_var += "log"

        sub_list = [data_hists[var][i], theory_hists[theory_var], pythia_hists[
            var][i], herwig_hists[var][i], sherpa_hists[var][i]]

        compilation.append(sub_list)

    return compilation


default_dir = "plots/Version 6/"
output_directory = sys.argv[1]
data_file = output_directory + "analyzed_data.dat"
print(output_directory)


start = time.time()


parsed_file = parse.parse_to_root_files(output_directory, data_file)
parsed_linear = parse.load_root_files_to_hist(log=False)
parsed_hists = compile_sources(parsed_linear)


pfc_parsed = pfc_parse.parse_pfc_to_root_files(output_directory, data_file)
parsed_pfc = pfc_parse.load_pfc_root_files_to_hist()
parsed_pfc_hists = compile_sources(parsed_pfc)


end = time.time()


print "Finished parsing all files in {} seconds. Now plotting them!".format(end - start)


start = time.time()

create_multi_page_plot(filename=default_dir + "hardest_jet_phi_all_linear.pdf", hists=compile_hists('hardest_phi', parsed_hists))

create_multi_page_plot(filename=default_dir + "pfc_eta.pdf",
                        hists=compile_hists('pfc_eta', parsed_pfc_hists))


end = time.time()
print "Finished all plotting in {} seconds.".format(end - start)
