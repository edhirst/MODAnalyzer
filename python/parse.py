from __future__ import division

import math
import time
import sys
import hists

import subprocess
import os.path

from MODPlot import *


from rootpy.io import File as TFile


import matplotlib.pyplot as plt
import rootpy.plotting.root2matplotlib as rplt



def parse_file(input_files, output_filename, data_type, all_hists, log_hists):

    # We read all th efiles in input_file line by line, and for each line, we fill the corresponding histograms.
    # This is desirable to creating lists of values since this will not hold
    # anything in memory.

    keywords = []
    keywords_set = False

    parse_count = 1

    if data_type == '2011':
        trigger_luminosity_dictionary = trigger_luminosity_dictionary_2011()
    #elif (data_type == 'sim_gen' or data_type == 'sim_pfc'):
    #    mod_file_to_crosssection = mod_file_to_crosssection_sim(data_type)

    for input_file in input_files:

        with open(input_file) as infile:

            print "Parsing {}".format(input_file)

            line_number = 0

            largest_pT = 0.0

            for line in infile:

                start = time.time()

                if len(line.strip()) == 0:
                    continue
                line_number += 1

                if line_number % 10000 == 0:
                    print "at line number: " + str(line_number)

                try:
                    # if True:
                    numbers = line.split()

                    if numbers[0] == "#" and (not keywords_set):
                        keywords = numbers[1:]
                        keywords_set = True

                        # Build a dictionary of keyword indices.
                        keywords_index_dictionary = {}
                        for keyword_index, keyword in enumerate(keywords):
                            keywords_index_dictionary[keyword] = keyword_index

                    elif numbers[0] == "Entry":

                        pT_of_this_event = float(numbers[keywords_index_dictionary['hardest_pT']])

                        prescale_to_use = 0.0
                        
                        if data_type == '2011':
                            trigger = numbers[keywords_index_dictionary['trigger_fired']]
                            prescale_to_use = 1000000.0/trigger_luminosity_dictionary[trigger]
                        elif (data_type == 'sim_pfc' or data_type == 'sim_gen'):
                            input_file = numbers[keywords_index_dictionary['filename']]
                            prescale_to_use = mod_file_to_crosssection[input_file]
                                                    
                        for some_hists in [all_hists, log_hists]:

                            for keyword in some_hists.keys():

                                for mod_hist in some_hists[keyword]:

                                    hist = mod_hist.hist()

                                    # Now, see whether or not the conditions are
                                    # met.
                                    condition_satisfied = True
                                    conditions = mod_hist.conditions()

                                    for condition_keyword, condition_boundaries in conditions:
                                        keyword_index = keywords_index_dictionary[
                                            condition_keyword]

                                        if condition_boundaries[0] != None and float(numbers[keyword_index]) < float(condition_boundaries[0]):
                                            condition_satisfied = False

                                        if condition_boundaries[1] != None and float(numbers[keyword_index]) > float(condition_boundaries[1]):
                                            condition_satisfied = False

                                    if condition_satisfied:

                                        if pT_of_this_event > largest_pT:
                                            largest_pT = pT_of_this_event

                                        parse_count += 1
                                        x = float(numbers[keywords_index_dictionary[keyword]])

                                        hist.fill_array([x], [prescale_to_use])


                except Exception as e:
                    pass
                    print "Some exception occured!",
                    print e.message

                end = time.time()

                # print "Took {} seconds for current line.".format(end - start)

    print "Total parsed count is", parse_count

    return all_hists, log_hists


def write_to_root_files(parsed_hists, parsed_log_hists, output_filename):
    start = time.time()

    f = TFile(output_filename[0], "RECREATE")

    for var in parsed_hists.keys():

        index = 0

        for mod_hist in parsed_hists[var]:
            hist = copy.deepcopy(mod_hist.hist())
            hist.SetName("{}#{}".format(var, index))
            hist.Write()

            index += 1

    f.Close()

    f = TFile(output_filename[1], "RECREATE")

    for var in parsed_log_hists.keys():

        index = 0

        for mod_hist in parsed_log_hists[var]:
            hist = copy.deepcopy(mod_hist.hist())
            hist.SetName("{}#{}".format(var, index))
            hist.Write()

            index += 1

    f.Close()

    end = time.time()

    # print "\n" * 3
    # print "Took {} seconds to write down to files.".format(end - start)


def parse_to_root_file(input_files, output_filename, data_type, hist_templates):

    # print "Parsing {} to {}".format(input_filename, output_filename)

    # First, get the already histogrammed hists.

    hists, log_hists = copy.deepcopy(hist_templates[0]), copy.deepcopy(hist_templates[1])

    parsed_hists, parsed_log_hists = parse_file(input_files, output_filename, data_type, hists, log_hists)

    write_to_root_files(parsed_hists, parsed_log_hists, output_filename)


def root_file_to_hist(input_filename, hist_templates):

    hists = copy.deepcopy(hist_templates)

    root_file = TFile(input_filename, "read")

    for var in hists.keys():

        index = 0

        for mod_hist in hists[var]:
            hist_name = "{}#{}".format(var, index)

            # Get hist from ROOT file.
            hist = root_file.Get(hist_name)

            mod_hist.replace_hist(hist)

            index += 1

    return hists


def parse_to_root_files():

    input_directory = sys.argv[1]
    data_files = os.listdir(input_directory)

    output_directory = sys.argv[2]

    data_files_2011 = []
    data_files_sim_pfc = []
    data_files_sim_gen = []

    for data_file in data_files:
        if data_file.endswith(".dat"):
            if 'pfc' in data_file:
                data_files_sim_pfc.append(input_directory + data_file)
            elif 'gen' in data_file:
                data_files_sim_gen.append(input_directory + data_file)
            else:
                data_files_2011.append(input_directory + data_file)
            
    hist_templates = hists.multi_page_plot_hist_templates()

    log_hist_templates = hists.multi_page_log_plot_hist_templates()
 
    parse_to_root_file(input_files=data_files_2011, output_filename=(output_directory + "data_2011.root",
                                                                  output_directory + "data_2011_log.root"), data_type = "2011", hist_templates=(hist_templates,
                                                                                                                       log_hist_templates))
    parse_to_root_file(input_files=data_files_sim_pfc, output_filename=(output_directory + "data_sim_pfc.root",
                                                                  output_directory + "data_sim_pfc_log.root"), data_type = "sim_pfc", hist_templates=(hist_templates,
                                                                                                                       log_hist_templates))
    parse_to_root_file(input_files=data_files_sim_gen, output_filename=(output_directory + "data_sim_gen.root",
                                                                  output_directory + "data_sim_gen_log.root"), data_type = "sim_gen", hist_templates=(hist_templates,
                                                                                                                       log_hist_templates))

def load_root_files_to_hist(log=False):

    output_directory = sys.argv[1]

    if not log:
        hist_templates = hists.multi_page_plot_hist_templates()
	filenames = ["data_2011.root", "data_sim_pfc.root", "data_sim_gen.root"]
	filenames = ["data_sim_all_new.root"]
    else:
        hist_templates = hists.multi_page_log_plot_hist_templates()
        filenames = ["data_2011_log.root", "data_sim_pfc_log.root", "data_sim_gen_log.root"]

    return [root_file_to_hist(output_directory + filename, hist_templates) for filename in filenames]


def file_merge(output, input_files):
    with open(output, 'w') as merged_file:
        for fname in input_files:
            with open(fname) as infile:
                for line in infile:
                    merged_file.write(line)


def trigger_luminosity_dictionary_2011():
    mod_file_with_trigger = [line.rstrip('\n') for line in open("effective_luminosity_by_trigger.csv")]
    mod_trigger_luminosities = {}
    for mod_trigger_lumi in mod_file_with_trigger:
        mod_trigger_luminosities[(mod_trigger_lumi.split(',')[0], mod_trigger_lumi.split(',')[1])] = mod_trigger_lumi.split(',')[2]
    trigger_luminosity_total = {}
    for mod_trigger in mod_trigger_luminosities:
        if mod_trigger[1] in trigger_luminosity_total:
            #current_effective_lumi = trigger_luminosity_total[mod_trigger[1]]
            trigger_luminosity_total[mod_trigger[1]] += mod_trigger_luminosities[mod_trigger]
        else:
            trigger_luminosity_total[mod_trigger[1]] = mod_trigger_luminosities[mod_trigger]
    return trigger_luminosity_total

def mod_file_to_crosssection_sim(flag):
    output_file_event_count = [line.rstrip('\n') for line in open("event_count_by_pythia_and_mod.csv")]
    pythia_cross_sections = {'QCD_Pt-30to50_TuneZ2_7TeV_pythia6': 53122370.0,
                             'QCD_Pt-50to80_TuneZ2_7TeV_pythia6': 6359119.0,
                             'QCD_Pt-80to120_TuneZ2_7TeV_pythia6': 784265.0,
                             'QCD_Pt-120to170_TuneZ2_7TeV_pythia6': 115134.0,
                             'QCD_Pt-170to300_TuneZ2_7TeV_pythia6': 24262.8,
                             'QCD_Pt-300to470_TuneZ2_7TeV_pythia6': 1168.49,
                             'QCD_Pt-470to600_TuneZ2_7TeV_pythia6': 70.2242,
                             'QCD_Pt-600to800_TuneZ2_7TeV_pythia6': 15.5537,
                             'QCD_Pt-800to1000_TuneZ2_7TeV_pythia6': 1.84369,
                             'QCD_Pt-1000to1400_TuneZ2_7TeV_pythia6': 0.332105,
                             'QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6': 0.0108721}
    # Build dictionary of mod file, pythia_set and associated number of events
    # We do it this in way as protection in case the csv file has been appended
    # with same data set twice accidentally. keep it fresh. 
    mod_file_and_pythia_set_count = {}
    for line in output_file_event_count:
        output_file = line.split(',')[0]
        event_count = line.split(',')[1]
        mod_name_start_index = output_file.find("pythia6") + len(pythia6)
        mod_name_end_index = output_file.find("_"+flag)
        mod_name = outputfile[mod_name_start_index: mod_name_end_index]+".mod"
        pythia_set = outputfile[:mod_name_start_index]
        mod_file_and_pythia_set_count[(mod_name, pythia_set)] = event_count

    for mod_file, pythia_set in mod_file_and_pythia_set_count:
        if pythia_set in pythia_set_to_event_count:
            pythia_set_to_event_count[pythia_set] += event_count
        else:
            pythia_set_to_event_count[pythia_set] = event_count
    mod_file_to_cross_section = {}
    for mod_file, pythia_set in mod_file_and_pythia_set_count:
        mod_file_to_cross_section[mod_file] = pythia_cross_sections[pythia_set]/pythia_set_to_event_count[pythia_set]        
        
    return mod_file_to_cross_section

if __name__ == "__main__":

    parse_to_root_files()

    pass
