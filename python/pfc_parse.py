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



def parse_file(input_files, output_filename, data_type, all_hists):

    # We read all th efiles in input_file line by line, and for each line, we fill the corresponding histograms.
    # This is desirable to creating lists of values since this will not hold
    # anything in memory.

    keywords = []
    keywords_set = False

    parse_count = 1

    #print(data_type)

    if data_type == '2011':
        trigger_luminosity_dictionary = trigger_luminosity_dictionary_2011(input_files)
        print(trigger_luminosity_dictionary)
    elif (data_type == 'sim_gen' or data_type == 'sim_pfc'):
        mod_file_to_crosssection = mod_file_to_crosssection_sim(data_type, input_files)

        #print(mod_file_to_crosssection)

    for input_file in input_files:

        with open(input_file) as infile:

            #print "Parsing {}".format(input_file)

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

                        

                        prescale_to_use = 0.0
                        
                        if data_type == '2011':
                            trigger = numbers[keywords_index_dictionary['trigger_fired']]
                            prescale_to_use = 1000000.0/trigger_luminosity_dictionary[trigger]
                        elif (data_type == 'sim_pfc' or data_type == 'sim_gen'):
                            short_name = input_file.split('/')[-1]
                            #input_file = numbers[keywords_index_dictionary['filename']]
                            prescale_to_use = mod_file_to_crosssection[short_name]
                                                    
                        for some_hists in [all_hists]:

                            for keyword in some_hists.keys():
                            #for keyword in ['track_mass_pre_SD']:

                                for mod_hist in some_hists[keyword]:

                                    #print(some_hists.keys())

                                    hist = mod_hist.hist()
                                    conditions = mod_hist.conditions()

                                    # Now, see whether or not the conditions are
                                    # met.
                                    condition_satisfied = 1
                                    for condition_keyword, condition_boundaries in conditions:
                                        keyword_index = keywords.index(condition_keyword)


                                        #print(condition_keyword, condition_boundaries, numbers[keyword_index])

                                        if len(condition_boundaries) <= 2:
                                            if condition_boundaries[0] == None and condition_boundaries[1] != None:
                                                condition_satisfied *= int(
                                                    float(numbers[keyword_index]) < float(condition_boundaries[1]))
                                            elif condition_boundaries[0] != None and condition_boundaries[1] == None:
                                                condition_satisfied *= int(
                                                    float(numbers[keyword_index]) > float(condition_boundaries[0]))
                                            elif condition_boundaries[0] == None and condition_boundaries[1] == None:
                                                condition_satisfied *= 1
                                            elif condition_boundaries[0] != None and condition_boundaries[1] != None:
                                                condition_satisfied *= int(float(numbers[keyword_index]) > float(condition_boundaries[
                                                                           0]) and float(numbers[keyword_index]) < float(condition_boundaries[1]))
                                        else:

                                            if condition_boundaries[0] == "in":
                                                condition_satisfied *= int(
                                                    int(numbers[keyword_index]) in list(condition_boundaries[1:]))
                                            else:
                                                condition_satisfied *= int(
                                                    int(numbers[keyword_index]) not in list(condition_boundaries[1:]))

                                    condition_satisfied = bool(condition_satisfied)

                                    if condition_satisfied and not float(numbers[keywords.index('hardest_pT')]) >= 220:
                                        print("WE HAVE A PROBLEM")
                                        break

                                    if condition_satisfied:

                                        #print(condition_satisfied)

                                        parse_count += 1

                                        x = float(numbers[keywords_index_dictionary[keyword]])

                                        #print(x, prescale_to_use)

                                        if keyword == 'hardest_eta' or keyword == 'hardest_phi':

                                            hist.fill_array([x], [prescale_to_use/(1000.0)])
                                        else:
                                            hist.fill_array([x], [prescale_to_use])

                except Exception as e:
                    pass
                    print "Some exception occured!",
                    print e.message

                end = time.time()

                # print "Took {} seconds for current line.".format(end - start)

    print "Total parsed count is", parse_count

    return all_hists


def write_to_root_files(parsed_hists,  output_filename):
    start = time.time()

    f = TFile(output_filename, "RECREATE")

    for var in parsed_hists.keys():

        index = 0

        for mod_hist in parsed_hists[var]:
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

    hists = copy.deepcopy(hist_templates)

    parsed_hists = parse_file(input_files, output_filename, data_type, hists)

    write_to_root_files(parsed_hists, output_filename)


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
            
    hist_templates = hists.get_pfc_hists()

    parse_to_root_file(input_files=data_files_2011, output_filename=(output_directory + "data_pfc_2011.root"), data_type = "2011", hist_templates=(hist_templates))
    parse_to_root_file(input_files=data_files_sim_pfc, output_filename=(output_directory + "data_pfc_sim_pfc.root"), data_type = "sim_pfc", hist_templates=(hist_templates))
    parse_to_root_file(input_files=data_files_sim_gen, output_filename=(output_directory + "data_pfc_sim_gen.root"), data_type = "sim_gen", hist_templates=(hist_templates))

def load_root_files_to_hist():

    output_directory = sys.argv[1]

    hist_templates = hists.get_pfc_hists()
    filenames = ["data_pfc_2011.root", "data_pfc_sim_pfc.root", "data_pfc_sim_gen.root"]

    return [root_file_to_hist(output_directory + filename, hist_templates) for filename in filenames]


def file_merge(output, input_files):
    with open(output, 'w') as merged_file:
        for fname in input_files:
            with open(fname) as infile:
                for line in infile:
                    merged_file.write(line)


def trigger_luminosity_dictionary_2011(input_files):
    input_files_short = []
    for input_file in input_files:
        input_files_short.append(input_file.split('/')[-1])
    print(input_files_short)
    mod_file_with_trigger = [line.rstrip('\n') for line in open("/home/preksha/Documents/mengproject/MODAnalyzer/effective_luminosity_by_trigger.csv")]
    mod_trigger_luminosities = {}
    for mod_trigger_lumi in mod_file_with_trigger:
        if mod_trigger_lumi.replace(" ",""):
            #print(mod_trigger_lumi)
            mod_file = mod_trigger_lumi.split(',')[0].replace('.mod', '.dat')
            trigger = mod_trigger_lumi.split(',')[1]
            lumi = mod_trigger_lumi.split(',')[2]
            if mod_file in input_files_short:
                mod_trigger_luminosities[(mod_file, trigger)] = lumi
    trigger_luminosity_total = {}
    for mod_trigger in mod_trigger_luminosities:
        if mod_trigger[1] in trigger_luminosity_total:
            #current_effective_lumi = trigger_luminosity_total[mod_trigger[1]]
            trigger_luminosity_total[mod_trigger[1]] += float(mod_trigger_luminosities[mod_trigger])
        else:
            trigger_luminosity_total[mod_trigger[1]] = float(mod_trigger_luminosities[mod_trigger])
    return trigger_luminosity_total

def mod_file_to_crosssection_sim(flag, input_files):
    input_files_short = []
    for input_file in input_files:
        input_files_short.append(input_file.split('/')[-1])
    output_file_event_count = [line.rstrip('\n') for line in open("event_count_by_pythia_and_mod.csv")]
    pythia_cross_sections = {'QCD_Pt-15to30_TuneZ2_7TeV_pythia6': 815912800.0, 
                             'QCD_Pt-30to50_TuneZ2_7TeV_pythia6': 53122370.0,
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

    #print(input_files_short)
    # Build dictionary of mod file, pythia_set and associated number of events
    # We do it this in way as protection in case the csv file has been appended
    # with same data set twice accidentally.
    mod_file_and_pythia_set_count = {}
    for line in output_file_event_count:
        output_file = line.split(',')[0]
        #print(output_file)
        if output_file in input_files_short:
            event_count = line.split(',')[1]
            mod_name_start_index = output_file.find("pythia6") + len("pythia6")
            #mod_name_end_index = output_file.find("_"+flag)
            #mod_name = output_file[mod_name_start_index: mod_name_end_index]+".mod"
            pythia_set = output_file[:mod_name_start_index]
            mod_file_and_pythia_set_count[(output_file, pythia_set)] = event_count

    #print(mod_file_and_pythia_set_count)
    pythia_set_to_event_count = {}
    for mod_file_pythia_set in mod_file_and_pythia_set_count:
        pythia_set = mod_file_pythia_set[1]
        event_count = float(mod_file_and_pythia_set_count[mod_file_pythia_set])
        if pythia_set in pythia_set_to_event_count:
            pythia_set_to_event_count[pythia_set] += event_count
        else:
            pythia_set_to_event_count[pythia_set] = event_count

    #print(pythia_set_to_event_count)
    mod_file_to_cross_section = {}
    for mod_file, pythia_set in mod_file_and_pythia_set_count:
        #print(pythia_cross_sections[pythia_set], pythia_set_to_event_count[pythia_set], pythia_set)
        mod_file_to_cross_section[mod_file] = pythia_cross_sections[pythia_set]/pythia_set_to_event_count[pythia_set]        
        
    return mod_file_to_cross_section

if __name__ == "__main__":

    parse_to_root_files()

    pass
