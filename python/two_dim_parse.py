from __future__ import division

import math
import time
import sys
import hists


from MODPlot import *


from rootpy.io import File as TFile



import matplotlib.pyplot as plt
import rootpy.plotting.root2matplotlib as rplt



def parse_file(input_files, all_hists, log_hists):
	
	# We read the file line by line, and for each line, we fill the corresponding histograms.
	# This is desirable to creating lists of values since this will not hold anything in memory. 

	keywords = []
	keywords_set = False

        trigger_luminosity_dictionary = trigger_luminosity_dictionary_2011(input_files)
        #mod_file_to_crosssection = mod_file_to_crosssection_sim("sim_pfc", input_files)

	for input_file in input_files:

                with open(input_file) as infile:

                        line_number = 0

                        for line in infile:

                                line_number += 1

                                if line_number % 100000 == 0:
                                        print "At line number {}".format(line_number)

                                try:
                                        numbers = line.split()

                                        if numbers[0] == "#" and (not keywords_set):
                                                keywords = numbers[2:]
                                                keywords_set = True

                                        elif numbers[0] == "Entry":

                                                trigger_index = keywords.index("trigger_fired") + 1
                                                pT_index = keywords.index("hardest_pT") + 1
                                                zg_10_index = keywords.index("zg_10") + 1
                                                track_zg_10_index = keywords.index("track_zg_10") + 1
                                                rg_10_index = keywords.index("rg_10") + 1
                                                track_rg_10_index = keywords.index("track_rg_10") + 1
                                                eta_index = keywords.index("hardest_eta") + 1

                                                pT_of_this_event = float(numbers[pT_index])
                                                zg_10 = float(numbers[zg_10_index])
                                                track_zg_10 = float(numbers[track_zg_10_index])
                                                rg_10 = float(numbers[rg_10_index])
                                                track_rg_10 = float(numbers[track_rg_10_index])
                                                #prescale = float(numbers[prescale_index])
                                                eta = float(numbers[eta_index])


                                                current_line_trig_name = numbers[trigger_index]
                                                #print(current_line_trig_name)

                                                prescale = 1000000.0/trigger_luminosity_dictionary[current_line_trig_name]
                                                #print(zg_10, rg_10, prescale)

                                                #short_name = input_file.split('/')[-1]
                                                #prescale = mod_file_to_crosssection[short_name]

                                                #print(prescale)

                                                
                                                for mod_hist in all_hists[('zg_10', 'rg_10')]:
                                                        hist = mod_hist.hist()
                                                        conditions = mod_hist.conditions()

                                                        try:
                                                                condition_satisfied = 1
                                                                for condition_keyword, condition_boundaries in conditions:
                                                                        keyword_index = keywords.index(condition_keyword) + 1

                                                                        if condition_boundaries[0] == None and condition_boundaries[1] != None:
                                                                                condition_satisfied *= int( float(numbers[keyword_index]) < condition_boundaries[1] ) 
                                                                        elif condition_boundaries[0] != None and condition_boundaries[1] == None:
                                                                                condition_satisfied *= int( float(numbers[keyword_index]) > condition_boundaries[0] ) 
                                                                        elif condition_boundaries[0] == None and condition_boundaries[1] == None:
                                                                                condition_satisfied *= 1 
                                                                        elif condition_boundaries[0] != None and condition_boundaries[1] != None:
                                                                                condition_satisfied *= int( float(numbers[keyword_index]) > condition_boundaries[0] and float(numbers[keyword_index]) < condition_boundaries[1] )

                                                                condition_satisfied = bool(condition_satisfied)
                                                        except Exception as e:
                                                                print "ASF", e
                                                        

                                                        if condition_satisfied:  

                                                                hist.fill_array( np.array([[zg_10, rg_10]]), [prescale] )



                                                for mod_hist in log_hists[('zg_10', 'rg_10')]:
                                                        hist = mod_hist.hist()
                                                        conditions = mod_hist.conditions()

                                                        try:
                                                                condition_satisfied = 1
                                                                for condition_keyword, condition_boundaries in conditions:
                                                                        keyword_index = keywords.index(condition_keyword) + 1

                                                                        if condition_boundaries[0] == None and condition_boundaries[1] != None:
                                                                                condition_satisfied *= int( float(numbers[keyword_index]) < condition_boundaries[1] ) 
                                                                        elif condition_boundaries[0] != None and condition_boundaries[1] == None:
                                                                                condition_satisfied *= int( float(numbers[keyword_index]) > condition_boundaries[0] ) 
                                                                        elif condition_boundaries[0] == None and condition_boundaries[1] == None:
                                                                                condition_satisfied *= 1 
                                                                        elif condition_boundaries[0] != None and condition_boundaries[1] != None:
                                                                                condition_satisfied *= int( float(numbers[keyword_index]) > condition_boundaries[0] and float(numbers[keyword_index]) < condition_boundaries[1] )

                                                                condition_satisfied = bool(condition_satisfied)
                                                        except Exception as e:
                                                                print "ASF", e
                                                        

                                                        if condition_satisfied:        

                                                                hist.fill_array( np.array([[zg_10, rg_10]]), [prescale] )


                                                for mod_hist in all_hists[('track_zg_10', 'track_rg_10')]:
                                                        hist = mod_hist.hist()
                                                        conditions = mod_hist.conditions()

                                                        try:
                                                                condition_satisfied = 1
                                                                for condition_keyword, condition_boundaries in conditions:
                                                                        keyword_index = keywords.index(condition_keyword) + 1

                                                                        if condition_boundaries[0] == None and condition_boundaries[1] != None:
                                                                                condition_satisfied *= int( float(numbers[keyword_index]) < condition_boundaries[1] ) 
                                                                        elif condition_boundaries[0] != None and condition_boundaries[1] == None:
                                                                                condition_satisfied *= int( float(numbers[keyword_index]) > condition_boundaries[0] ) 
                                                                        elif condition_boundaries[0] == None and condition_boundaries[1] == None:
                                                                                condition_satisfied *= 1 
                                                                        elif condition_boundaries[0] != None and condition_boundaries[1] != None:
                                                                                condition_satisfied *= int( float(numbers[keyword_index]) > condition_boundaries[0] and float(numbers[keyword_index]) < condition_boundaries[1] )

                                                                condition_satisfied = bool(condition_satisfied)
                                                        except Exception as e:
                                                                print "ASF", e
                                                        

                                                        if condition_satisfied:

                                                                hist.fill_array( np.array([[track_zg_10, track_rg_10]]), [prescale] )



                                                for mod_hist in log_hists[('track_zg_10', 'track_rg_10')]:
                                                        hist = mod_hist.hist()
                                                        conditions = mod_hist.conditions()

                                                        try:
                                                                condition_satisfied = 1
                                                                for condition_keyword, condition_boundaries in conditions:
                                                                        keyword_index = keywords.index(condition_keyword) + 1

                                                                        if condition_boundaries[0] == None and condition_boundaries[1] != None:
                                                                                condition_satisfied *= int( float(numbers[keyword_index]) < condition_boundaries[1] ) 
                                                                        elif condition_boundaries[0] != None and condition_boundaries[1] == None:
                                                                                condition_satisfied *= int( float(numbers[keyword_index]) > condition_boundaries[0] ) 
                                                                        elif condition_boundaries[0] == None and condition_boundaries[1] == None:
                                                                                condition_satisfied *= 1 
                                                                        elif condition_boundaries[0] != None and condition_boundaries[1] != None:
                                                                                condition_satisfied *= int( float(numbers[keyword_index]) > condition_boundaries[0] and float(numbers[keyword_index]) < condition_boundaries[1] )

                                                                condition_satisfied = bool(condition_satisfied)
                                                        except Exception as e:
                                                                print "ASF", e
                                                        

                                                        if condition_satisfied:
                                                                        
                                                                hist.fill_array( np.array([[track_zg_10, track_rg_10]]), [prescale] )	

                                except Exception as e:
                                        print "Some exception occured!",
                                        print e


	return all_hists, log_hists



def parse_to_root_file(input_files, output_filename, hist_templates):

    # print "Parsing {} to {}".format(input_filename, output_filename)

    # First, get the already histogrammed hists.

    hists, log_hists = copy.deepcopy(hist_templates[0]), copy.deepcopy(hist_templates[1])

    parsed_hists, parsed_log_hists = parse_file(input_files, hists, log_hists)

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

	hist_templates = hists.two_dim_hists()
	log_hist_templates = hists.two_dim_log_hists()

	parse_to_root_file(input_files=data_files_2011, output_filename=(output_directory + "2d_data.root", output_directory + "2d_data_log.root"), hist_templates=(hist_templates, log_hist_templates))


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

def load_root_files_to_hist(log=False):
	
	# parse_to_root_files()

	output_directory = sys.argv[1]

	if not log:
		hist_templates = hists.two_dim_hists()
		filenames = ["2d_data.root"]
	else:
		hist_templates = hists.two_dim_log_hists()
		filenames = ["2d_data_log.root"]

	return  [ root_file_to_hist(output_directory + filename, hist_templates) for filename in filenames ] 


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

	# load_root_files_to_hist()

	# parse_pfc_to_root_files()

	pass
