from __future__ import division

import math
import time
import sys
import hists


from MODPlot import *


from rootpy.io import File as TFile



import matplotlib.pyplot as plt
import rootpy.plotting.root2matplotlib as rplt

satisfied_sets = {}



def parse_file(input_files, all_hists):

	
	# We read the file line by line, and for each line, we fill the corresponding histograms.
	# This is desirable to creating lists of values since this will not hold anything in memory. 

	keywords = []
	keywords_set = False

	print(input_files)

	for input_file in input_files:

                print(input_file)

                mod_file_to_crosssection = mod_file_to_crosssection_sim("sim_pfc", input_files)

                #print(mod_file_to_crosssection)

                with open(input_file) as infile:
                        
                        line_number = 0

                        current_line_trig_name = ""

                        for line in infile:


                                # if line_number > 100000:	# Ideal length.
                                # if line_number > 100000000:	# Very Huge.
                                # if line_number > 10000000:	# Huge.
                                # if line_number > 1000000:	# Big enough.
                                # if line_number > 1000:		# Small tests.
                                # if line_number > 30000:		# Small tests.
                                if False:
                                        break

                                line_number += 1

                                #print "At line number {}".format(line_number)

                                try:
                                        numbers = line.split()

                                        if len(numbers) > 0 and numbers[0] == "#" and (not keywords_set):
                                                keywords = numbers[2:]
                                                keywords_set = True

                                        elif len(numbers) > 0 and numbers[0] == "Entry":

                                                #print(line)


                                                pT_index = keywords.index("hardest_pT") + 1
                                                #trig_name_index = keywords.index("trigger_name") + 1
                                                eta_index = keywords.index("hardest_eta") + 1

                                                #print(pT_index)
                                                pT_of_this_event = float(numbers[pT_index])
                                                #print(pT_of_this_event)
                                                eta_of_this_event = float(numbers[eta_index])
                                                #current_line_trig_name = numbers[trig_name_index]
                                                #trigger = numbers[keywords_index_dictionary['trigger_fired']]


                                                if abs(eta_of_this_event) > 2.4:
                                                        continue
                                                
                                                short_name = input_file.split('/')[-1].replace('sim_pfc', 'sim_gen')
                                                prescale_to_use = mod_file_to_crosssection[short_name]


                                                mod_name_start_index = short_name.find("pythia6") + len("pythia6")
                                                pythia_set = short_name[:mod_name_start_index]

                                                mod_hist = all_hists[pythia_set]
                                                hist = mod_hist.hist()


                                                conditions = mod_hist.conditions()

                                                condition_satisfied = True


                                                for condition_keyword, condition_func in conditions:
                                                        keyword_index = keyword_index = keywords.index(condition_keyword[0]) + 1
                                                        condition_func_param = numbers[keyword_index]

                                                        if (not condition_func(condition_keyword[1], int(condition_func_param))):
                                                                condition_satisfied = False


                                                if condition_satisfied:
                                                        #print(pythia_set, pT_of_this_event, prescale_to_use)
                                                        satisfied_sets[pythia_set] = 1
                                                        all_hists[pythia_set].hist().fill_array([pT_of_this_event], [prescale_to_use])


                                except Exception as exc:
                                        pass
                                        #print "Some exception occured.", exc
                                        #print line
                                        #print "\n"
                                                

        #print(all_hists['HLT_Jet190'])
        print(satisfied_sets)
        return all_hists



def parse_to_root_file(input_filename, output_filename, hist_templates):

	print "Parsing {} to {}".format(input_filename, output_filename)
	
	parsed_hists = parse_file( input_filename, copy.deepcopy( hist_templates ) )

	f = TFile(output_filename, "RECREATE")

	for var in parsed_hists.keys():
		mod_hist = parsed_hists[var]
		hist = copy.deepcopy( mod_hist.hist() )
		hist.SetName("{}".format(var))
		hist.Write()

	f.Close()


def root_file_to_hist(input_filename, hist_templates):

	hists = copy.deepcopy( hist_templates )

	root_file = TFile(input_filename, "read")

	for var in hists.keys():
		mod_hist = hists[var]
		hist = root_file.Get(var)
		mod_hist.replace_hist(hist)

		
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

	hist_templates = hists.trigger_hists()

	parse_to_root_file(input_filename=data_files_sim_pfc, output_filename=output_directory + "trig_pfc.root", hist_templates=hist_templates)

def load_root_files_to_hist(log=False):

        output_directory = sys.argv[1]
	
	hist_templates = hists.trigger_hists()

	filenames = ["trig_pfc.root"]

	return  [ root_file_to_hist(output_directory + filename, hist_templates) for filename in filenames ] 

def trigger_luminosity_dictionary_2011(input_files):
    input_files_short = []
    for input_file in input_files:
        input_files_short.append(input_file.split('/')[-1])
    print(input_files_short)
    mod_file_with_trigger = [line.rstrip('\n') for line in open("effective_luminosity_by_trigger.csv")]
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
        input_files_short.append(input_file.split('/')[-1].replace('sim_pfc', 'sim_gen'))
    output_file_event_count = [line.rstrip('\n') for line in open("event_count_by_pythia_and_mod.csv")]
    pythia_cross_sections = {'QCD_Pt-0to5_TuneZ2_7TeV_pythia6': 48444950000.0,
                             'QCD_Pt-5to15_TuneZ2_7TeV_pythia6': 36745720000.0,
                             'QCD_Pt-15to30_TuneZ2_7TeV_pythia6': 815912800.0, 
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
                             'QCD_Pt-1400to1800_TuneZ2_7TeV_pythia6': 0.0108721,
                             'QCD_Pt-1800_TuneZ2_7TeV_pythia6': 0.000357463000}

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

	#load_root_files_to_hist()

	pass
