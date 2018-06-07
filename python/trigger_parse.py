from __future__ import division

import math
import time
import sys
import hists


from MODPlot import *


from rootpy.io import File as TFile



import matplotlib.pyplot as plt
import rootpy.plotting.root2matplotlib as rplt





def parse_file(input_files, all_hists):

	
	# We read the file line by line, and for each line, we fill the corresponding histograms.
	# This is desirable to creating lists of values since this will not hold anything in memory. 

	keywords = []
	keywords_set = False

	print(input_files)

	for input_file in input_files:

                print(input_file)
                trigger_luminosity_dictionary = trigger_luminosity_dictionary_2011()


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


                                                pT_index = keywords.index("corr_hardest_pT") + 1
                                                trig_name_index = keywords.index("trigger_name") + 1
                                                eta_index = keywords.index("hardest_eta") + 1

                                                #print(pT_index)
                                                pT_of_this_event = float(numbers[pT_index])
                                                eta_of_this_event = float(numbers[eta_index])
                                                current_line_trig_name = numbers[trig_name_index]
                                                #trigger = numbers[keywords_index_dictionary['trigger_fired']]

                                                #print(current_line_trig_name)

                                                if abs(eta_of_this_event) > 2.5:
                                                        # print "oops, pseudorapidity too large.", eta_of_this_event
                                                        continue

                                                #print(pT_of_this_event, prescale)

                                                #print(all_hists.keys())

                                                #print(pT_of_this_event, prescale)

                                                        #print(current_line_trig_name)

                                                #print(trigger_luminosity_dictionary)

                                                prescale = 1000000.0/trigger_luminosity_dictionary[current_line_trig_name]
                                                #print(current_line_trig_name, prescale)

                                                mod_hist = all_hists[current_line_trig_name]
                                                hist = mod_hist.hist()


                                                conditions = mod_hist.conditions()

                                                condition_satisfied = True

                                                #print(conditions)

                                                for condition_keyword, condition_func in conditions:
                                                        keyword_index = keyword_index = keywords.index(condition_keyword[0]) + 1
                                                        condition_func_param = numbers[keyword_index]

                                                        if (not condition_func(condition_keyword[1], int(condition_func_param))):
                                                                condition_satisfied = False


                                                if condition_satisfied:
                                                        all_hists[current_line_trig_name].hist().fill_array([pT_of_this_event], [prescale])

                                                """

                                                for key in all_hists.keys():
						
                                                        # See whether condition/s are satisfied or not.
                                                        mod_hist = all_hists[key]
                                                	hist = mod_hist.hist()


                                                        conditions = mod_hist.conditions()

                                                        try:
                                                                condition_satisfied = True

                                                                for condition_keyword, condition_func in conditions:
                                                                        keyword_index = keyword_index = keywords.index(condition_keyword[0]) + 1
                                                                        condition_func_param = numbers[keyword_index]

                                                                        if (not condition_func(condition_keyword[1], int(condition_func_param))):
                                                                                condition_satisfied = False

                                                        except Exception as excp:
                                                                print "Some exception occured while processing conditions.", excp

                                                        print(condition_satisfied)

                                                """


                                                """

                                                for key in all_hists.keys():
                                                        
                                                        # See whether condition/s are satisfied or not.
                                                        mod_hist = all_hists[key]
                                                        hist = mod_hist.hist()

                                                        print(key)

                                                        print(current_line_trig_name)

                                                        if 'NoJetID' not in current_line_trig_name and (key+"_" in current_line_trig_name or ("prescale" in key and key.split("_prescale")[0]+"_" in current_line_trig_name)):

                                                                conditions = mod_hist.conditions()

                                                                try:
                                                                        condition_satisfied = True

                                                                        for condition_keyword, condition_func in conditions:
                                                                                keyword_index = keyword_index = keywords.index(condition_keyword[0]) + 1
                                                                                condition_func_param = numbers[keyword_index]

                                                                                if (not condition_func(condition_keyword[1], int(condition_func_param))):
                                                                                        condition_satisfied = False

                                                                except Exception as excp:
                                                                        print "Some exception occured while processing conditions.", excp

                                                                x = pT_of_this_event
                                                                y = prescale

                                                                print("condition satisfied")

                                                                print("prescale_and_x" , x, prescale)
                                                                
                                                                all_hists[key].hist().fill_array([x], [y])

                                                """

                                except Exception as exc:
                                        pass
                                        #print "Some exception occured.", exc
                                        #print line
                                        #print "\n"
                                                

        #print(all_hists['HLT_Jet190'])
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

	parse_to_root_file(input_filename=data_files_2011, output_filename=output_directory + "trig.root", hist_templates=hist_templates)

def load_root_files_to_hist(log=False):

        output_directory = sys.argv[1]
	
	hist_templates = hists.trigger_hists()

	filenames = ["trig.root"]

	return  [ root_file_to_hist(output_directory + filename, hist_templates) for filename in filenames ] 


def trigger_luminosity_dictionary_2011():
    mod_file_with_trigger = [line.rstrip('\n') for line in open("/home/preksha/Documents/mengproject/MODAnalyzer/effective_luminosity_by_trigger.csv")]
    mod_trigger_luminosities = {}
    for mod_trigger_lumi in mod_file_with_trigger:
        mod_trigger_lumi.replace(" ","")
        if not mod_trigger_lumi == "":
            #print(mod_trigger_lumi)
            mod_trigger_luminosities[(mod_trigger_lumi.split(',')[0], mod_trigger_lumi.split(',')[1])] = mod_trigger_lumi.split(',')[2]
    trigger_luminosity_total = {}
    for mod_trigger in mod_trigger_luminosities:
        if mod_trigger[1] in trigger_luminosity_total:
            #current_effective_lumi = trigger_luminosity_total[mod_trigger[1]]
            trigger_luminosity_total[mod_trigger[1]] += float(mod_trigger_luminosities[mod_trigger])
        else:
            trigger_luminosity_total[mod_trigger[1]] = float(mod_trigger_luminosities[mod_trigger])
    return trigger_luminosity_total



if __name__ == "__main__":

	
	parse_to_root_files()

	#load_root_files_to_hist()

	pass
