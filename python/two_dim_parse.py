from __future__ import division

import math
import time
import sys
import hists


from MODPlot import *


from rootpy.io import File as TFile



import matplotlib.pyplot as plt
import rootpy.plotting.root2matplotlib as rplt



output_directory = "/home/preksha/Documents/mengproject/"

data_file = "/home/preksha/Documents/mengproject/analyzed_data10.dat"



average_prescales = {}


average_prescales[(250, None)] = 1.0
average_prescales[(200, 250)] = 1.933420103
average_prescales[(150, 200)] = 5.361922609
average_prescales[(115, 150)] =  100.3122906
average_prescales[(85, 115)] =  851.3943491


def parse_file(input_file, all_hists, log_hists):

	print "Parsing {}".format(input_file)
	
	# We read the file line by line, and for each line, we fill the corresponding histograms.
	# This is desirable to creating lists of values since this will not hold anything in memory. 

	keywords = []
	keywords_set = False

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

					prescale_index = keywords.index("prescale") + 1
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
					prescale = float(numbers[prescale_index])
					eta = float(numbers[eta_index])


					
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
								
							if input_file == data_file:	# For data file only.

								if not mod_hist.use_prescale():
									hist.fill_array( np.array([[zg_10, rg_10]]) )
								else:
									
									#  Average prescale. 

									# To find which prescale to use, we need to find which trigger fired. 
									# To do that, we need to find the pT of the hardest jet.
									
									prescale_to_use = 0.0

									if pT_of_this_event > 250.:
										prescale_to_use = 1.0
									else:
										for pT_boundaries, prescale in average_prescales.items():
											
											lower, upper = pT_boundaries

											if upper != None:
												if pT_of_this_event > float(lower) and pT_of_this_event < float(upper):
													prescale_to_use = prescale
													break

									# hist.fill_array( np.array([[zg_10, rg_10]]), [prescale_to_use] )	 
									hist.fill_array( np.array([[zg_10, rg_10]]), [float(numbers[prescale_index])] )	 
							
							else:	# MC so always use prescales.

								hist.fill_array( np.array([[zg_10, rg_10]]), [float(numbers[prescale_index])] )



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
								
							if input_file == data_file:	# For data file only.

								if not mod_hist.use_prescale():
									hist.fill_array( np.array([[zg_10, rg_10]]) )
								else:
									
									#  Average prescale. 

									# To find which prescale to use, we need to find which trigger fired. 
									# To do that, we need to find the pT of the hardest jet.
									
									prescale_to_use = 0.0

									if pT_of_this_event > 250.:
										prescale_to_use = 1.0
									else:
										for pT_boundaries, prescale in average_prescales.items():
											
											lower, upper = pT_boundaries

											if upper != None:
												if pT_of_this_event > float(lower) and pT_of_this_event < float(upper):
													prescale_to_use = prescale
													break

									# hist.fill_array( np.array([[zg_10, rg_10]]), [prescale_to_use] )	 
									hist.fill_array( np.array([[zg_10, rg_10]]), [float(numbers[prescale_index])] )	 
							
							else:	# MC so always use prescales.

								hist.fill_array( np.array([[zg_10, rg_10]]), [float(numbers[prescale_index])] )


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
								
							if input_file == data_file:	# For data file only.

								if not mod_hist.use_prescale():
									hist.fill_array( np.array([[track_zg_10, track_rg_10]]) )
								else:
									
									#  Average prescale. 

									# To find which prescale to use, we need to find which trigger fired. 
									# To do that, we need to find the pT of the hardest jet.
									
									prescale_to_use = 0.0

									if pT_of_this_event > 250.:
										prescale_to_use = 1.0
									else:
										for pT_boundaries, prescale in average_prescales.items():
											
											lower, upper = pT_boundaries

											if upper != None:
												if pT_of_this_event > float(lower) and pT_of_this_event < float(upper):
													prescale_to_use = prescale
													break

									# hist.fill_array( np.array([[zg_10, rg_10]]), [prescale_to_use] )	 
									hist.fill_array( np.array([[track_zg_10, track_rg_10]]), [float(numbers[prescale_index])] )	 
							
							else:	# MC so always use prescales.

								hist.fill_array( np.array([[track_zg_10, track_rg_10]]), [float(numbers[prescale_index])] )



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
								
							if input_file == data_file:	# For data file only.

								if not mod_hist.use_prescale():
									hist.fill_array( np.array([[track_zg_10, track_rg_10]]) )
								else:
									
									#  Average prescale. 

									# To find which prescale to use, we need to find which trigger fired. 
									# To do that, we need to find the pT of the hardest jet.
									
									prescale_to_use = 0.0

									if pT_of_this_event > 250.:
										prescale_to_use = 1.0
									else:
										for pT_boundaries, prescale in average_prescales.items():
											
											lower, upper = pT_boundaries

											if upper != None:
												if pT_of_this_event > float(lower) and pT_of_this_event < float(upper):
													prescale_to_use = prescale
													break

									# hist.fill_array( np.array([[zg_10, rg_10]]), [prescale_to_use] )	 
									hist.fill_array( np.array([[track_zg_10, track_rg_10]]), [float(numbers[prescale_index])] )	 
							
							else:	# MC so always use prescales.

								hist.fill_array( np.array([[track_zg_10, track_rg_10]]), [float(numbers[prescale_index])] )	

			except Exception as e:
				print "Some exception occured!",
				print e


	return all_hists, log_hists



def parse_to_root_file(input_filename, output_filename, hist_templates):

	print "Parsing {} to {}".format(input_filename, output_filename)
	
	parsed_hists, parsed_log_hists = parse_file( input_filename, copy.deepcopy( hist_templates[0] ), copy.deepcopy(hist_templates[1]) )

	f = TFile(output_filename[0], "RECREATE")

	for var in parsed_hists.keys():
		
		index = 0

		for mod_hist in parsed_hists[var]:
			hist = copy.deepcopy( mod_hist.hist() )
			hist.SetName("{}#{}".format(var, index))
			hist.Write()

			index += 1

	f.Close()


	f = TFile(output_filename[1], "RECREATE")

	for var in parsed_log_hists.keys():
		
		index = 0

		for mod_hist in parsed_log_hists[var]:
			hist = copy.deepcopy( mod_hist.hist() )
			hist.SetName("{}#{}".format(var, index))
			hist.Write()

			index += 1

	f.Close()


def root_file_to_hist(input_filename, hist_templates):

	hists = copy.deepcopy( hist_templates )

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
	hist_templates = hists.two_dim_hists()
	log_hist_templates = hists.two_dim_log_hists()

	parse_to_root_file(input_filename=data_file, output_filename=(output_directory + "2d_data.root", output_directory + "2d_data_log.root"), hist_templates=(hist_templates, log_hist_templates))
	# parse_to_root_file(input_filename=pythia_file, output_filename=(output_directory + "2d_pythia.root", output_directory + "2d_pythia_log.root"), hist_templates=(hist_templates, log_hist_templates))
	# parse_to_root_file(input_filename=herwig_file, output_filename=(output_directory + "2d_herwig.root", output_directory + "2d_herwig_log.root"), hist_templates=(hist_templates, log_hist_templates))
	# parse_to_root_file(input_filename=sherpa_file, output_filename=(output_directory + "2d_sherpa.root", output_directory + "2d_sherpa_log.root"), hist_templates=(hist_templates, log_hist_templates))

	# parse_to_root_file(input_filename=data_file, output_filename=output_directory + "2d_data_log.root", hist_templates=log_hist_templates)
	# parse_to_root_file(input_filename=pythia_file, output_filename=output_directory + "2d_pythia_log.root", hist_templates=log_hist_templates)
	# parse_to_root_file(input_filename=herwig_file, output_filename=output_directory + "2d_herwig_log.root", hist_templates=log_hist_templates)
	# parse_to_root_file(input_filename=sherpa_file, output_filename=output_directory + "2d_sherpa_log.root", hist_templates=log_hist_templates)


def load_root_files_to_hist(log=False):
	
	# parse_to_root_files()

	if not log:
		hist_templates = hists.two_dim_hists()
		filenames = ["2d_data.root"]
	else:
		hist_templates = hists.two_dim_log_hists()
		filenames = ["2d_data_log.root"]

	return  [ root_file_to_hist(output_directory + filename, hist_templates) for filename in filenames ] 




if __name__ == "__main__":

	parse_to_root_files()

	# load_root_files_to_hist()

	# parse_pfc_to_root_files()

	pass
