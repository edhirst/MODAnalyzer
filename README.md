# MIT Open Data Analyzer

This repository has code to analyze MOD (MIT Open Data) files produced using [MODProducer](https://github.com/rmastand/MODProducer "MODProducer" for latest version).

## Workflow

We adopt the following workflow for extracting MOD files out of the given AOD files.

1. Download all the ROOT files and arrange them in the same directory structure as they live on the CMS server.

2. Create a registry that maps each event and run number to a certain ROOT file. This is done so that things can be done one file at a time as we're dealing with multiple TeraBytes of data here and it's a pain to have to do everything at once. 

3. Run the Producer on those AOD files. This reads the download directory and processes only the files in there. This produces N MOD files. 

4. Read in the output files one by one and perform calculation to produce the relevant dat files corresponding to each input file. In the case of 2011 data, the number of output files is 1. In the the case of simulated data, the number of output files is 2, corresponding to reconstructed and truth data.  

5. Produce plots using the dat files produced in step (4)

This repository is concerned with steps (4) to (5) only. Steps (1) to (3) are carried out by the [MODProducer](https://github.com/rmastand/MODProducer/ "MODProducer") package.

## Usage Instruction

### Preparation

 - Install [FastJet](http://www.fastjet.fr/ "FastJet") (in addition to fastjet-contrib). Note the path to FastJet config. 
 
 - Open `./Makefile` and paste the path to FastJet config to the **PATH\_TO\_FASTJET** variable on line 1. You will likely end up with a path ending with 'fastjet-install/bin/fastjet-config'
 
 - Install [Pythia](http://home.thep.lu.se/~torbjorn/pythia81html/Welcome.html "Pythia 8") Note the path to the pythia config. 
 
 - Open `./Makefile` and paste the path to Pythia config to the **PATH\_TO\_PYTHIA** variable on line 3. You will likely end up with a path ending with 'pythia8235/bin/pythia8-config'

 -  Note the path to the directory that contains the MOD files you've produced from the [MODProducer](https://github.com/rmastand/MODProducer/ "MODProducer") package.
 
 - Compile everything with `make`.

### Workflow Instructions

 - First, we run the analyzer. We use the Python script `utilities/analyze.py` for general data analysis.  This script will run the executable `bin/analyze` M times for M MOD files. This script has three required arguments and two optional arguments. 

   1. [required] If running 2011 data, a path to the directory of where all the MOD files. If running simulated data, a path to the directory where the QCD data sets are(e.g. 'All_QCD_data_sets' which contains 'QCD_Pt-15to30_TuneZ2_7TeV_pythia6', 'QCD_Pt-30to50_TuneZ2_7TeV_pythia6', etc.). Make sure the directory label corresponds to the pythia set. The files can be nested within the directory somewhere(e.g. 'QCD_Pt-30to50_TuneZ2_7TeV_pythia6 --> MODFILE --> SIMFiles --> x.mod') but make sure the top directory is labeled correctly.
   2. [required] Folder path to output the dat files into. 
   3. [required] '2011' if running 2011 data and 'sim' if running simulated data
   4. [optional] To limit number of arguments, use flag '-number_events'. By default, the analysis is run on all the events in the file. 
   5. [optional] If your 2011lumibyls.csv file (which is necessary to get luminosity information for 2011 data) lives somewhere else besides the directory you are running the command from, use flag '-2011_lumi_file' to specify its location. 
   
   ```
   python ./utilities/analyze.py /home/opendata/eos/opendata/cms/eos/opendata/cms/Run2011A/Jet/MOD/12Oct2013-v1/20000/ /home/prekshan/Documents/opendata_output/ 2011 -number_events 1000 -2011_lumi_file /home/prekshan/Downloads/2011lumibyls.csv
   ```
   The file skips analysis of anything that is already in the output directory even if partially downloaded so make sure to    delete those. For the simulated data, take care to delete both the sim_pfc version and the sim_gen versions of the file. 
   
- The output of running this script will be the following files:
   
   1. N or 2N files depending on whether you ran 2011 or simulated data
   2. If running 2011 data, you will get a 'effective_luminosity_by_jet.csv' file that will be used in the plotting. 
   3. If running simulated data, you will get a 'event_count_by_pythia_and_mod.csv' that will be used in the plotting. 
     
 - Finally, we are ready to produce plots. Note that the plotting framework makes use of [matplotlib](http://matplotlib.org/ "matplotlib") and [rootpy](http://rootpy.org/ "rootpy") so please make sure these are installed. 
 
To run the plotting, first run the Python script `python/parse.py`. This script takes two arguments:
 
   1. path to where all your dat files are located.
   2. path to an output directory where your plots will be produced 
   
  Make sure the csv files produced in the analysis from the previous step are in the same directory as where you are running your command from. 

     ```
     python ./python/parse.py /home/preksha/Documents/cms-opendata-dat/ /home/preksha/Documents/MODAnalyzer/plots/
    
     ```
   Next, to produce the actual plots, run the script `python/plots.py` using the path where you want your plots produced.  
   
     ```
     python ./python/plots.py /home/preksha/Documents/mengproject/MODAnalyzer/plots/
    
     ```
   
 Furthermore, to generate plots for additional variables, you can modify or add in create_multi_page_plot commands to `python/plots.py`. See `python/plots.py` for examples of linear, log, and PFC plot outputs. Plots currently get output into the plots/version 1 subdirectory. Change the default_dir in `python/plots.py` if you would like to modify this, making sure the directory exists. For whatever variable you are plotting, make sure in hists.py, the functions multi_page_log_hist_templates or multi_page_hist_templates contain an histogram for that variable. 

## TODO
- [ ] Fix plot formatting.
- [ ] Weighted MC events. 
- [ ] Luminosity conversion factor. Right now our total luminosity is ~55 pb-1 and that's with just ~40% of all data.
- [ ] Analytic Calculations.
- [ ] Run code on entire dataset.



####config_jthaler.mk:
- [x] It seems like you can delete this file, since you've set PATH_TO_FASTJET the same as I have on my computer.

####examples:
- [x] This folder contains both example analyses and essential tools.  You should probably call it "exec" instead of examples.

####interface:
- [x] All .hh/.cc files should be CamelCase capitalized.  Only .cc files for main routines should be lowercase with underscores.  (Basically, the name of the file should match the name of the class.

####mod_logo.eps:
- [ ] The bounding box seems to be a bit too small on the right-hand side, which will make the conversion to png a bit weird.  I can fix that if you don't know how.

####README.md:
- [ ] You should emphasize that the skimming workflow only works for single jet triggers from the Jet Primary Data set, assuming the user agrees with our trigger boundaries.

####src:
- [x] Same comment as for interface.  Should be CamelCase file names.

####analyze_data.cc:
- [ ] Line 67:  You should say how much time has elapsed and the number of seconds per event.  Do you really need to write this out every 100 events?  (You should have the time elapsed in all of the files.)
- [ ] Linee 111-112:  Does the order of entries matter?  If not, trigger name and prescale factor should come before jet pT information.
- [ ] This code only looks at the kinematics of the hardest jet, right?  So I would call this analyze_hardest_jet_kinematics.cc

####analyze_pair_pfc.cc:
- [ ] I guess we aren't going to use this, so can eventually go in obsolete.

####analyze_pfc.cc:
- [ ] Why aren't you checking that the assigned trigger fired, or using the hardest jet to trigger jet?  It seems that you are assuming that you are running this on the skimmed files, but you should not have to skim to make this analysis work.  Rather, skimmed files should have all of the flags set to give "true" as needed.  (Same comment for other files.)
- [ ] You define "jet_def_cambridge", but never use it.
- [ ] You should *always* output run number, event number, and name of the trigger used.  Otherwise, things are impossible to debug.
- [ ] I'm a bit confused about the workflow.  Does hardest_jet() have the JEC factors applied?  As far as I can tell, it does not.  If the JEC factor is applied, does it apply to the constituents as well?  I think you always need to output the original jet pT, the JEC factor, and the corrected jet pT to make sure you don't have an issue.  If I were you, I would always apply JEC factors as part of the analysis, even if the JEC factors are just 1 in the case of MC.
- [ ] The code only looks at the kinematics of the PFCs, so I would call it analyze_pfc_kinematics.cc.

####analyze.cc:
- [ ] I am going to stop looking at your code for the evening right now.  But I think you have the potential for lots of confusion with the JEC factor, especially with line 121 craziness which makes no sense to me.  It is unclear whether hardest_jet() means corrected or not, so you should probably be explicit about this and have separate corrected_hardest_jet() and uncorrected_hardest_jet() methods to deal with this.  You should also have a hardest_jet_JEC() method.
- [ ] The more confusing technical issue is whether multiplying the JEC factor affects the constituents.  I think it does not.  So on line 136, this means that "frac_pT_loss" gives you garbage at the moment, since hardest_jet.pt() has the JEC factor applied, but soft_drop( hardest_jet ).pt() does not, since soft_drop rebuilds the jet from the constituents.
- [ ] I think this means, in fact, that you should never multiply JEC times a PseudoJet, since it is ambiguous what it does.  Rather, you should only multiply JEC factors times doubles explicitly as needed.  This way, we always know exactly what we are doing.  I guess this means that you don't want a corrected_hardest_jet() method.  Rather, you should be doing a lot of your analysis on uncorrected jets, and then explicitly rescaling by the JEC factor later.
- [ ] Last comment, and then I'll stop.  You should make sure to read through the PRD draft at some point and get clarification on any definition that you are unfamiliar with.  In particular, your function for angularity_lambda does not agree with what we have in the text, since in the text, I say that you calculate the angle to a recoil-free axis.  But don't change it, since we might want to use your definition since it is simpler to explain.


####analyze.cc:
- [ ] Line 152 and 157:  Beware of magic numbers.  Define beta = 0.0 and R0 = 0.5 before using them.
- [ ] Line 164:  Instead of pow, define an inline sq() command that does multiplication rg*rg
- [x] Line 173:  We are no longer doing SoftKiller in this way, so get rid of that part of the code.
- [ ] Line 219:  What does "pT_after_SD" mean?  Does it include soft killer?  JEC factors?  Also, you keep recalling, soft_drop(hardest_jet), which is computationally expense.
- [x] In angularity_lambda(), do not hard code in the jet radius, since we might want to use different values later.  
- [ ] Also, I am now thinking that we should do WTA axes for define the angularities.  I can tell you how to do that in person.

####convert_to_pristine.cc:
- [x] Change name of this program to convert_to_one_jet().
- [ ] **Partially Done**.Similarly, get rid of all references to "pristine" everywhere in the code.

####luminosity.cc:
- [x] Change name to analyze_lumi.cc (or something like that)

####turn_on.cc:
- [x] Change name of this program to analyze_triggers.cc (or something like that)

####event.h:
- [ ] When you return a reference (i.e. with a &) it makes sense to label it as const.  But when you return an object, then best not to have the output be const.  So "int event_number() const;"  The function is const, but the int that is returned can be changed by the user if desired.
- [x] Why is closest_fastjet_jet_to_trigger_jet() not const?  It shouldn't change the event.
- [ ] Line 90:  I don't think that is valid C++ code using the old standard.  I'll see if it compiles on my machine, but the standard is to define these things in the constructor.

####event.cc:
- [ ] It seems to me that stringify_jet and stringify_pfc should really be part of InfoCalibratedJet and InfoPFC.  You are allowed to store a reference to the PseudoJet in those info classes so you can extract the px, py, etc information as needed to make the strings.
- [ ] Line 196:  Your comment about "pristine" makes me think that _data_type needs to have three strings.  That way the user can know if it is original, skimmed, or hardest jet.  Probably the key words will be something like "All_Events", "Matched_Trigger_Selection", "Hardest_Jet_Restriction"
- [ ] Line 359:  It seems like you treat the hardest jet case different from the ordinary case in terms of applying or not applying JEC factors.  As I said before, you should probably never rescale a PseudoJet by the JEC factor.  Only apply JECs to doubles.
- [ ] Line 394:  I'm a bit nervous about you having the weight in the "BeginEvent" line.  The data format of the first line should be identical between MC and data.  But this might have to wait for Version 6 if you already generated a ton of MC files with the weight put in the first line.
- [ ] In read_event(), if the MOD file is formatted incorrectly with a combination of RPFC and PDPFC, etc., then no error is thrown.  Rather, the data source is just changed.  The data source should be unambiguous from the "BeginEvent" line, and then cross checked for each particle input.  In fact, there really isn't a reason to have separate PDPFC, RPFC, etc., since they all do the same thing, once the data source is centrally determined.  This might have to wait until Version 6, though.
- [x] Line 561.  We don't use HLT_Jet15U or HLT_L1Jet6U in our analysis, so you shouldn't assign a trigger for those ranges.
- [x] Line 680.  Change 10 and -3 to doubles, since sometimes pow is defined on integers to return integers.
