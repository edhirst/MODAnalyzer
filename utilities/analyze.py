from subprocess import call
import os
from time import time
import sys
from collections import defaultdict

input_path = sys.argv[1]
output_path = sys.argv[2]
type_of_data = sys.argv[3]



def get_existing_files(output_path):
    existing_files = {}
    for path, subdirs, files in os.walk(output_path):
        for name in files:
            existing_files[(os.path.join(path, name)).split("/")[-1]] = 1
    return existing_files



def run_analyzer(input_path, output_path, type_of_data):
  existing_files = get_existing_files(output_path)
  to_analyze_and_output = []
  if type_of_data == "2011":
    for file in os.listdir(input_path):
      if file.endswith("mod"):
        if file.replace(".mod", ".dat") not in existing_files:
          to_analyze_and_output.append((input_path + "/" + file, output_path + "/" + file.replace(".mod", ".dat")))
  elif type_of_data == "sim":
    all_pythia_runs = os.listdir(input_path)
    for pythia_run in all_pythia_runs:
        print(pythia_run)
        for dirName, subdirList, files in os.walk(input_path + pythia_run):
            for file in files:
                if file.endswith("mod"):
                    if (pythia_run + file.replace(".mod", "_sim_pfc.dat")) not in existing_files:
                        to_analyze_and_output.append((dirName + "/" + file, output_path + "/" + pythia_run + file.replace(".mod", ".dat")))
  for input, output in to_analyze_and_output:
    args = ['./bin/analyze', input, output ] + sys.argv[3:]
    call(args)

start = time()

run_analyzer(input_path, output_path, type_of_data)

end = time()

print "Everything done in " + str(end - start) + " seconds!"
