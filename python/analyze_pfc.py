from subprocess import call
import os
from time import time
import sys
from collections import defaultdict

input_path = sys.argv[1]
output_path = sys.argv[2]
type_of_data = sys.argv[3]

def run_analyzer(input_path, output_path, type_of_data):
  to_analyze_and_output = []
  if type_of_data == "2011":
    for file in os.listdir(input_path):
      if file.endswith("mod"):
        to_analyze_and_output.append((input_path + "/" + file, output_path + "/" + file.replace(".mod", ".dat")))
  elif type_of_data == "sim":
    all_pythia_runs = os.listdir(input_path)
    for pythia_run in all_pythia_runs:
        print(pythia_run)
        for dirName, subdirList, files in os.walk(input_path + pythia_run):
            for file in files:
                if file.endswith("mod"):
                    to_analyze_and_output.append((dirName + "/" + file, output_path + "/" + pythia_run + file.replace(".mod", ".dat")))
  for input, output in to_analyze_and_output:
    args = ['./bin/analyze_pfc', input, output ] + sys.argv[3:]
    call(args)

start = time()

run_analyzer(input_path, output_path, type_of_data)

end = time()

print "Everything done in " + str(end - start) + " seconds!"
