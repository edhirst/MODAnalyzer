import os

mod_path = '/Users/prekshanaik/Documents/MEng_Project'

all_files = []

for filename in os.listdir(mod_path):
    if filename.endswith('.mod'):
        all_files.append(mod_path + '/' + filename)

for filename in all_files:
    current_events = [line.rstrip('\n') for line in open(filename)]
    ak5count = 0
    PFCcount = 0
    for line in current_events:
        if 'Trig' in line:
            split_line = line.split()
        ak5count = 0
        if 'AK5' in line:
            ak5count += 1
        if 'PFC' in line:
            PFCcount += 1
        print(ak5count, PFCcount)
