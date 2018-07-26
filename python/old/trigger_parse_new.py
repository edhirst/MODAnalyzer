import matplotlib.pyplot as plt

vals = ['#', 'Entry', 'event_number', 'run_number', 'trig_jet_matched', 'jet_quality', 'hardest_pT', 'corr_hardest_pT', 'hardest_eta', 'prescale', 'trigger_name']

trigger_values_pT = {}

print(trigger_values_pT)
with open('trigger_out_test') as f:
    lines = f.readlines()
    print(lines[0].split())
    for line in lines[1:]:
        current = line.split()
        if current[0] != '#':
            trigger_name = current[-1]
            pT_value = current[6]
            hardest_eta = current[7]
            if abs(float(hardest_eta)) <= 2.4:
                if trigger_name not in trigger_values_pT:
                    trigger_values_pT[trigger_name] = []
                else:
                    current = trigger_values_pT[trigger_name]
                    currnet = current.append(float(pT_value))
                    trigger_values_pT[trigger_name] = current
print(sorted(trigger_values_pT['HLT_Jet70U']))

n, bins, patches = plt.hist(sorted(trigger_values_pT['HLT_Jet70U'])[:-10], 50, facecolor='g')


plt.xlabel('Trigger_jet_pT')
plt.ylabel('Events')
plt.grid(True)
plt.show()

