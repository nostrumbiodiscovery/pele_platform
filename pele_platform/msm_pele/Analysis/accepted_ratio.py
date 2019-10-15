import os
import glob
import pandas as pd
import numpy as np

FILES = glob.glob("report_*")

acceptances_rates = np.empty(len(FILES))
total_steps_l = np.empty(len(FILES))
accepted_steps_l  = np.empty(len(FILES))

for i, f in enumerate(FILES):
    data = pd.read_csv(f, sep='    ', engine='python')
    accepted_steps = float(data.iat[-1, 2])
    total_steps = float(data.iat[-1, 1])
    acceptance_ratio = accepted_steps/total_steps*100
    total_steps_l[i] = total_steps
    accepted_steps_l[i] = accepted_steps
    acceptances_rates[i] = acceptance_ratio
print("Steps made {} +- {}".format(np.mean(total_steps_l), np.std(total_steps_l)))
print("Steps accepted {} +- {}".format(np.mean(accepted_steps_l), np.std(accepted_steps_l)))
print("Acceptance ratio {} +- {}".format(np.mean(acceptances_rates), np.std(acceptances_rates)))
    

        
