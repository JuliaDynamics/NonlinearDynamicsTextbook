"""
This is a python file. Why? Because I found no `Julia` package to read `.edf` files
and using `pyedflib` via `PythonCall` doesn't work very nicely.

This file is basically just converting an EEG recording of a patient with an epileptic
seizure from a data set containing such EEG recordings [1, 2] available at `physionet` [3]
into a `.csv` file so that you can  nicely read it in.

All needed python packages are found in this directory in `python_requirements.txt`.
I am using python 3.8.10.

## References

[1] Detti, P. (2020). Siena Scalp EEG Database (version 1.0.0). PhysioNet. https://doi.org/10.13026/5d4a-j060.

[2] Detti, Paolo, Giampaolo Vatti, and Garazi Zabalo Manrique de Lara. “EEG Synchronization Analysis for Seizure Prediction: A Study on Data of Noninvasive Recordings.” Processes 8, no. 7 (July 16, 2020): 846. https://doi.org/10.3390/pr8070846.

[3] Goldberger, A., Amaral, L., Glass, L., Hausdorff, J., Ivanov, P. C., Mark, R., ... & Stanley, H. E. (2000). PhysioBank, PhysioToolkit, and PhysioNet: Components of a new research resource for complex physiologic signals. Circulation [Online]. 101 (23), pp. e215–e220.

"""
import os

import numpy as np
import wget
import pandas as pd
from pyedflib import highlevel


URL = 'https://physionet.org/files/siena-scalp-eeg/1.0.0/PN00/'

if not os.path.exists('Seizures-list-PN00.txt'):
        r = wget.download(URL+'Seizures-list-PN00.txt')

# we will first read out all the relevant time stamps from the info file
# (we actually only need the first one but it's easier to just iterate over all)
reg_start_times = []
seizure_start_times = []
seizure_end_times = []
with open('Seizures-list-PN00.txt', 'r', encoding='UTF-8') as f:
    for line in f.readlines():
        if 'Registration start' in line:
            reg_start_times.append(pd.to_datetime(line.split()[-1], format='%H.%M.%S'))
        elif 'Seizure start' in line:
            seizure_start_times.append(pd.to_datetime(line.split()[-1], format='%H.%M.%S'))
        elif 'Seizure end' in line:
            seizure_end_times.append(pd.to_datetime(line.split()[-1], format='%H.%M.%S'))

# we now fill two dictionaries with the time series and the time stamps

file = 'PN00-1.edf'
if not os.path.exists(file):
    r = wget.download(URL+file)
    print(r)
# We get channel F8 as one of the right hemisphere channels
signal, *_ = highlevel.read_edf(file, ch_names='EEG F8')
time_stamps = [
    (seizure_start_times[0] - reg_start_times[0]).total_seconds(),
    (seizure_end_times[0] - reg_start_times[0]).total_seconds()
]
signal = signal.flatten()
# here, we create a list of time stamps of each point in the time series
time = np.arange(0, len(signal)/512, 1/512)

print(f'seizure started at {time_stamps[0]}s and ended at {time_stamps[1]}s')

result = np.vstack([time, signal])
np.savetxt('exercise_data/18.csv', result.T)