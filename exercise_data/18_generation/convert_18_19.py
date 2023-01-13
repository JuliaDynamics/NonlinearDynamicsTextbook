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
import numpy as np
import wget
import pandas as pd
from pyedflib import highlevel
import os

URL = 'https://physionet.org/files/siena-scalp-eeg/1.0.0/PN00/'

downsampled_rate = 64
sampling_rate = 512

if not os.path.exists('Seizures-list-PN00.txt'):
        r = wget.download(URL+'Seizures-list-PN00.txt')

# we will first read out all the relevant time stamps from the info file
reg_start_times = []
seizure_start_times = []
seizure_end_times = []
with open('Seizures-list-PN00.txt', 'r') as f:
    for line in f.readlines():
        if 'Registration start' in line:
            reg_start_times.append(pd.to_datetime(line.split()[-1], format='%H.%M.%S'))
        elif 'Seizure start' in line:
            seizure_start_times.append(pd.to_datetime(line.split()[-1], format='%H.%M.%S'))
        elif 'Seizure end' in line:
            seizure_end_times.append(pd.to_datetime(line.split()[-1], format='%H.%M.%S'))

# we now fill two dictionaries with the time series and the time stamps
time_series_dict = {}
time_stamp_dict = {}
mx = -1e8
# we ignore recording three because it has an impossible seizure end documented (days after end of recording)
for recording in [1, 4]:
    file = f'PN00-{recording}.edf'
    if not os.path.exists(file):
        r = wget.download(URL+file)
        print(r)
    # We get channel F8 as one of the right hemisphere channels
    signal, *_ = highlevel.read_edf(file, ch_names='EEG F8')
    print(f'Seizure {recording} started at '
        f'{(seizure_start_times[recording-1] - reg_start_times[recording-1]).total_seconds()}s '
        f'index {(seizure_start_times[recording-1] - reg_start_times[recording-1]).total_seconds()*sampling_rate} '
        'and ended at '
        f'{(seizure_end_times[recording-1] - reg_start_times[recording-1]).total_seconds()}s '
        f'index {(seizure_end_times[recording-1] - reg_start_times[recording-1]).total_seconds()*sampling_rate}'
    )
    if mx < len(signal.flatten()):
        mx = len(signal.flatten())
    # downsample by taking only every sampling_rate//downsampled_rate step
    time_series_dict[f'{file}'] = signal.flatten()[::sampling_rate//downsampled_rate]


# to put them in a csv-File, we need to make them all the same length-> fill with NaN
for rec, signal in time_series_dict.items():
    a = np.full((mx//sampling_rate*downsampled_rate,), np.nan)
    a[:len(signal)] = signal
    time_series_dict[rec] = a

# we finally create a dataframe and save it as a csv-file
df = pd.DataFrame({'PN00-1': time_series_dict['PN00-1.edf']})
df.to_csv('exercise_data/19.csv', index=False, header=False, sep='\t')
df = pd.DataFrame({'PN00-4': time_series_dict['PN00-4.edf']})
df.to_csv('exercise_data/19.csv', index=False, header=False, sep='\t')