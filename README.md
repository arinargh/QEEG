# QEEG
Houses code created during my Master of Philosophy (Science) candidature at the University of Sydney

## Routine EEG Sleep Data Preparation Script
**Code:** [`all_fxns_dataprep.py`](https://github.com/arinargh/QEEG/blob/main/all_fxns_dataprep.py)

**Description**: A key responsibility during my role as a Medical Research Assistant was to prepare quantitative electroencephalography (QEEG) sleep data for collaborative use between the Healthy Brain Ageing Clinic (Brain and Mind Centre, the University of Sydney) and the Sleep and Circadian Research Group (Woolcock Institute of Medical Research).

This involved processing each overnight sleep study using internally-developed software, identifying and removing poor quality data through the visual assessment of EEG signal quality, identifying and removing statistical outliers and consolidating all data in a standardised format. I collated data from >200 sleep studies in total. 

I created this program to automate these processes, saving time and reducing human error rates.

## Slow Oscillation-Sleep Spindle Coupling
My empirical project utilised overnight sleep data collected through 256-channel high-density EEG systems to

1. explore the synchronised activity between slow oscillations and sleep spindles, two distinct features of non-rapid eye movement (NREM) sleep, and 
2. assess the associations between slow oscillation-sleep spindle coupling strength and overnight memory consolidation processes

in older adults presenting with versus without a risk of dementia development.

### Coupling Frequency: Processing Pipeline
**Code (Module):** [`all_fxns_aac.py`](https://github.com/arinargh/QEEG/blob/main/all_fxns_aac.py)
 | **(Script):** [`run_aac.py`](https://github.com/arinargh/QEEG/blob/main/run_aac.py)

**Description:** This pipeline was built to assess slow oscillation-sleep spindle coupling through the parameter of coupling frequency (_ _operationalised as the ratio of general co-occurrence of sleep spindles to slow oscillations_ _)

The criteria for co-occurrence follows that of established methodology; a sleep spindle is defined as co-occurring (i.e. coupled) to a slow oscillation if at least half of the sleep spindle's activity fell within a time window of |1.2| seconds around the trough of that slow oscillation. Example:

![example_coupling_freq](https://user-images.githubusercontent.com/129021696/231706307-6a9836a9-a4c0-40f7-a8c5-29ff1f2b2ffb.jpg)


### Coupling Precision: Data Preparation Script
**Code (Module):** [`all_fxns_pac.py`](https://github.com/arinargh/QEEG/blob/main/all_fxns_pac.py)
 | **(Script):** [`run_pac.py`](https://github.com/arinargh/QEEG/blob/main/run_pac.py)

**Description:** An existing pipeline was built by the research team from Concordia University to assess slow oscillation-sleep spindle coupling through the parameter of coupling precision (_ _operationalised as the preferred phase angle of co-occurring sleep spindles to slow oscillations_ _)

The pipeline was used to analyse data recorded through routine EEG (which comprise up to 25 channels), but was subsequently extended to allow for the analysis in high-density EEG. 

I wrote this script to automate the organisation and cleaning of the large volumes of output data to quickly prepare analysis-ready datasets.
