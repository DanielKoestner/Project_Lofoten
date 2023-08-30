# Project_Lofoten
This repository contains both the OnArgo Matlab package for extracting and analyzing BGC-ARGO float data,
processing BBP700 and chl-a data to separate the signal associated with small and large particles, and converting
BBP700/CHLA to POC. OneArgo is extracted from the repository here: https://github.com/NOAA-PMEL/OneArgo-Mat

## Non-OneArgo Directories:

i. Scipts: Contains example code for extracting and visualizing BGC-ARGo data (Lofoten.m), 
processing BBP700/CHLA data (process_BGCARGO_BBP_CHL.m), visualizing BBP700/CHLA
associated with small/large particles (Koest23_Lofoten_Map_bs_bl.m), and converting and visualizing
BBP700/CHLA to POC (process_BGCARGO_POC.m). 

ii. matfiles: Directories containing matlab data files for use by processing and visualization 
scripts. These are empty by default and need to be populated "manually". See the readme in 
the "Scripts" directory for the structure of these files and how they can be accessed.
