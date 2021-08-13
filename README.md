# OPT2k

Ocean Potential Temperature of the last 2k years (15 - 2015 CE)

These codes are associated with the paper, "The Little Ice Age and 20th Century Deep Pacific Cooling."
The inversion was identified as OPT-0015 in the paper. 

* Directory structure

`scripts`: MATLAB programs (start here) \
`src`: MATLAB source code 
`data`: input data files for MATLAB code \
`output`: output files, user-generated directory

* Data files

Challenger_WOCE_Temperature_basinwide_avg.mat (20 July 2018 version) \
Challenger_WOCE_Temperature_list.mat 

* Dependencies

OPT2k requires the TMI package v.0.8.2 or greater. You can clone it
manually with `git clone https://github.com/ggebbie/TMI` or follow the
MATLAB scripts. If TMI already exists in your configuration, the code
will not check for the latest version, but it may be required for
consistency.

Code has been tested for MATLAB version 9.10.0.1602886 (R2021a).
