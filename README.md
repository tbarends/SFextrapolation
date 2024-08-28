# SFextrapolation
EXTRAPOLATION OF STRUCTURE FACTORS

THIS SOFTWARE IS DISTRIBUTED UNDER THE GNU GPLv3 LICENSE

1. SYSTEM REQUIREMENTS
-Linux (tested on 3.10.0-1160.119.1.e17.x86-64 
-CCP4 (tested for v8.0.017
-PHENIX (tested for 1.19.2-4158)
-python2.7 with numpy, matplotlib, scipy 

2. INSTALLATION INSTRUCTIONS
Copy all files to a stand-alone directory

3. and 4. DEMO ON ACTUAL DATA
The directory includes the original data associated with the paper. 
Run the extrapolation script by entering:

./extrapolate_auto

the script will first calculate Flight-Fdark map coefficients. As provided,
these will be Q-weighted. The coefficients will be in

Qweighted_difference_occ0.25_dm0.6.mtz    

There are plots of the delta-F amplitudes vs. resolution as well as of the 
Q-weighting factors vs. resolution in

Diffs_occ0.25_dm0.6.png 
and
Qs_occ0.25_dm0.6.png 
, respectively

The script has been set up to then run extrapolations with assumed occupancies 
between 0.1 and 0.5, in steps of 0.05. The resulting Q-weighted extrapolated Fs are
in the files  

Qweighted_extrapolated_occ0.25_dm0.6_#.###.mtz 
with #.### being the occupancies used. There is also a file called 

Qweighted_extrapolated_occ0.25_dm0.6_0.250autoocc.mtz
                                 
which contains the coefficients for a preset occupancy of 0.25. These can be
read by e.g. COOT for map inspection or be used for refinement.

THESE ARE THE FILES USED FOR THE ANALYSIS IN THE ASSOCIATED PAPER.
The other files are made by experimental features of the script 
such as automatic occupancy determination, density modification
etc. and should not be used.

To change settings, edit the "CHANGE PARAMETERS HERE" section at the top of
the script. Setting qweight to FALSE wil switch off the Q-weighting, 
changing the jobname parameter will change the filename root for the output
files. preset_occ controls the final occupancy used for extrapolation.
