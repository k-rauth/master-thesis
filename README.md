# master-thesis
Code from my master thesis to calculate common envelope consequences from MESA models

SN_to_CE.py: 

Takes MESA binary models where each component ended with helium depletion (orbital period set to infinity after first SN).
Calculates effect of SN kick on orbital period.
Determines whether secondary initiates common envelope evolution and whether the binary survives it.
Also calculates needed ejection efficiency for minimum orbital period (where core radius equals Roche-lobe radius).

CE_ejection_range.py:

Uses the same code as SN_to_CE.py, modified to take a range of pre-CE orbital periods as input instead of using the value calculated from SN kick.
Calculates the needed ejection efficiency for the minimum orbital period each step and plots it.
