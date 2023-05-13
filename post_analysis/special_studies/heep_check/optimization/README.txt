Brief:
This directory contains hydrogen elastic optimization studies for
CaFe using the singles data (SHMS: 8.55 GeV/c,  6.8, 7.495, 8.295 deg)
taken on Aug 08 2022 and the coincidence run (SHMS: 8.55, 8.3 deg) taken
on Sep 2022. As of May 12, Holly optimized the SHMS delta and removed the
W vs. yptar dependence which was the cause of the W not being aligned
between data and SIMC for the 8.3 deg setting. On directory ./step1,
this W offset will be checked without momentum corrections to make
sure the offset is the same for single and coincidence before proceding
to next step


Directory Structure

------------------------
./current_status_may12:
------------------------
slides with data/simc comparison of cafe h(e,e') singles
and h(e,e'p) coincidence using the newly optimized SHMS delta
matrix by Holly on May 11. This is mainly to check where we stand
in terms of the kinematics (W, Pm, Em, Pmx,y,z) alignment between data
and SIMC, and determine the course of action to take for improvements.


--------
./step1
--------
removed SHMS central momentum correction factor to check/compare W offsets
between the hydrogen elastic singles and coincidence runs
Ideally, for the particular case of the 8.3 deg setting for both singles
and coincidence, the W should be offset by the same amount as the
SHMS (e) spectrometer was set at the same angle/momentum
