# Scripts used to reconstruct & fit $T_2^*$ maps from N-periodic SSFP data
MATLAB scripts and .dat files used to generate $T_2^*$ maps from raw data.
This is accompanies the paper 'SSFP as an alternative to multi-echo gradient echo' authored by Lally PJ, Jin Y, Huo Z, Beitone C, Chiew M, Matthews PM, Miller KL and Bangerter NK.

Zenodo DOI for this release: [![DOI](https://zenodo.org/badge/869096073.svg)](https://doi.org/10.5281/zenodo.13902686)


Also check out the interactive tutorial on N-periodic sequences: https://github.com/petelally/nperiodic_tutorial


## Notes on running the scripts
- First download the raw data from Zenodo and place the .dat files in the 'data' folder: https://doi.org/10.5281/zenodo.13900056 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13900056.svg)](https://doi.org/10.5281/zenodo.13900056)
- Then run each script in 'scripts' in turn: 1) extract individual echoes from raw data, 2) run the least squares fit, 3) apply the ROIs and compare
- **<ins>Each of the 5 images takes ~3h to fit</ins>** if using 50x multistart and a quad-core i7 CPU as done in the manuscript. Reduce the number of multistart points to speed this up or avoid it altogether.

## External dependencies
It also includes existing software from other groups (in the 'external' folder): 
- Philip Ehses' mapVBVD library (https://github.com/pehses/mapVBVD) for reading in .dat files
- some helper functions (sos, fft2c, ifft2c) from the ESPIRiT toolbox (https://people.eecs.berkeley.edu/~mlustig/Software.html)

