## Noise Estimation scripts

Author : Martin Deudon (2017)

Matlab scripts to estimate the noise level in micro-electrode recordings (they might also be used with Macro-electrodes)

4 different measures of the noise are available : 

* RMS (Root-Mean-Square) value of the filtered signal
* Median of the filtered signal's enveloppe : `getenvelopemedian`
* Distribution spread of the filtered signal
* Spectral peaks in the raw signal : `findspectralpeaks.m`

The main scripts is `noiseEvolution.m`

