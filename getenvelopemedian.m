function [medEnv] = getenvelopemedian (x)
%[medEnv] = GETENVELOPEMEDIAN (x)
%   Calculate the envelope of the signal x using Hilbert's transform. Then
%   calculate the median of the envelope giving a measure of the width of
%   the signal. 
%
% INPUTS : 
%   - x             : input signal
%
% OUTPUTS :
%   - medEnv        : median of the envelope
%
%
% Author(s) : Martin Deudon (2017)
 

%- Calcul the mean of the signal
xMean   = mean(x);
%- Compute the hilbert transform of the signal minus the mean
env     = abs(hilbert(x-xMean));
%- Add the mean to the envelope
yUp     = env+xMean;
yLow    = xMean-env;
%- Median of the envelope gives an estimate of the width of the signal
medEnv  = median(abs([yUp,yLow]));



end

