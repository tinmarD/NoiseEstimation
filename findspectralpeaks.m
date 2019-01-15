function [peaksdB,fPeaks] = findspectralpeaks(x, Fe, nfft, minPeakProminence, fMax, showFig)
% [peaksdB,fPeaks] = FINDSPECTRALPEAKS  (x, Fe, nfft, ...
%                                       minPeakProminence, fMax, showFig)
% Finds the peaks in the power spectral density estimated using Welch's
% method. 
%
% INPUTS : 
%   - x                 : input raw signal (1D)
%   - Fe                : sampling frequency (Hz)
%   - nfft              : number of points for the fft
%   - minPeakProminence : minimum peak prominence
%   - fMax              : frequency higher than fMax are not considered
%   - showFig           : show figure if 1
% 
% OUTPUTS : 
%   - peaksdB           : peaks gain in dB
%   - fPeaks            : peaks frequencies
%
%
% Author(s) : Martin Deudon (2017)

if nargin<6; showFig=0; end;

%% Inner Parameters
freqRangeSmooth     = 150;  % Hz

%% PSD estimation 
[pxx,f]     = pwelch(double(x),[],[],nfft,Fe, 'onesided', 'power');
idx_fmax    = find(f >= fMax, 1, 'first');
pxx         = pxx(1:idx_fmax); 
f           = f(1:idx_fmax);
pxxdB       = 10*log10(pxx);

%% PSD smoothing :
smoothWinLength = round(idx_fmax*freqRangeSmooth/fMax);
pxxdBSmooth     = movingaverage1d(pxxdB,smoothWinLength);    

%% Find peaks on the difference
pxxdBDiff   = pxxdB-pxxdBSmooth;
[~,locs]    = findpeaks(pxxdBDiff, 'MINPEAKHEIGHT', minPeakProminence);
peaksdB     = pxxdB(locs);
fPeaks      = f(locs);

%% Figure
if showFig;figure;
ax(1) = subplot(4,1,1:3);  hold on;
plot(f,pxxdB);
plot(f,pxxdBSmooth,'r','linewidth',2);
plot(f(locs),pxxdB(locs),'k^','markerfacecolor',[1 0 0]);
ylabel('Amplitude (dB)');
title('Spectral peaks');
ax(2) = subplot(414);  hold on;
plot(f,pxxdBDiff);
plot([f(1),f(end)],[minPeakProminence,minPeakProminence],'k');
plot(f(locs),pxxdBDiff(locs),'k^','markerfacecolor',[1 0 0]);
xlabel('Frequency (Hz)');
ylabel('Amplitude (dB)');
title('Difference');
linkaxes(ax,'x');
axis tight;

end;



end

