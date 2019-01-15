% noiseEvolution
% Compute 4 measure for estimating the noise level in micro electrodes
% recordings when searching for action potentials. Filtered signal refers
% to the bandpass filtered signal between 300Hz and 3kHz 
% Measures are :
%   - Spectral peaks in the raw signal 
%   - Median of the filtered signal's enveloppe 
%   - Distribution spread of the filtered signal
%   - RMS value of the filtered signal
%
% See also : findspectralpeaks, getenvelopemedian
%
% Author(s) : Martin Deudon (2017)

%% Parameters
inputFilename   = '20170303-115128-001.ns5';
inputFilepath   = 'C:\Users\deudon\Desktop\SpikeSorting\_Data\CP15-Screening_20170303-115128';
resultsDirOri   = 'C:\Users\deudon\Desktop\SpikeSorting\_Results\NoiseEstimation';
%- Frame parameters
frameDuration   = 20;       % Sec
overlap         = 80;       % Percent
fMinHz          = 300;      % Hz
fMaxHz          = 3000;     % Hz
%- Spectral peaks parameters
nfft                = 2048;
minPeakProminence   = 6;
fMaxPeaks           = 8000;
%- Figure plots parameters
outFigH             = 700;

%% Create Results dir
resultsDir = createuniquedir(fullfile(resultsDirOri,inputFilename(1:end-4)));

%% Read file
fileExt = regexp(inputFilename,'\..*','match');
fileExt = fileExt{end};
switch fileExt
    case '.edf'
        EEG     = pop_biosig(fullfile(inputFilepath,inputFilename));
        %- Remove non eeg channels
        EEG     = removenoneegchannels(EEG);
        Fe      = EEG.srate;
        data    = EEG.data;
        tMax    = EEG.tmax;
        
    case '.ns5'
        NSX     = openNSx(fullfile(inputFilepath,inputFilename));
        Fe      = NSX.MetaTags.SamplingFreq;
        %- Remove last channel (trigger) if names 'ainp1'
        if strcmp(deblank(NSX.ElectrodesInfo(end).Label),'ainp1')
            data    = NSX.Data(1:end-1,:);
        else
            data    = NSX.Data;
        end
        tMax        = NSX.MetaTags.DataDurationSec;
        channames   = {NSX.ElectrodesInfo.Label};
        channames   = cellfun(@(x)deblank(x),channames,'UniformOutput',0);
end
channames = regexprep(channames,'''','p');

%%
nChan       = size(data,1);
nPnts       = size(data,2);

frameSize   = round(frameDuration*Fe);
overlapSize = (1-(overlap/100))*frameSize;
nFrames     = fix((nPnts-frameSize+overlapSize)/overlapSize);

%% Filtering 
disp('Filtering...');
[b,a]           = butter(8,(2/Fe)*[fMinHz,fMaxHz]);
dataFiltered    = zeros(nChan,nPnts);
for c=1:nChan
    dataFiltered(c,:) = filtfilt(b,a,double(data(c,:)));
    disp([num2str(c),'/',num2str(nChan)]);
end


%% Noise measurements
rmsVal          = zeros(nChan,nFrames);
envMed          = zeros(nChan,nFrames);
distribSpread   = zeros(nChan,nFrames);
peaksdB         = cell(nChan,nFrames);
fPeaksdB        = cell(nChan,nFrames);

frameInd    = 1:frameSize;

tic;
% try matlabpool open; catch; end;

for i=1:nFrames
    
    for c=1:nChan
        xFrame              = data(c,frameInd);
        xFrameFiltered      = dataFiltered(c,frameInd);
        [peaksdB{c,i},fPeaksdB{c,i}] = findspectralpeaks(xFrame,Fe,nfft,minPeakProminence,fMaxPeaks);
        rmsVal(c,i)         = rms(xFrameFiltered);        
        envMed(c,i)         = getenvelopemedian(xFrameFiltered);
        distribSpread(c,i)  = norminv(0.95, mean(xFrameFiltered), std(xFrameFiltered));        
    end   

    % Update frame index
    frameInd = int32(frameInd+overlapSize);
    % Progression bar
    disp([num2str(i),'/',num2str(nFrames)]);
    
end
% try matlabpool close; catch; end;
toc;

%% Results
tVect = linspace(0,tMax,nFrames);

% RMS
figure;
% surf(tVect,1:nChan,rmsVal,'edgecolor','none'); view([0,0,90]); axis tight;
imagesc(tVect,1:nChan,rmsVal); axis tight; axis xy;
xlabel('Time (s)'); ylabel('Channel'); colorbar;
title(['Noise Measure - RMS Value - frame duration: ',num2str(frameDuration),...
    ' s - overlap: ',num2str(overlap),'%']);
saveas(gca,fullfile(resultsDir,['rmsValue_frameDuration_',num2str(frameDuration),'s_overlap_',num2str(overlap),'.png']));

% Envelope Median
figure;
% surf(tVect,1:nChan,envMed,'edgecolor','none'); view([0,0,90]); axis tight;
imagesc(tVect,1:nChan,envMed); axis tight; axis xy;
xlabel('Time (s)'); ylabel('Channel'); colorbar;
title(['Noise Measure - Envelope Median - frame duration: ',num2str(frameDuration),...
    ' s - overlap: ',num2str(overlap),'%']);
saveas(gca,fullfile(resultsDir,['envelopeMediane_frameDuration_',num2str(frameDuration),'s_overlap_',num2str(overlap),'.png']));

% Distribution Spread
figure;
% surf(tVect,1:nChan,distribSpread,'edgecolor','none'); view([0,0,90]); axis tight;
imagesc(tVect,1:nChan,distribSpread); axis tight; axis xy;
xlabel('Time (s)'); ylabel('Channel'); colorbar;
title(['Noise Measure - Distribution Spread- frame duration: ',num2str(frameDuration),...
    ' s - overlap: ',num2str(overlap),'%']);
saveas(gca,fullfile(resultsDir,['DistributionSpread_frameDuration_',num2str(frameDuration),'s_overlap_',num2str(overlap),'.png']));

%% Results - spectrum peaks
iChan       = 12;
% Spectrum Peaks
baseVal     = -20; % in dB
fPrecision  = 30;    % Hz;
nFreqs      = 1+round(fMaxPeaks/fPrecision);
fVect       = linspace(0,fMaxPeaks,nFreqs);

peaksImAll  = baseVal*ones(nFreqs,nFrames,nChan);
for c=1:nChan
    for i=1:nFrames
        peaksdB_i   = peaksdB{c,i};
        fPeaksdB_i  = fPeaksdB{c,i};
        fIndex    	= 1+round(fPeaksdB_i/fPrecision);
        peaksImAll(fIndex,i,c) = peaksdB_i;
    end
end
peakImMean  = mean(peaksImAll,3);
peaksImage  = peaksImAll(:,:,iChan);


figure;
surf(tVect,fVect,peaksImage,'edgecolor','none'); view([0,0,90]); axis tight;
shading interp; colorbar;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title(['Noise Measure - Spectrum Peaks - Channel ',channames{iChan},' - frame duration: ',num2str(frameDuration),...
    ' s - overlap: ',num2str(overlap),'%']);
figNamePng = fullfile(resultsDir,['SpectralPeaks_channel_',channames{iChan},...
    '_frameDuration_',num2str(frameDuration),'s_overlap_',num2str(overlap),'.png']);
figNameFig = [figNamePng(1:end-3),'fig'];
pos = get(gcf,'position'); zoomFactor = outFigH/pos(4); 
set(gcf,'position',[100 100 zoomFactor*pos(3) zoomFactor*pos(4)]);
export_fig(figNamePng,'-transparent');
saveas(gca,figNameFig);


figure;
surf(tVect,fVect,peakImMean,'edgecolor','none'); view([0,0,90]); axis tight;
imagesc(tVect,fVect,peakImMean); axis tight; axis xy;
shading interp; colorbar;
xlabel('Time (s)'); ylabel('Frequency (Hz)');
title(['Noise Measure - Spectrum Peaks - Mean over all channels - frame duration: ',num2str(frameDuration),...
    ' s - overlap: ',num2str(overlap),'%']);
figNamePng = fullfile(resultsDir,['SpectralPeaks_meanChannel_',...
    '_frameDuration_',num2str(frameDuration),'s_overlap_',num2str(overlap),'.png']);
figNameFig = [figNamePng(1:end-3),'fig'];
pos = get(gcf,'position'); zoomFactor = outFigH/pos(4); 
set(gcf,'position',[100 100 zoomFactor*pos(3) zoomFactor*pos(4)]);
export_fig(figNamePng,'-transparent');
saveas(gca,figNameFig);

%% Results - spectrum peaks summed for every channel
peaksSum = zeros(nChan,nFrames);
for c=1:nChan
    for i=1:nFrames
        peaksdB_i       = peaksdB{c,i};
        peaksSum(c,i)   = sum(peaksdB_i);        
    end
end

figure;
imagesc(tVect,1:nChan,peaksSum);axis tight;
xlabel('Time (s)'); ylabel('Channel'); colorbar;
title(['Noise Measure - Spectrum peaks sum - frame duration: ',num2str(frameDuration),...
    ' s - overlap: ',num2str(overlap),'%']);
axis xy;
saveas(gca,fullfile(resultsDir,['SpectrumSpeakSum',num2str(frameDuration),'s_overlap_',num2str(overlap),'.png']));





