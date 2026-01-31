function acqResults = acquisitionL1(settings)
% Function performs cold start acquisition on the collected "data". It
% searches for GPS signals of all satellites, which are listed in field
% "acqSatelliteList" in the settings structure. Function saves code phase
% and frequency of the detected signals in the "acqResults" structure.
%
% acqResults = acquisitionL1(settings)
%
%   Inputs:
%       settings      - Receiver settings. Provides information about
%                       sampling and intermediate frequencies and other
%                       parameters including the list of the satellites to
%                       be acquired.
%   Outputs:
%       acqResults    - Function saves code phases and frequencies of the
%                       detected signals in the "acqResults" structure. The
%                       field "carrFreq" is set to 0 if the signal is not
%                       detected for the given PRN number.

%--------------------------------------------------------------------------
%                             MJTL_SDR v1.0
% Copyright (C) 2026 J. Sun, Z. Tang, J. Wei and J. Zhou.
% Developed for GPS L1/L5 multi-frequency joint tracking by J. Sun, 
% Z. Tang, J. Wei and J. Zhou.
% -------------------------------------------------------------------------

%% Read GPS L1 signal
[fid, ~] = fopen(settings.fileRoute, 'rb');
fseek(fid, settings.fileType*settings.skipNumberOfBytes, 'bof'); 
% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq /(settings.codeFreqBasis / settings.codeLength));

% At least 42ms of signal are needed for fine frequency estimation
codeLen = max(42,settings.acqNonCohTime+2);
% Read data for acquisition.
data  = fread(fid, settings.fileType*codeLen*samplesPerCode, settings.dataType)';

if (settings.fileType==2)    
    data1=data(1:2:end);    
    data2=data(2:2:end);    
    longSignal=data1 + 1i .* data2;
else 
    longSignal = data;
end

%% Initialization
% Varaibles for coarse acquisition 
% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));
% Find sampling period
ts = 1 / settings.samplingFreq;
% Find phase points of 2ms local carrier wave (1ms for local duplicate,
% the other 1ms for zero padding)
phasePoints = (0 : (samplesPerCode * 2 -1)) * 2 * pi * ts;
% Number of the frequency bins for the specified search band
numberOfFreqBins = round(settings.acqSearchBand * 2 / settings.acqSearchStep) + 1;
% Carrier frequency bins to be searched
coarseFreqBin = zeros(1, numberOfFreqBins);

% Initialize acqResults
% Carrier frequencies of detected signals
acqResults.carrFreq     = zeros(1, 32);
% C/A code phases of detected signals
acqResults.codePhase    = zeros(1, 32);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 32);

% Varaibles for fine acquisition
% Carrier frequency search step for fine acquisition
fineSearchStep = 25;
% Number of the frequency bins for fine acquisition
numOfFineBins = round(settings.acqSearchStep/ fineSearchStep) + 1;
% Carrier frequencies of the fine frequency bins
fineFreqBins = zeros(1, numOfFineBins);
% Correlation values for all fine frequency bins
fineResult = zeros(1,numOfFineBins);
% Coherent integration for each of 40 codes
sumPerCode = zeros(1,40);
% Phase points of the local carrier wave
finePhasePoints = (0 : (40*samplesPerCode-1)) * 2 * pi * ts;

% Input signal power for GLRT statistic calculation
sigPower = sqrt(var(longSignal(1:samplesPerCode)) * samplesPerCode);

% Perform search for all listed PRN numbers ...
fprintf('GPS L1 C/A signal satellite acquisition result:\n');
fprintf('(');
for PRN = settings.acqSatelliteList
    %% Coarse acquisition
    % Generate C/A codes and sample them according to the sampling freq.
    caCodesTable = makeCaTable(PRN,settings);
    % Add zero-padding samples
    caCodes2ms = [caCodesTable zeros(1,samplesPerCode)];
    % Search results of all frequency bins and code shifts (for one satellite)
    results = zeros(numberOfFreqBins, samplesPerCode*2);
    % Perform DFT of C/A code
    caCodeFreqDom = conj(fft(caCodes2ms));
    
    % Make the correlation for all frequency bins
    for freqBinIndex = 1:numberOfFreqBins
        % Generate carrier wave frequency grid
        coarseFreqBin(freqBinIndex) = settings.IF + settings.acqSearchBand - ...
            settings.acqSearchStep * (freqBinIndex - 1);
        % Generate local sine and cosine
        sigCarr = exp(-1i * coarseFreqBin(freqBinIndex) * phasePoints);
        
        % Do correlation
        for nonCohIndex = 1: settings.acqNonCohTime
            % Take 2ms vectors of input data to do correlation
            signal = longSignal((nonCohIndex - 1) * samplesPerCode + ...
                1 : (nonCohIndex + 1) * samplesPerCode);
            % "Remove carrier" from the signal
            I      = real(sigCarr .* signal);
            Q      = imag(sigCarr .* signal);
            % Convert the baseband signal to frequency domain
            IQfreqDom = fft(I + 1i*Q);
            % Multiplication in the frequency domain (correlation in time domain)
            convCodeIQ = IQfreqDom .* caCodeFreqDom;
            % Perform inverse DFT and store correlation results
            cohRresult = abs(ifft(convCodeIQ));
            % Non-coherent integration
            results(freqBinIndex, :) = results(freqBinIndex, :) + cohRresult;
        end % nonCohIndex = 1: settings.acqNonCohTime
    end % frqBinIndex = 1:numberOfFreqBins
    
    %% Look for correlation peaks for coarse acquisition
    % Find the correlation peak and the carrier frequency
    [~, acqCoarseBin] = max(max(results, [], 2));
    % Find code phase of the same correlation peak
    [peakSize, codePhase] = max(max(results));
    % Store GLRT statistic
    acqResults.peakMetric(PRN) = peakSize/sigPower/settings.acqNonCohTime;
    
    %% Fine carrier frequency search
    % Do fine acquisition
    if acqResults.peakMetric(PRN) > settings.acqThreshold
        
        % Indicate PRN number of the detected signal
        fprintf('%02d ', PRN);
        
        % Prepare 20ms code, carrier and input signals
        % C/A code with 10230 chips
        caCode = generateCAcode(PRN);
        % C/A code sample index
        codeValueIndex = floor((ts * (0 : 40*samplesPerCode -1)) / ...
            (1/settings.codeFreqBasis));
        % C/A code samples
        caCode40ms = caCode(rem(codeValueIndex, settings.codeLength) + 1);
        
        % Take 40cm incoming signal for fine acquisition
        sig40cm = longSignal(codePhase:codePhase + 40*samplesPerCode -1);
        
        % Search different fine freq bins
        for fineBinIndex = 1 : numOfFineBins
            %--- Correlation for each code
            % Carrier frequencies of the frequency bins
            fineFreqBins(fineBinIndex) = coarseFreqBin(acqCoarseBin) + ...
                settings.acqSearchStep/2 - fineSearchStep * (fineBinIndex - 1);
            % Local carrier signal
            sigCarr40cm = exp(-1i * fineFreqBins(fineBinIndex) * finePhasePoints);
            % Wipe off code and carrier from incoming signals
            basebandSig = sig40cm .* caCode40ms .* sigCarr40cm;
            
            % Coherent integration for each code
            for index = 1:40
                sumPerCode(index) = sum(basebandSig( samplesPerCode * ...
                    (index - 1) + 1 : samplesPerCode * index ));
            end
            
            % Search Nav bit edge for 20 cases of Nav bit edge
            maxPower = 0;
            for comIndex = 1:20
                % Power for 20ms coherent integration
                comPower = abs(sum(sumPerCode(comIndex:comIndex+19)));
                % Maximal integration power
                maxPower = max(maxPower,comPower);
            end % Search different NH code combiniations
            fineResult(fineBinIndex) = maxPower;
        end % for numOfFineBins
        
        % Find the fine carrier freq.
        [~, maxFinBin] = max(fineResult);
        acqResults.carrFreq(PRN) = fineFreqBins(maxFinBin);
        % Save code phase acquisition result
        acqResults.codePhase(PRN) = codePhase;
        %signal found, if IF =0 just change to 1 Hz to allow processing
        if(acqResults.carrFreq(PRN) == 0)
            acqResults.carrFreq(PRN) = 1;
        end
        
    else
        % No signal with this PRN 
        fprintf('. ');
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold
    
end    % for PRN = satelliteList

%% Acquisition is over
fclose(fid);
fprintf(')\n');
