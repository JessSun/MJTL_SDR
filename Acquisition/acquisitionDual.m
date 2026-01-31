function [acqResultsL1,acqResultsL5] = acquisitionDual(settingsL1,settingsL5)
% Function performs cold start acquisition on the collected "data". It
% searches for GPS signals of all satellites, which are listed in field
% "acqSatelliteList" in the settings structure. Function saves code phase
% and frequency of the detected signals in the "acqResults" structure.
%
% acqResults = acquisitionL5(settings)
%
%   Inputs:
%       settingsL1    - L1 signal paremeters settings. Provides information
%                       about sampling and intermediate frequencies and 
%                       otherparameters including the list of the satellites
%                       to be acquired.
%       settingsL5    - L5 signal paremeters settings. Provides information
%                       about sampling and intermediate frequencies and 
%                       otherparameters including the list of the satellites
%                       to be acquired.
%   Outputs:
%       acqResultsL1  - Function saves code phases and frequencies of the
%                       detected L1 signals in the "acqResults" structure. 
%                       The field "carrFreq" is set to 0 if the signal is 
%                       not detected for the given PRN number.
%       acqResultsL5  - Function saves code phases and frequencies of the
%                       detected L5 signals in the "acqResults" structure. 
%                       The field "carrFreq" is set to 0 if the signal is 
%                       not detected for the given PRN number.

%--------------------------------------------------------------------------
%                             MJTL_SDR v1.0
% Copyright (C) 2026 J. Sun, Z. Tang, J. Wei and J. Zhou.
% Developed for GPS L1/L5 multi-frequency joint tracking by J. Sun, 
% Z. Tang, J. Wei and J. Zhou.
% -------------------------------------------------------------------------

%% Read GPS L1 signal
% Move the starting point of processing. 
[fidL1, ~] = fopen(settingsL1.fileRoute, 'rb');
fseek(fidL1, settingsL1.fileType*settingsL1.skipNumberOfBytes, 'bof'); 
% Find number of samples per spreading code
samplesPerCode = round(settingsL1.samplingFreq /(settingsL1.codeFreqBasis / settingsL1.codeLength));

% At least 42ms of signal are needed for fine frequency estimation
codeLenL1 = max(42,settingsL1.acqNonCohTime+2);
% Read data for acquisition.
dataL1  = fread(fidL1, settingsL1.fileType*codeLenL1*samplesPerCode, settingsL1.dataType)';

if (settingsL1.fileType==2)    
    data1=dataL1(1:2:end);    
    data2=dataL1(2:2:end);    
    longSignalL1=data1 + 1i .* data2;
end

%% Read GPS L5 signal
% Move the starting point of processing. 
[fidL5, ~] = fopen(settingsL5.fileRoute, 'rb');
fseek(fidL5, settingsL5.fileType*settingsL5.skipNumberOfBytes, 'bof'); 
% Find number of samples per spreading code
samplesPerCodeL5 = round(settingsL5.samplingFreq / (settingsL5.codeFreqBasis / settingsL5.codeLength));

% At least 42ms of signal are needed for fine frequency estimation
codeLenL5 = max(42,settingsL5.acqNonCohTime+2);
% Read data for acquisition.
dataL5  = fread(fidL5, settingsL5.fileType*codeLenL5*samplesPerCodeL5, settingsL5.dataType)';

if (settingsL5.fileType==2)    
    data1=dataL5(1:2:end);    
    data2=dataL5(2:2:end);    
    longSignalL5=data1 + 1i .* data2;
end
sigPowerL5 = sqrt(var(longSignalL5(1:samplesPerCodeL5)) * samplesPerCodeL5);


%% Initialization ===================================================
% Varaibles for coarse acquisition
% Find number of samples per spreading code
samplesPerCode = round(settingsL1.samplingFreq / ...
    (settingsL1.codeFreqBasis / settingsL1.codeLength));
% Find sampling period
ts = 1 / settingsL1.samplingFreq;
% Find phase points of 2ms local carrier wave (1ms for local duplicate,the other 1ms for zero padding)
phasePoints = (0 : (samplesPerCode * 2 -1)) * 2 * pi * ts;
% Number of the frequency bins for the specified search band
numberOfFreqBins = round(settingsL1.acqSearchBand * 2 / settingsL1.acqSearchStep) + 1;
% Carrier frequency bins to be searched
coarseFreqBin = zeros(1, numberOfFreqBins);

%% Initialize acqResults
% Carrier frequencies of detected signals
acqResultsL1.carrFreq     = zeros(1, 32);
% C/A code phases of detected signals
acqResultsL1.codePhase    = zeros(1, 32);
% Correlation peak ratios of the detected signals
acqResultsL1.peakMetric   = zeros(1, 32);
% Carrier frequencies of detected signals
acqResultsL5.carrFreq     = zeros(1, 32);
% C/A code phases of detected signals
acqResultsL5.codePhase    = zeros(1, 32);
% Correlation peak ratios of the detected signals
acqResultsL5.peakMetric   = zeros(1, 32);

% Varaibles for fine acquisition
% Carrier frequency search step for fine acquisition
fineSearchStep = 25;
% Number of the frequency bins for fine acquisition
numOfFineBins = round(settingsL1.acqSearchStep/ fineSearchStep) + 1;
% Carrier frequencies of the fine frequency bins
fineFreqBins = zeros(1, numOfFineBins);
% Correlation values for all fine frequency bins
fineResult = zeros(1,numOfFineBins);
% Coherent integration for each of 40 codes
sumPerCode = zeros(1,40);
% Phase points of the local carrier wave
finePhasePoints = (0 : (40*samplesPerCode-1)) * 2 * pi * ts;

% Input signal power for GLRT statistic calculation
sigPower = sqrt(var(longSignalL1(1:samplesPerCode)) * samplesPerCode);

% Perform search for all listed PRN numbers ...
fprintf('GPS L1 C/A signal satellite acquisition result:\n');
fprintf('(');
for PRN = settingsL1.acqSatelliteList
    %% Coarse acquisition
    % Generate C/A codes and sample them according to the sampling freq.
    caCodesTable = makeCaTable(PRN,settingsL1);
    % Add zero-padding samples
    caCodes2ms = [caCodesTable zeros(1,samplesPerCode)];
    % Search results of all frequency bins and code shifts (for one satellite)
    results = zeros(numberOfFreqBins, samplesPerCode*2);
    % Perform DFT of C/A code
    caCodeFreqDom = conj(fft(caCodes2ms));
    
    %--- Make the correlation for all frequency bins
    for freqBinIndex = 1:numberOfFreqBins
        % Generate carrier wave frequency grid
        coarseFreqBin(freqBinIndex) = settingsL1.IF + settingsL1.acqSearchBand - ...
            settingsL1.acqSearchStep * (freqBinIndex - 1);
        % Generate local sine and cosine
        sigCarr = exp(-1i * coarseFreqBin(freqBinIndex) * phasePoints);
        
        %--- Do correlation -----------------------------------------------
        for nonCohIndex = 1: settingsL1.acqNonCohTime
            % Take 2ms vectors of input data to do correlation
            signal = longSignalL1((nonCohIndex - 1) * samplesPerCode + ...
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
    
    %% Look for correlation peaks for coarse acquisition ============
    % Find the correlation peak and the carrier frequency
    [~, acqCoarseBin] = max(max(results, [], 2));
    % Find code phase of the same correlation peak
    [peakSize, codePhase] = max(max(results));
    % Store GLRT statistic
    acqResultsL1.peakMetric(PRN) = peakSize/sigPower/settingsL1.acqNonCohTime;
    
    %% Fine carrier frequency search ================================
    % Do fine acquisition
    if acqResultsL1.peakMetric(PRN) > settingsL1.acqThreshold
        
        % Indicate PRN number of the detected signal
        fprintf('%02d ', PRN);
        
        % Prepare 20ms code, carrier and input signals
        % C/A code with 10230 chips
        caCode = generateCAcode(PRN);
        % C/A code sample index
        codeValueIndex = floor((ts * (0 : 40*samplesPerCode -1)) / ...
            (1/settingsL1.codeFreqBasis));
        % C/A code samples
        caCode40ms = caCode(rem(codeValueIndex, settingsL1.codeLength) + 1);
        
        % Take 40cm incoming signal for fine acquisition
        sig40cm = longSignalL1(codePhase:codePhase + 40*samplesPerCode -1);
        
        %% First get the acquisition result of GPS L1 frequency
        % Search different fine freq bins
        for fineBinIndex = 1 : numOfFineBins
            % Correlation for each code
            % Carrier frequencies of the frequency bins
            fineFreqBins(fineBinIndex) = coarseFreqBin(acqCoarseBin) + ...
                settingsL1.acqSearchStep/2 - fineSearchStep * (fineBinIndex - 1);
            % Local carrier signal
            sigCarr40cm = exp(-1i * fineFreqBins(fineBinIndex) * finePhasePoints);
            % Wipe off code and carrier from incoming signals
            basebandSig = sig40cm .* caCode40ms .* sigCarr40cm;
            
            % Coherent integration for each code
            for index = 1:40
                sumPerCode(index) = sum(basebandSig( samplesPerCode * ...
                    (index - 1) + 1 : samplesPerCode * index ));
            end
            
            %  Search Nav bit edge for
            % 20 cases of Nav bit edge
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
        acqResultsL1.carrFreq(PRN) = fineFreqBins(maxFinBin);
        % Save code phase acquisition result
        acqResultsL1.codePhase(PRN) = codePhase;
        % signal found, if IF =0 just change to 1 Hz to allow processing
        if(acqResultsL1.carrFreq(PRN) == 0)
            acqResultsL1.carrFreq(PRN) = 1;
        end

        %% Second calculate the acquisition result of GPS L5 frequency
        doppler = acqResultsL1.carrFreq(PRN) - settingsL1.IF; 
        acqResultsL5.carrFreq(PRN) = settingsL5.IF + doppler*settingsL5.carrFreqBasis/settingsL1.carrFreqBasis;

        L5ICodesTable = makeL5ITable(PRN,settingsL5);
        L5QCodesTable = makeL5QTable(PRN,settingsL5);    

        % generate local code duplicate to do correlate
        localL5ICode = [L5ICodesTable, zeros(1,samplesPerCodeL5)];
        localL5QCode = [L5QCodesTable, zeros(1,samplesPerCodeL5)];
        L5ICodeFreqDom = conj(fft(localL5ICode));
        L5QCodeFreqDom = conj(fft(localL5QCode));
        resultsL5 = zeros(1, samplesPerCodeL5*2);

        % Generate local sine and cosine
        sigCarrL5 = exp(-1i * acqResultsL5.carrFreq(PRN) * phasePoints);

        % Do correlation
        for nonCohIndex = 1: settingsL5.acqNonCohTime
            % Take 2ms vectors of input data to do correlation
            signalL5 = longSignalL5((nonCohIndex - 1) * samplesPerCodeL5 + 1 : (nonCohIndex + 1) * samplesPerCodeL5);
            % "Remove carrier" from the signal
            I      = real(sigCarrL5 .* signalL5);
            Q      = imag(sigCarrL5 .* signalL5);

            % Convert the baseband signal to frequency domain
            IQfreqDom = fft(I + 1i*Q);

            % Multiplication in the frequency domain (correlation in time domain)
            convL5I = IQfreqDom .* L5ICodeFreqDom;
            convL5Q = IQfreqDom .* L5QCodeFreqDom;

            % Perform inverse DFT and store correlation results
            cohRresult = abs(ifft(convL5I)) + abs(ifft(convL5Q));
            % Non-coherent integration
            resultsL5 = resultsL5 + cohRresult;            
        end % nonCohIndex = 1: settings.acqNonCohTime

        [peakSize, codePhase] = max(resultsL5);
        peakMetric = peakSize/sigPowerL5/settingsL5.acqNonCohTime;
        if peakMetric>settingsL5.acqThreshold
            acqResultsL5.codePhase(PRN) = codePhase;
            acqResultsL5.peakMetric(PRN) = peakMetric;
        else
            acqResultsL5.carrFreq(PRN) = 0;
        end
        
    else
        % No signal with this PRN 
        fprintf('. ');
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold
    
end    % for PRN = satelliteList
fprintf(')\n');
fprintf('GPS L5 signal satellite acquisition result:\n');

%% Output L5 signal satellite acquisition result
fprintf('(');
for PRN = settingsL5.acqSatelliteList
    if acqResultsL5.carrFreq(PRN)
        fprintf('%02d ', PRN);
    else
        fprintf('. ');
    end
end
fprintf(')\n');

%% Acquisition is over
fclose(fidL1);

