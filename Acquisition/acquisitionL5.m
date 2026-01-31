function acqResults = acquisitionL5(settings)
% Function performs cold start acquisition on the collected "data". It
% searches for GPS signals of all satellites, which are listed in field
% "acqSatelliteList" in the settings structure. Function saves code phase
% and frequency of the detected signals in the "acqResults" structure.
%
% acqResults = acquisitionL5(settings)
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

%% Read GPS L5 signal
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
end

%% Acquisition initialization
%  Varaibles for coarse acquisition
% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / (settings.codeFreqBasis / settings.codeLength));
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
% PRN code phases of detected signals
acqResults.codePhase    = zeros(1, 32);
% Correlation peak ratios of the detected signals
acqResults.peakMetric   = zeros(1, 32);

% Varaibles for fine acquisition
% Antipodal form of Neuman-Hoffman (NH) codes
NHcode = [1 1 1 1 1 -1 1 1 -1 -1 1 -1 1 -1 1 1 -1 -1 -1 1];
% Carrier frequency search step for fine acquisition
fineSearchStep = 25;
% Number of the frequency bins for fine acquisition
NumOfFineBins = round(settings.acqSearchStep / fineSearchStep) + 1;
% Carrier frequencies of the frequency bins
FineFreqBins     = zeros(1, NumOfFineBins);
% Search results of all frequency bins
FineResult = zeros(1,NumOfFineBins);
% Coherent integration for each code
sumPerCode = zeros(1,20);
% Find phase points of the local carrier wave
finePhasePoints = (0 : (20*samplesPerCode-1)) * 2 * pi * ts;
% Input signal power for GLRT statistic calculation
sigPower = sqrt(var(longSignal(1:samplesPerCode)) * samplesPerCode);

% Perform search for all listed PRN numbers ...
fprintf('GPS L5 signal satellite acquisition result:\n');
fprintf('(');
for PRN = settings.acqSatelliteList

    %% Coarse acquisition ===========================================
    % Generate all L5I + L5Q code and sample them according to the sampling freq.
    L5ICodesTable = makeL5ITable(PRN,settings);
    L5QCodesTable = makeL5QTable(PRN,settings);

    % generate local code duplicate to do correlate
    localL5ICode = [L5ICodesTable, zeros(1,samplesPerCode)];
    localL5QCode = [L5QCodesTable, zeros(1,samplesPerCode)];
    % Search results of all frequency bins and code shifts (for one satellite)
    results = zeros(numberOfFreqBins, samplesPerCode*2);

    % Perform DFT of PRN code
    L5ICodeFreqDom = conj(fft(localL5ICode));
    L5QCodeFreqDom = conj(fft(localL5QCode));

    % Make the correlation for whole frequency band (for all freq. bins)
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
            convL5I = IQfreqDom .* L5ICodeFreqDom;
            convL5Q = IQfreqDom .* L5QCodeFreqDom;

            % Perform inverse DFT and store correlation results
            cohRresult = abs(ifft(convL5I)) + abs(ifft(convL5Q));
            % Non-coherent integration
            results(freqBinIndex, :) = results(freqBinIndex, :) + cohRresult;
        end % nonCohIndex = 1: settings.acqNonCohTime
    end % freqBinIndex = 1:numberOfFreqBins

    %% Look for correlation peaks for coarse acquisition
    % Find the correlation peak and the carrier frequency
    [~, acqCoarseBin] = max(max(results, [], 2));
    % Find code phase of the same correlation peak
    [peakSize, codePhase] = max(max(results));
    % Store GLRT statistic
    acqResults.peakMetric(PRN) = peakSize/sigPower/settings.acqNonCohTime;

    %% Fine resolution frequency search
    % If the result is above threshold, then there is a signal ...
    if acqResults.peakMetric(PRN) > settings.acqThreshold
        % Indicate PRN number of the detected signal
        fprintf('%02d ', PRN);

        % Generate 20msec long L5Q codes sequence for given PRN
        L5QCode = generateL5Qcode(PRN,settings);
        codeValueIndex = floor((ts * (1:20*samplesPerCode)) / (1/settings.codeFreqBasis));
        longL5QCode = L5QCode((rem(codeValueIndex, settings.codeLength) + 1));

        % 10cm incoming signal
        sig20cm = longSignal(codePhase:codePhase + 20*samplesPerCode -1);

        % Search different frequency bins
        for FineBinIndex = 1 : NumOfFineBins

            % Carrier frequencies of the frequency bins
            FineFreqBins(FineBinIndex) = coarseFreqBin(acqCoarseBin) + settings.acqSearchStep/2 - fineSearchStep * (FineBinIndex - 1);
            % Generate local sine and cosine
            sigCarr20cm = exp(-1i*FineFreqBins(FineBinIndex) * finePhasePoints);

            % Wipe off L5I code and carrier from incoming signals to
            % produce baseband signal
            basebandSig = longL5QCode .* sigCarr20cm .* sig20cm;

            % Coherent integration for each code
            for index = 1:20
                sumPerCode(index) = sum( basebandSig( samplesPerCode*(index-1)+1:samplesPerCode*index ) );
            end

            % Initialize maximal power for 10 NH code combiniations
            maxPower = 0;
            %--- Search different NH code combinations -------------------
            for comIndex = 1:20
                % Wipe off NH code
                sumPerCodeNH = sumPerCode .* NHcode;

                % Maximal coherent power for different NH code combiniations
                maxPower = max(maxPower,abs(sum(sumPerCodeNH)));

                % Shift NH code for next  combiniation
                NHcode = circshift(NHcode',1)';
            end % Search different NH code combiniations

            FineResult(FineBinIndex) = maxPower;

        end % FineBinIndex = 1 : NumOfFineBins

        % Find the fine carrier freq. 
        % Corresponding to the largest noncoherent power
        [~,maxFinBin] = max(FineResult);
        acqResults.carrFreq(PRN) = FineFreqBins(maxFinBin);

        % Code phase acquisition result
        acqResults.codePhase(PRN) = codePhase;

        %signal found, if IF = 0 just change to 1 Hz to allow processing
        if(acqResults.carrFreq(PRN) == 0)
            acqResults.carrFreq(PRN) = 1;
        end

    else
        % No signal with this PRN
        fprintf('. ');
    end   % if (peakSize/secondPeakSize) > settings.acqThreshold

end    % for PRN = satelliteList

%% Acquisition is over
fprintf(')\n');
