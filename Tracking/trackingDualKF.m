function [trackResultsL1,trackResultsL5]= trackingDualKF(channel, settings)
% Kalman-based joint tracking function for GPS L1/L5 signals.
%
% [trackResultsL1,trackResultsL5]= trackingDualKF(channel, settings)
%
%   Inputs:
%       channel         - Acquisition results from GPS L1 and L5 signals
%       settings        - Receiver settings.
%   Outputs:
%       trackResultsL1  - Tracking results for GPS L1 signals
%       trackResultsL5  - Tracking results for GPS L5 signals

%--------------------------------------------------------------------------
%                             MJTL_SDR v1.0
% Copyright (C) 2026 J. Sun, Z. Tang, J. Wei and J. Zhou.
% Developed for GPS L1/L5 multi-frequency joint tracking by J. Sun, 
% Z. Tang, J. Wei and J. Zhou.
% -------------------------------------------------------------------------

%% Initialize result structure ============================================
[fidL1, ~] = fopen(settings.fileRoute1, 'rb');
[fidL5, ~] = fopen(settings.fileRoute5, 'rb');

%% Intialisation of L1 tracking variables =================================
trackResultsL1.status         = '-';      % No tracked signal, or lost lock
trackResultsL1.absoluteSample = zeros(1, settings.msToProcess);
trackResultsL1.codeFreq       = inf(1, settings.msToProcess);
trackResultsL1.carrFreq       = inf(1, settings.msToProcess);
trackResultsL1.I_P            = zeros(1, settings.msToProcess);
trackResultsL1.I_E            = zeros(1, settings.msToProcess);
trackResultsL1.I_L            = zeros(1, settings.msToProcess);
trackResultsL1.Q_E            = zeros(1, settings.msToProcess);
trackResultsL1.Q_P            = zeros(1, settings.msToProcess);
trackResultsL1.Q_L            = zeros(1, settings.msToProcess);
% Loop discriminators
trackResultsL1.dllDiscr       = inf(1, settings.msToProcess);
trackResultsL1.dllDiscrFilt   = inf(1, settings.msToProcess);
trackResultsL1.pllDiscr       = inf(1, settings.msToProcess);
trackResultsL1.pllDiscrFilt   = inf(1, settings.msToProcess);
% Remain code and carrier phase
trackResultsL1.remCodePhase       = inf(1, settings.msToProcess);
trackResultsL1.remCarrPhase       = inf(1, settings.msToProcess);
%C/No
trackResultsL1.CNo.VSMValue = zeros(1,floor(settings.msToProcess/settings.CNo.VSMinterval));
trackResultsL1.CNo.VSMIndex = zeros(1,floor(settings.msToProcess/settings.CNo.VSMinterval));
% Copy initial settings for all channels
trackResultsL1 = repmat(trackResultsL1, 1, length(channel));

%% Intialisation of L5 tracking variables =================================
trackResultsL5.status         = '-';      % No tracked signal, or lost lock
trackResultsL5.absoluteSample = zeros(1, settings.msToProcess);
trackResultsL5.codeFreq       = inf(1, settings.msToProcess);
trackResultsL5.carrFreq       = inf(1, settings.msToProcess);
trackResultsL5.I_P            = zeros(1, settings.msToProcess);
trackResultsL5.I_E            = zeros(1, settings.msToProcess);
trackResultsL5.I_L            = zeros(1, settings.msToProcess);
trackResultsL5.Q_E            = zeros(1, settings.msToProcess);
trackResultsL5.Q_P            = zeros(1, settings.msToProcess);
trackResultsL5.Q_L            = zeros(1, settings.msToProcess);
% for pilot signal
if (settings.pilotTRKflag == 1)
    trackResultsL5.Pilot_I_P  = zeros(1, settings.msToProcess);
    trackResultsL5.Pilot_Q_P  = zeros(1, settings.msToProcess);
end
% Loop discriminators
trackResultsL5.dllDiscr       = inf(1, settings.msToProcess);
trackResultsL5.dllDiscrFilt   = inf(1, settings.msToProcess);
trackResultsL5.pllDiscr       = inf(1, settings.msToProcess);
trackResultsL5.pllDiscrFilt   = inf(1, settings.msToProcess);
% Remain code and carrier phase
trackResultsL5.remCodePhase   = inf(1, settings.msToProcess);
trackResultsL5.remCarrPhase   = inf(1, settings.msToProcess);
% C/No
trackResultsL5.CNo.VSMValue = zeros(1,floor(settings.msToProcess/settings.CNo.VSMinterval));
trackResultsL5.CNo.VSMIndex = zeros(1,floor(settings.msToProcess/settings.CNo.VSMinterval));
% Copy initial settings for all channels 
trackResultsL5 = repmat(trackResultsL5, 1, length(channel));

% Processing time
codePeriods = settings.msToProcess;     % For GPS one C/A code is one ms
% DLL variables
earlyLateSpcL1 = settings.dllCorrelatorSpacing1;
earlyLateSpcL5 = settings.dllCorrelatorSpacing5;

% Number of acqusired signals
TrackedNr =0 ;
for channelNr = 1:length(channel)
    if ((channel(channelNr).statusL1 == 'T') && (channel(channelNr).statusL5 == 'T'))
        TrackedNr = TrackedNr+1;
    end
end

% Start waitbar
hwb = waitbar(0,'Tracking...');

%Adjust the size of the waitbar to insert text
CNoPos=get(hwb,'Position');
set(hwb,'Position',[CNoPos(1),CNoPos(2),CNoPos(3),90],'Visible','on');

if (settings.fileType==1)
    dataAdaptCoeff=1;
else
    dataAdaptCoeff=2;
end


%% Kalman Filter settting =================================================
% Dimension of the KF state vector
num_state   = 6; 
% Frequency scaling factor
beta1 = settings.codeFreqBasis1/settings.carrFreqBasis5;
beta2 = settings.codeFreqBasis5/settings.carrFreqBasis5;
alpha = settings.carrFreqBasis1/settings.carrFreqBasis5;
% Transition Matrix
Transistion_Matrix = [1,0,0,0,beta1*settings.intTime,0;...
                      0,1,0,0,beta2*settings.intTime,0;...
                      0,0,1,0,alpha*settings.intTime,alpha*settings.intTime*settings.intTime/2;...
                      0,0,0,1,settings.intTime,settings.intTime*settings.intTime/2;...
                      0,0,0,0,1,settings.intTime;...
                      0,0,0,0,0,1];
% Measure Matrix
Measure_Matrix  = [1,0,0,0,beta1*settings.intTime/2,0;...
                   0,0,1,0,alpha*settings.intTime/2,alpha*settings.intTime*settings.intTime/6;...
                   0,1,0,0,beta2*settings.intTime/2,0;...
                   0,0,0,1,settings.intTime/2,settings.intTime*settings.intTime/6;];
% Process noise covariance matrix
process_noise = diag([10^-9,10^-9,10^-6,10^-6,0.0001,0.01]);
% Initial estimate error covariance matrix
state_cov = diag([10,10,10,10,1000,10]);
% Measurement noise covariance matrix
CNoValueL1 = 35;
CNoValueL5 = 35;
tempCN0L1 = 10^(CNoValueL1/10);
tau_varL1 = (earlyLateSpcL1/(tempCN0L1*4*settings.intTime))*(1+2/((2-earlyLateSpcL1)*tempCN0L1*settings.intTime));
phi_varL1 = (1/(tempCN0L1*2*settings.intTime))*(1+1/(2*tempCN0L1*settings.intTime));
tempCN0L5 = 10^(CNoValueL5/10);
tau_varL5 = (earlyLateSpcL5/(tempCN0L5*4*settings.intTime))*(1+2/((2-earlyLateSpcL5)*tempCN0L5*settings.intTime));
phi_varL5 = (1/(tempCN0L5*2*settings.intTime))*(1+1/(2*tempCN0L5*settings.intTime));
mesurement_noise = diag([tau_varL1,phi_varL1,tau_varL5,phi_varL5]);

%% Start processing channels ==============================================
for channelNr = 1:length(channel)

    if (channel(channelNr).PRN ~= 0)
        % Move the starting point of processing. 
        fseek(fidL1, dataAdaptCoeff*(settings.skipNumberOfBytes + channel(channelNr).codePhaseL1-1),'bof');
        fseek(fidL5, dataAdaptCoeff*(settings.skipNumberOfBytes + channel(channelNr).codePhaseL5-1),'bof');

        %% L1 code generation
        % Get a vector with the C/A code sampled 1x/chip
        caCode = generateCAcode(channel(channelNr).PRN);
        % Then make it possible to do early and late versions
        caCode = [caCode(1023) caCode caCode(1)];

        %% L5 code generation
        % Get a vector with the C/A code sampled 1x/chip
        L5ICode = generateL5Icode(channel(channelNr).PRN,settings.st5);
        % Then make it possible to do early and late versions
        L5ICode = [L5ICode(settings.codeLength5) L5ICode L5ICode(1)];    
        if (settings.pilotTRKflag == 1)
            % Get a vector with the C/A code sampled 1x/chip
            L5QCode = generateL5Qcode(channel(channelNr).PRN,settings);
            % Then make it possible to do early and late versions
            L5QCode = [L5QCode(settings.codeLength) L5QCode L5QCode(1)]; 
        end

        %% --- Perform various initializations for L1 and L5 Channels =====
        % L1 acquired initialization
        % define initial code frequency basis of NCO
        codeFreqL1      = settings.codeFreqBasis1;
        codeNcoL1       = channel(channelNr).acquiredFreqL1*settings.codeFreqBasis1/settings.carrFreqBasis1;
        % Define residual code phase (in chips)
        remCodePhaseL1  = 0.0;
        % Define carrier frequency which is used over whole tracking period
        carrFreqL1      = channel(channelNr).acquiredFreqL1;
        carrFreqBasisL1 = 0.0;      
        carrNcoL1       = channel(channelNr).acquiredFreqL1;
        % Define residual carrier phase
        remCarrPhaseL1  = 0.0;
        % For C/No computation
        vsmCntL1  = 0;CNoL1 = 0;
        % L5 acquired initialization
        % define initial code frequency basis of NCO
        codeFreqL5      = settings.codeFreqBasis5;
        codeNcoL5       = channel(channelNr).acquiredFreqL5*settings.codeFreqBasis5/settings.carrFreqBasis5;
        % Define residual code phase (in chips)
        remCodePhaseL5  = 0.0;
        % Define carrier frequency which is used over whole tracking period
        carrFreqL5      = channel(channelNr).acquiredFreqL5;
        carrFreqBasisL5 = 0.0;      
        carrNcoL5       = channel(channelNr).acquiredFreqL5;
        % Define residual carrier phase
        remCarrPhaseL5  = 0.0;
        % For C/No computation
        vsmCntL5  = 0;CNoL5 = 0;
        % Kalman filter state vector initialization 
        error_state = zeros(num_state,1);
        error_state(5) = channel(channelNr).acquiredFreqL5;

        % Process the number of specified code periods
        for loopCnt =  1:codePeriods

            %% GUI update -------------------------------------------------------------
            % The GUI is updated every 200ms. This way Matlab GUI is still
            % responsive enough. At the same time Matlab is not occupied
            % all the time with GUI task.
            if (rem(loopCnt, 200) == 0)

                Ln = newline;
                trackingStatus=['Tracking: Ch ', int2str(channelNr), ...
                    ' of ', int2str(TrackedNr),Ln ...
                    'PRN: ', int2str(channel(channelNr).PRN),Ln ...
                    'Completed ',int2str(loopCnt), ...
                    ' of ', int2str(codePeriods), ' msec',Ln...
                    'L1 C/No: ',CNoL1,' (dB-Hz),',Ln...
                    'L5 C/No: ',CNoL5,' (dB-Hz)'];

                try
                    waitbar(loopCnt/codePeriods,hwb,trackingStatus);
                catch %#ok<CTCH>
                    % The progress bar was closed. It is used as a signal
                    % to stop, "cancel" processing. Exit.
                    disp('Progress bar closed, exiting...');
                    return
                end
            end

            %% Read next block of data ------------------------------------------------
            % Record sample number (based on 8bit samples)
            trackResultsL1(channelNr).absoluteSample(loopCnt) =(ftell(fidL1))/dataAdaptCoeff;
            trackResultsL5(channelNr).absoluteSample(loopCnt) =(ftell(fidL5))/dataAdaptCoeff;


            %% L1 signal tracking loop
            % Update the phasestep based on code freq (variable) and sampling frequency (fixed)
            codePhaseStepL1 = codeFreqL1 / settings.samplingFreq;            
            % Find the size of a "block" or code period in whole samples
            blksizeL1 = ceil((settings.codeLength1-remCodePhaseL1) / codePhaseStepL1);
            % Read in the appropriate number of samples to process this interation
            [rawSignalL1, samplesReadL1] = fread(fidL1, dataAdaptCoeff*blksizeL1, settings.dataType);
            rawSignalL1 = rawSignalL1';
            % For complex data 
            if (dataAdaptCoeff==2)
                rawSignal1=rawSignalL1(1:2:end);
                rawSignal2=rawSignalL1(2:2:end);
                rawSignalL1 = rawSignal1 + 1i .* rawSignal2;  % transpose vector
            end
            % If did not read in enough samples, then could be out of data - better exit
            if (samplesReadL1 ~= dataAdaptCoeff*blksizeL1)
                disp('Not able to read the specified number of samples  for tracking, exiting!')
                delete(hwb);
                return
            end

            %% Set up all the code phase tracking information -------------------------
            % Save remCodePhase for current correlation
            trackResultsL1(channelNr).remCodePhase(loopCnt) = remCodePhaseL1;
            
            % Define index into early code vector
            tcode       = (remCodePhaseL1+earlyLateSpcL1) : codePhaseStepL1 : ((blksizeL1-1)*codePhaseStepL1+remCodePhaseL1+earlyLateSpcL1);
            tcode2      = ceil(tcode) + 1;
            earlyCode   = caCode(tcode2);
            % Define index into late code vector
            tcode       = (remCodePhaseL1-earlyLateSpcL1) : codePhaseStepL1 : ((blksizeL1-1)*codePhaseStepL1+remCodePhaseL1-earlyLateSpcL1);
            tcode2      = ceil(tcode) + 1;
            lateCode    = caCode(tcode2);
            % Define index into prompt code vector
            tcode       = remCodePhaseL1 : codePhaseStepL1 : ((blksizeL1-1)*codePhaseStepL1+remCodePhaseL1);
            tcode2      = ceil(tcode) + 1;
            promptCode  = caCode(tcode2);

            % Remaining code phase for each tracking update
            remCodePhaseL1 = (tcode(blksizeL1) + codePhaseStepL1) - settings.codeLength1;
            
            %% Generate the carrier frequency to mix the signal to baseband -----------
            % Save remCarrPhase for current correlation
            trackResultsL1(channelNr).remCarrPhase(loopCnt) = remCarrPhaseL1;

            % Get the argument to sin/cos functions
            time    = (0:blksizeL1) ./ settings.samplingFreq;
            trigarg = ((carrFreqL1 * 2.0 * pi) .* time) + remCarrPhaseL1;
            % Remaining carrier phase for each tracking update
            remCarrPhaseL1 = rem(trigarg(blksizeL1+1), (2 * pi));
            % Finally compute the signal to mix the collected data to bandband
            carrsig = exp(-1i .* trigarg(1:blksizeL1));
            %% Do correlation to Generate the six standard accumulated values -----------
            % First mix to baseband
            iBasebandSignal = real(carrsig .* rawSignalL1);
            qBasebandSignal = imag(carrsig .* rawSignalL1);
            % Now get early, late, and prompt values for each
            I_E = sum(earlyCode  .* iBasebandSignal);
            Q_E = sum(earlyCode  .* qBasebandSignal);
            I_P = sum(promptCode .* iBasebandSignal);
            Q_P = sum(promptCode .* qBasebandSignal);
            I_L = sum(lateCode   .* iBasebandSignal);
            Q_L = sum(lateCode   .* qBasebandSignal);
            %% Calculate PLL and DLL error --------------------------------------------
            carrErrorL1 = atan(Q_P / I_P) / (2.0 * pi);
            trackResultsL1(channelNr).carrFreq(loopCnt) = carrFreqL1;    
            codeErrorL1 = (1-earlyLateSpcL1)*(sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) / ...
                (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));
            trackResultsL1(channelNr).codeFreq(loopCnt) = codeFreqL1;

            %% Record various measures to show in postprocessing ----------------------
            trackResultsL1(channelNr).PRN     = channel(channelNr).PRN;
            trackResultsL1(channelNr).dllDiscr(loopCnt)       = codeErrorL1;
            trackResultsL1(channelNr).dllDiscrFilt(loopCnt)   = codeNcoL1;
            trackResultsL1(channelNr).pllDiscr(loopCnt)       = carrErrorL1;
            trackResultsL1(channelNr).pllDiscrFilt(loopCnt)   = carrNcoL1;
            trackResultsL1(channelNr).I_E(loopCnt) = I_E;
            trackResultsL1(channelNr).I_P(loopCnt) = I_P;
            trackResultsL1(channelNr).I_L(loopCnt) = I_L;
            trackResultsL1(channelNr).Q_E(loopCnt) = Q_E;
            trackResultsL1(channelNr).Q_P(loopCnt) = Q_P;
            trackResultsL1(channelNr).Q_L(loopCnt) = Q_L;

            %% L5 signal tracking loop -----------------------------------------------
            % Update the phasestep based on code freq (variable) and sampling frequency (fixed)
            codePhaseStepL5 = codeFreqL5 / settings.samplingFreq;
            % Find the size of a "block" or code period in whole samples
            blksizeL5 = ceil((settings.codeLength5-remCodePhaseL5) / codePhaseStepL5);
            % Read in the appropriate number of samples to process this interation
            [rawSignalL5, samplesReadL5] = fread(fidL5,dataAdaptCoeff*blksizeL5, settings.dataType);
            rawSignalL5 = rawSignalL5';
            if (dataAdaptCoeff==2)
                rawSignal1=rawSignalL5(1:2:end);
                rawSignal2=rawSignalL5(2:2:end);
                rawSignalL5 = rawSignal1 + 1i .* rawSignal2;  % transpose vector
            end
            % If did not read in enough samples, then could be out of data - better exit
            if (samplesReadL5 ~= dataAdaptCoeff*blksizeL5)
                disp('Not able to read the specified number of samples  for tracking, exiting!')
                delete(hwb);
                return
            end
            %% Set up all the code phase tracking information -------------------------
            % Save remCodePhase for current correlation
            trackResultsL5(channelNr).remCodePhase(loopCnt) = remCodePhaseL5;
            % Define index into early code vector
            tcode       = (remCodePhaseL5+earlyLateSpcL5) : codePhaseStepL5 : ((blksizeL5-1)*codePhaseStepL5+remCodePhaseL5+earlyLateSpcL5);
            tcode2      = ceil(tcode) + 1;
            earlyCode   = L5ICode(tcode2);            
            % For pilot channel signal tracking
            if (settings.pilotTRKflag == 1)
                earlyCodeQ   = L5QCode(tcode2);
            end
            % Define index into late code vector
            tcode       = (remCodePhaseL5-earlyLateSpcL5) : codePhaseStepL5 : ((blksizeL5-1)*codePhaseStepL5+remCodePhaseL5-earlyLateSpcL5);
            tcode2      = ceil(tcode) + 1;
            lateCode    = L5ICode(tcode2);           
            % For pilot channel signal tracking
            if (settings.pilotTRKflag == 1)
                lateCodeQ   = L5QCode(tcode2);
            end
            % Define index into prompt code vector
            tcode       = remCodePhaseL5 : codePhaseStepL5 : ((blksizeL5-1)*codePhaseStepL5+remCodePhaseL5);
            tcode2      = ceil(tcode) + 1;
            promptCode  = L5ICode(tcode2);          
            % For pilot channel signal tracking
            if (settings.pilotTRKflag == 1)
                promptCodeQ   = L5QCode(tcode2);
            end

            remCodePhaseL5 = (tcode(blksizeL5) + codePhaseStepL5) - settings.codeLength5;

            %% Generate the carrier frequency to mix the signal to baseband -----------           
            % Save remCarrPhase for current correlation
            trackResultsL5(channelNr).remCarrPhase(loopCnt) = remCarrPhaseL5;
            % Get the argument to sin/cos functions
            time    = (0:blksizeL5) ./ settings.samplingFreq;
            trigarg = ((carrFreqL5 * 2.0 * pi) .* time) + remCarrPhaseL5;
            remCarrPhaseL5 = rem(trigarg(blksizeL5+1), (2 * pi));
            % Finally compute the signal to mix the collected data to bandband
            carrsig = exp(-1i .* trigarg(1:blksizeL5));
            %% Generate the six standard accumulated values ---------------------------
            % First mix to baseband
            iBasebandSignal = real(carrsig .* rawSignalL5);
            qBasebandSignal = imag(carrsig .* rawSignalL5);

            % Now get early, late, and prompt values for each
            I_E = sum(earlyCode  .* iBasebandSignal);
            Q_E = sum(earlyCode  .* qBasebandSignal);
            I_P = sum(promptCode .* iBasebandSignal);
            Q_P = sum(promptCode .* qBasebandSignal);
            I_L = sum(lateCode   .* iBasebandSignal);
            Q_L = sum(lateCode   .* qBasebandSignal);
            
            % For pilot channel signal tracking
            if (settings.pilotTRKflag == 1)
                I_EQ = sum(earlyCodeQ  .* iBasebandSignal);
                Q_EQ = sum(earlyCodeQ  .* qBasebandSignal);
                I_PQ = sum(promptCodeQ .* iBasebandSignal);
                Q_PQ = sum(promptCodeQ .* qBasebandSignal);
                I_LQ = sum(lateCodeQ   .* iBasebandSignal);
                Q_LQ = sum(lateCodeQ   .* qBasebandSignal);
            end

            %% Calculate PLL and DLL error --------------------------------------------
            carrErrorL5 = atan(Q_P / I_P) / (2.0 * pi);
            if (settings.pilotTRKflag == 1)                
                QI = (I_PQ + 1i*Q_PQ) * exp(-1i*pi/2);
                carrErrorQ = atan(imag(QI)/real(QI))/ (2.0 * pi);
                carrErrorL5 = (carrErrorL5 + carrErrorQ)/2;
            end
            trackResultsL5(channelNr).carrFreq(loopCnt) = carrFreqL5;
            codeErrorL5 = (1-earlyLateSpcL5)*(sqrt(I_E * I_E + Q_E * Q_E) - sqrt(I_L * I_L + Q_L * Q_L)) / ...
                (sqrt(I_E * I_E + Q_E * Q_E) + sqrt(I_L * I_L + Q_L * Q_L));
            if (settings.pilotTRKflag == 1)
                codeErrorQ = (1-earlyLateSpcL5)*(sqrt(I_EQ^2 + Q_EQ^2) - sqrt(I_LQ^2 + Q_LQ^2)) / ...
                             (sqrt(I_EQ^2 + Q_EQ^2) + sqrt(I_LQ^2 + Q_LQ^2));
                codeErrorL5 = (codeErrorL5 + codeErrorQ)/2;
            end
            trackResultsL5(channelNr).codeFreq(loopCnt) = codeFreqL5;


            %% Kalman filter predict and update
            Z = [codeErrorL1;carrErrorL1;codeErrorL5;carrErrorL5]; % Measurement vector
            state_cov = Transistion_Matrix * state_cov * transpose(Transistion_Matrix) + process_noise;  
            kalman_gain = state_cov * transpose(Measure_Matrix) * inv(Measure_Matrix * state_cov * transpose(Measure_Matrix) +  mesurement_noise);         
            error_state = error_state + kalman_gain*(Z-Measure_Matrix*error_state+settings.intTime/2*[codeNcoL1;carrNcoL1;codeNcoL5;carrNcoL5]);
            state_cov = (eye(num_state) - kalman_gain * Measure_Matrix) * state_cov;
            error_state = Transistion_Matrix * error_state - [settings.intTime;settings.intTime;settings.intTime;settings.intTime;0;0].*[codeNcoL1;codeNcoL5;carrNcoL1;carrNcoL5;0;0] ;

            %% Frequency backfeed to NCO
            codeNcoL1 = error_state(1)/settings.intTime+beta1*error_state(5) ;
            carrNcoL1 = error_state(3)/settings.intTime+alpha*error_state(5)+alpha*error_state(6)*settings.intTime/2;
            codeFreqL1 = settings.codeFreqBasis1 + codeNcoL1;
            carrFreqL1 = carrFreqBasisL1 + carrNcoL1;
            codeNcoL5 = error_state(2)/settings.intTime + beta2*error_state(5) ;
            carrNcoL5 = error_state(4)/settings.intTime + error_state(5) + error_state(6)*settings.intTime/2;
            codeFreqL5 = settings.codeFreqBasis5 + codeNcoL5;
            carrFreqL5 = carrFreqBasisL5 + carrNcoL5;

            %% Record various measures to show in postprocessing ----------------------
            trackResultsL5(channelNr).PRN                     = channel(channelNr).PRN;
            trackResultsL5(channelNr).dllDiscr(loopCnt)       = codeErrorL5;
            trackResultsL5(channelNr).dllDiscrFilt(loopCnt)   = codeNcoL5;
            trackResultsL5(channelNr).pllDiscr(loopCnt)       = carrErrorL5;
            trackResultsL5(channelNr).pllDiscrFilt(loopCnt)   = carrNcoL5;
            trackResultsL5(channelNr).I_E(loopCnt) = I_E;
            trackResultsL5(channelNr).I_P(loopCnt) = I_P;
            trackResultsL5(channelNr).I_L(loopCnt) = I_L;
            trackResultsL5(channelNr).Q_E(loopCnt) = Q_E;
            trackResultsL5(channelNr).Q_P(loopCnt) = Q_P;
            trackResultsL5(channelNr).Q_L(loopCnt) = Q_L;
            if (settings.pilotTRKflag == 1)
                trackResultsL5(channelNr).Pilot_I_P(loopCnt) = I_PQ ;
                trackResultsL5(channelNr).Pilot_Q_P(loopCnt) = Q_PQ;
            end


            %% CNo calculation --------------------------------------------------------
            if (rem(loopCnt,settings.CNo.VSMinterval)==0)
                % L1
                vsmCntL1=vsmCntL1+1;
                CNoValueL1=CNoVSM(trackResultsL1(channelNr).I_P(loopCnt-settings.CNo.VSMinterval+1:loopCnt),...
                    trackResultsL1(channelNr).Q_P(loopCnt-settings.CNo.VSMinterval+1:loopCnt),settings.CNo.accTime);
                trackResultsL1(channelNr).CNo.VSMValue(vsmCntL1)=CNoValueL1;
                trackResultsL1(channelNr).CNo.VSMIndex(vsmCntL1)=loopCnt;
                CNoL1=int2str(CNoValueL1);
                % update measure noise vector
                tempCN0L1 = 10^(CNoValueL1/10);
                tau_varL1 = (earlyLateSpcL1/(tempCN0L1*4*settings.intTime))*(1+2/((2-earlyLateSpcL1)*tempCN0L1*settings.intTime));
                phi_varL1 = (1/(tempCN0L1*2*settings.intTime))*(1+1/(2*tempCN0L1*settings.intTime));
                % L5
                vsmCntL5=vsmCntL5+1;
                CNoValueL5=CNoVSM(trackResultsL5(channelNr).I_P(loopCnt-settings.CNo.VSMinterval+1:loopCnt),...
                    trackResultsL5(channelNr).Q_P(loopCnt-settings.CNo.VSMinterval+1:loopCnt),settings.CNo.accTime);
                trackResultsL5(channelNr).CNo.VSMValue(vsmCntL5)=CNoValueL5;
                trackResultsL5(channelNr).CNo.VSMIndex(vsmCntL5)=loopCnt;
                CNoL5=int2str(CNoValueL5);
                % update measure noise vector
                tempCN0L5 = 10^(CNoValueL5/10);
                tau_varL5 = (earlyLateSpcL5/(tempCN0L5*4*settings.intTime))*(1+2/((2-earlyLateSpcL5)*tempCN0L5*settings.intTime));
                phi_varL5 = (1/(tempCN0L5*2*settings.intTime))*(1+1/(2*tempCN0L5*settings.intTime));
                mesurement_noise = diag([tau_varL1,phi_varL1,tau_varL5,phi_varL5]);
            end

        end % for loopCnt

        % If we got so far, this means that the tracking was successful
        % Now we only copy status, but it can be update by a lock detector
        % if implemented
        trackResultsL1(channelNr).status  = channel(channelNr).statusL1;
        trackResultsL5(channelNr).status  = channel(channelNr).statusL5;

    end % if a PRN is assigned
end % for channelNr

% Close the waitbar
close(hwb)
