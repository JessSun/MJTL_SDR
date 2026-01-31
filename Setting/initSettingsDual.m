function settings = initSettingsDual(settings1,settings5)
% Functions initializes and saves settings. Settings can be edited inside of
% the function, updated from the command line or updated using a dedicated
% GUI - "setSettings".  
%
% settings = initSettingsL1()
%
%   Inputs: none
%       settings1    - Receiver settings of L1 signal processing.
%       settings5    - Receiver settings of L5 signal processing.
%   Outputs:
%       settings     - Receiver settings of dual-frequency signals processing(a structure). 

%--------------------------------------------------------------------------
%                             MJTL_SDR v1.0
% Copyright (C) 2026 J. Sun, Z. Tang, J. Wei and J. Zhou.
% Developed for GPS L1/L5 multi-frequency joint tracking by J. Sun, 
% Z. Tang, J. Wei and J. Zhou.
% -------------------------------------------------------------------------

%% Processing settings ====================================================
settings.st1 = settings1;
settings.st5 = settings5;
settings.msToProcess        = settings1.msToProcess;                % [ms]
settings.numberOfChannels   = settings1.numberOfChannels;
settings.skipNumberOfBytes  = settings1.skipNumberOfBytes;          % [byte]

%% Raw data settings ======================================================
settings.fileName           = settings1.fileName(ismember(settings1.fileName, settings5.fileName));
settings.fileRoute1         = settings1.fileRoute;
settings.fileRoute5         = settings5.fileRoute;
settings.dataType           = settings1.dataType;          
settings.fileType           = settings1.fileType;

%% Signal settings ======================================================== 
settings.IF1                = settings1.IF;                         % [Hz]
settings.IF5                = settings5.IF;                         % [Hz]
settings.samplingFreq       = settings1.samplingFreq;               % [Hz]
settings.codeLength1        = settings1.codeLength;
settings.codeLength5        = settings5.codeLength;
settings.codeFreqBasis1     = settings1.codeFreqBasis;              % [Hz]
settings.codeFreqBasis5     = settings5.codeFreqBasis;              % [Hz]
settings.carrFreqBasis1     = settings1.carrFreqBasis;              % [Hz]
settings.carrFreqBasis5     = settings5.carrFreqBasis;              % [Hz]

%% Acquisition settings ===================================================
% Enable/disable use of pilot channel for tracking
settings.pilotTRKflag        = settings5.pilotTRKflag;              % 0 - Off
% Threshold for the signal presence decision rule
settings.acqThreshold1       = settings1.acqThreshold;
settings.acqThreshold5       = settings5.acqThreshold;

%% Tracking loops settings ================================================
% L1
settings.dllDampingRatio1         = settings1.dllDampingRatio;
settings.dllNoiseBandwidth1       = settings1.dllNoiseBandwidth;        % [Hz]
settings.dllCorrelatorSpacing1    = settings1.dllCorrelatorSpacing;     % [chips]
settings.pllDampingRatio1         = settings1.pllDampingRatio;
settings.pllNoiseBandwidth1       = settings1.pllNoiseBandwidth;        %[Hz]
settings.intTime                  = settings1.intTime;                  %[s]
% L5
settings.dllDampingRatio5         = settings5.dllDampingRatio;
settings.dllNoiseBandwidth5       = settings5.dllNoiseBandwidth;        %[Hz]  
settings.dllCorrelatorSpacing5    = settings5.dllCorrelatorSpacing;     %[chips]
settings.pllDampingRatio5         = settings5.pllDampingRatio;
settings.pllNoiseBandwidth5       = settings5.pllNoiseBandwidth;        %[Hz]
% CNo Settings
settings.CNo.VSMinterval = settings1.CNo.VSMinterval;
settings.CNo.accTime = settings1.CNo.accTime;

%% Navigation solution settings ===========================================
settings.navSolPeriod       = settings1.navSolPeriod;          %[ms]
settings.elevationMask      = settings1.elevationMask;         %[degrees 0 - 90]
settings.useTropCorr        = settings1.useTropCorr;
settings.truePosition.E     = settings1.truePosition.E;
settings.truePosition.N     = settings1.truePosition.N;
settings.truePosition.U     = settings1.truePosition.U;

%% Plot settings ==========================================================
% Enable/disable plotting of the tracking results for each channel
settings.plotTracking    = settings1.plotTracking;                   % 0 - Off; 1 - On
settings.plotAcquisition = settings1.plotAcquisition;
settings.plotNavigation  = settings1.plotAcquisition;

%% Constants ==============================================================
settings.c                  = 299792458;    % The speed of light, [m/s]
settings.startOffset        = settings1.startOffset;       %[ms] Initial sign. travel time


