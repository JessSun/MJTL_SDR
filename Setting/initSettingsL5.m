function settings = initSettingsL5()
% Functions initializes and saves settings. Settings can be edited inside of
% the function, updated from the command line or updated using a dedicated
% GUI - "setSettings".  
%
% settings = initSettingsL5()
%
%   Inputs: none
%
%   Outputs:
%       settings     - Receiver settings of L5 signal processing (a structure). 

%--------------------------------------------------------------------------
%                             MJTL_SDR v1.0
% Copyright (C) 2026 J. Sun, Z. Tang, J. Wei and J. Zhou.
% Developed for GPS L1/L5 multi-frequency joint tracking by J. Sun, 
% Z. Tang, J. Wei and J. Zhou.
% -------------------------------------------------------------------------

%% Processing settings ====================================================
settings.msToProcess        = 60000;            % [ms]
settings.numberOfChannels   = 8;
settings.skipNumberOfBytes  = 0;                % [byte]

%% Raw data settings ======================================================
settings.fileName       = 'L5_IQ_Round3'; 
settings.fileRoute      = ['./Data/',settings.fileName,'.bin']; 
settings.dataType       = 'int8';               % Raw data type: 'int8','int16'
settings.fileType       = 2;                    % File Types: 1 - I, 2 - I/Q                       

%% Signal settings ======================================================== 
settings.IF                 = 0e6;              % [Hz]
settings.samplingFreq       = 30690000;         % [Hz]
settings.carrFreqBasis      = 1176.45e6;        % [Hz]
settings.codeFreqBasis      = 10.23e6;          % [Hz]
settings.codeLength         = 10230;

%% Acquisition settings ===================================================
settings.acqSatelliteList   = [1:32];           % [PRN numbers]
settings.acqSearchBand      = 5000;             % [Hz]
settings.acqNonCohTime      = 25;               % [ms]
settings.acqThreshold       = 4.5;
settings.acqSearchStep      = 500;              % [Hz]

%% Tracking loops settings ================================================
settings.pilotTRKflag            = 0;           % pilot channel enable: 0 - Off, 1 - On
% Code tracking loop parameters
settings.dllDampingRatio         = 0.7;
settings.dllNoiseBandwidth       = 2;           % [Hz]  
settings.dllCorrelatorSpacing    = 0.5;         % [chips]
% Carrier tracking loop parameters
settings.pllDampingRatio         = 0.7;
settings.pllNoiseBandwidth       = 15;          % [Hz]
settings.intTime                 = 0.001;       % [s]
% CNo Settings
settings.CNo.accTime             = 0.001;       % [s]
settings.CNo.VSMinterval         = 400;         % [ms]

%% Navigation solution settings ===========================================
settings.navSolPeriod       = 500;              % [ms]
settings.elevationMask      = 5;                % [degrees 0 - 90]
settings.useTropCorr        = 1;                % tropospheric correction enable: 0 - Off, 1 - On
settings.truePosition.E     = nan;
settings.truePosition.N     = nan;
settings.truePosition.U     = nan;

%% Constants ==============================================================
settings.c                  = 299792458;        % [m/s]
settings.startOffset        = 68.802;           % [ms]

%% Plot settings ==========================================================
% Enable/disable plotting of the tracking results for each channel
settings.plotTracking    = 1;                   % 0 - Off; 1 - On
settings.plotAcquisition = 1;
settings.plotNavigation  = 1;

