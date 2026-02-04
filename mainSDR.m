%--------------------------------------------------------------------------
%                             MJTL_SDR v1.0
% Copyright (C) 2026 J. Sun, Z. Tang, J. Wei and J. Zhou.
% Developed for GPS L1/L5 multi-frequency joint tracking by J. Sun, 
% Z. Tang, J. Wei and J. Zhou.
% -------------------------------------------------------------------------

%% Clean up the environment first =========================================
clear; close all; clc;

format ('compact');
format ('long', 'g');

%--- Include folders with functions ---------------------------------------
addpath Common                % Common funSctions between differnt SDR receivers
addpath Setting               % Initialization functions
addpath Acquisition           % Acquisition functions
addpath Tracking              % Tracking functions
addpath Navigation            % Navigation decode functions
addpath Figplot               % Plot figure functions
addpath Data                  % GPS IF Data Folder

%% Print startup ==========================================================
fprintf('---------------------------------------------------------------------------------------');
fprintf(['\nWelcome to:  MJTL_SDR\n', ...
    'An open-source MATLAB code for GPS L1+L5 joint tracking\n', ...
    'The code was improved by Huazhong University of Science and Technology.\n']);
fprintf('---------------------------------------------------------------------------------------\n\n');

%% Initialize constants, settings =========================================
sys = input(['Choosing SDR mode: "1" -- L1 standard tracking Loop \n' ...
    '                   "5" -- L5 standard tracking Loop\n ' ...
    '                  "6" -- Doppler‑aided multi-frequency tracking Loop (DAJTL) \n'...
    '                   "7" -- Kalman filter-based multi-frequency tracking Loop (KFJTL):\n' ...
    'Input: ']);  % input
% sys = 6;   % 1:L1 5:L5 6:L1+L5 DAJTL 7:L1+L5 KFJTL 

switch sys
    case 1
        settings = initSettingsL1();
        % fprintf('Probing data (%s)...\n', settingsL1.fileName);
        % probeData(settings);  % Generate plot of raw data
        % disp('  Raw IF data plotted ');
    case 5
        settings = initSettingsL5();
        % fprintf('Probing data (%s)...\n', settingsL5.fileName);
        % probeData(settings);    % Generate plot of raw data
        % disp('  Raw IF data plotted ');
    case {6,7}
        settingsL1 = initSettingsL1();
        settingsL5 = initSettingsL5();
        settings   = initSettingsDual(settingsL1,settingsL5);
        % fprintf('Probing data (%s) and (%s)...\n', settingsL1.fileName, settingsL5.fileName);
        % probeData(settingsL1);
        % probeData(settingsL5);
        % disp('  Raw IF data plotted ');
    otherwise
        error('Unexpected input.');
end
    
%% Acquisition ============================================================
disp ('Starting acquiring satellites...');
if ~exist(['./Result/acqResults_',settings.fileName,'.mat'])
    switch sys
        case 1 
            acqResultsL1 = acquisitionL1(settings);
            save(['./Result/','acqResults_',settings.fileName,],'acqResultsL1');
        case 5
            acqResultsL5 = acquisitionL5(settings);
            save(['./Result/','acqResults_',settings.fileName,],'acqResultsL5');
        case {6,7}
            [acqResultsL1,acqResultsL5] = acquisitionDual(settingsL1,settingsL5);
            save(['./Result/','acqResults_',settings.fileName,],"acqResultsL1","acqResultsL5");
        case 8
            [acqResultsL1,acqResultsL5] = acquisitionDual(settingsL1,settingsL5);
            save(['./Result/','acqResults_',settings.fileName,],"acqResultsL1","acqResultsL5");
        otherwise
            disp('error warning');
    end
else
    load(['./Result/','acqResults_',settings.fileName,'.mat']);
end
    
%% Initialize channels and prepare for the run ============================
switch sys
    case 1 
        channelL1 = preRun(acqResultsL1, settings);
        showChannelStatus(channelL1, settings);
    case 5
        channelL5 = preRun(acqResultsL5, settings);
        showChannelStatus(channelL5, settings);
    case {6,7}
        channelL1 = preRun(acqResultsL1, settingsL1);
        channelL5 = preRun(acqResultsL5, settingsL5);
        vissat1 = [channelL1.PRN];vissat5 = [channelL5.PRN];
        % Get the common visiable satellite
        vissat = intersect(vissat1, vissat5);vissat = vissat(vissat ~= 0);
        for sn = 1:length(vissat)
            channel(sn).PRN = vissat(sn);
            Indices = find([channelL1.PRN] == vissat(sn));
            channel(sn).acquiredFreqL1 = channelL1(Indices).acquiredFreq;
            channel(sn).codePhaseL1 = channelL1(Indices).codePhase;
            channel(sn).statusL1 = channelL1(Indices).status;
            channel(sn).codeFreqL1 = channelL1(Indices).codeFreq;
            Indices = find([channelL5.PRN] == vissat(sn));
            channel(sn).acquiredFreqL5 = channelL5(Indices).acquiredFreq;
            channel(sn).codePhaseL5 = channelL5(Indices).codePhase;
            channel(sn).statusL5 = channelL5(Indices).status;
            channel(sn).codeFreqL5 = channelL5(Indices).codeFreq;
        end
        showDualChannelStatus(channel, settings);
    otherwise
        disp('error warning');
end

    
%% Track the signal =======================================================
startTime = now;
disp (['   Tracking started at ', datestr(startTime)]);
if (~exist(['./Result/trkResults_',settings.fileName,'.mat']))    
    switch sys
        case 1 
            % L1 single-frequency tracking loop
            [trkResultsL1] = trackingL1(channelL1, settings);
            save(['./Result/','trkResults_',settings.fileName,],'trkResultsL1');
        case 5
            % L5 single-frequency tracking loop
            [trkResultsL5] = trackingL5(channelL5, settings);
            save(['./Result/','trkResults_',settings.fileName,],'trkResultsL5');
        case 6
            % Doppler-aided L1/L5 joint tracking loop
            [trkResultsL1_DA,trkResultsL5_DA] = trackingDualDopAid(channel, settings);
            save(['./Result/','trkResults_DA_',settings.fileName,],'trkResultsL1_DA','trkResultsL5_DA');
        case 7
            % KF-based L1/L5 joint tracking loop
            [trkResultsL1_KF,trkResultsL5_KF] = trackingDualKF(channel, settings);
            save(['./Result/','trkResults_KF_',settings.fileName,],'trkResultsL1_KF','trkResultsL5_KF');
        otherwise
            disp('error warning');
    end
else
    load(['./Result/','trkResults_',settings.fileName,'.mat']);
end
disp(['   Tracking is over (elapsed time ', datestr(now - startTime, 13), ')'])                      
    
%% Calculate navigation solutions =========================================
disp('   Calculating navigation solutions...');
switch sys
    case 1 
        [navResultsL1, ephL1] = postNavigationL1(trkResultsL1, settings);
    case 5
        [navResultsL5, ephL5] = postNavigationL5(trkResultsL5, settings);
    case 6
        [navResultsL1_DA, ephL1] = postNavigationL1(trkResultsL1_DA, settingsL1);
        [navResultsL5_DA, ephL5] = postNavigationL5(trkResultsL5_DA, settingsL5);
    case 7
        [navResultsL1_KF, ephL1] = postNavigationL1(trkResultsL1_KF, settingsL1);
        [navResultsL5_KF, ephL5] = postNavigationL5(trkResultsL5_KF, settingsL5);
    otherwise
        disp('error warning');
end
disp('Post processing of the signal is over.');
    
%% Plot all results ===================================================
disp ('   Ploting results...');

% Plot acquisition results
if settings.plotAcquisition && (sys == 1|| sys==6 || sys==7)
    plotAcquisition(acqResultsL1);
end
if settings.plotAcquisition && (sys == 5|| sys==6 || sys==7)
    plotAcquisition(acqResultsL5);
end

% Plot tracking results
if settings.plotTracking && sys == 1
    plotTracking(trkResultsL1, settings);
end
if settings.plotTracking && sys == 5
    plotTracking(trkResultsL5, settings);
end
if settings.plotTracking && sys == 6
    plotTracking(trkResultsL1_DA, settingsL1);
    plotTracking(trkResultsL5_DA, settingsL5);
end
if settings.plotTracking && sys == 7
    plotTracking(trkResultsL1_KF, settingsL1);
    plotTracking(trkResultsL5_KF, settingsL5);
end

% Plot navigation results
if settings.plotNavigation && sys == 1
    plotNavigation(navResultsL5, settings);
end 
if settings.plotNavigation && sys == 5
    plotNavigation(navResultsL5, settings);
end
if settings.plotNavigation && sys == 6
    plotNavigation(navResultsL1_DA, settingsL1);
    plotNavigation(navResultsL5_DA, settingsL5);
end
if settings.plotNavigation && sys == 7
    plotNavigation(navResultsL1_KF, settingsL1);
    plotNavigation(navResultsL5_KF, settingsL5);
end

disp('Post processing of the signal is over.');