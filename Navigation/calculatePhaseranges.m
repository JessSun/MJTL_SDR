function [Phaseranges] = calculatePhaseranges(settings,integrated_doppler,subFrameStart,currMeasNr,channelList)
% CalculatePhaseranges, generates the ambiguous phase ranges from the integrated doppler measurements.
%
% [Phaseranges] = calculatePhaseranges(settings,integrated_doppler,msOfTheSignal,channelList,navSolutionPeriod)
%
% Inputs:
%       settings            - receiver settings
%       integrated_doppler  - Integrated Doppler Measurements
%       subFrameStart       - the array contains positions of the first
%                           preamble in each channel. The position is ms count 
%                           since start of tracking. Corresponding value will
%                           be set to 0 if no valid preambles were detected in
%                           the channel: 
%                           1 by settings.numberOfChannels
%       currMeasNr          - current measurement sample location
%       channelList         - number of active channels
%
% Outputs:
%       Phaseranges         - Ambiguous phase ranges of the detected satellite signals 

%--------------------------------------------------------------------------
%                             MJTL_SDR v1.0
% Copyright (C) 2026 J. Sun, Z. Tang, J. Wei and J. Zhou.
% Developed for GPS L1/L5 multi-frequency joint tracking by J. Sun, 
% Z. Tang, J. Wei and J. Zhou.
% -------------------------------------------------------------------------

% Pseudorange measurement point (millisecond) in the trackResults structure
msOfTheSignal = subFrameStart + (settings.navSolPeriod/(settings.intTime*1000)) * (currMeasNr-1);
% Signal wavelength
lambda = settings.c/settings.carrFreqBasis;

% For all channels in the list ...
for channelNr = channelList
    % Ambiguous phase ranges
    Phaseranges(channelNr) = integrated_doppler(channelNr,msOfTheSignal(channelNr))*lambda*settings.navSolPeriod*1e-3;   % in meters/sec 

end

