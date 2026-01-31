function caCodesTable = makeCaTable(PRN,settings)
% Function generates CA codes for all 32 satellites based on the settings
% provided in the structure "settings". The codes are digitized at the
% sampling frequency specified in the settings structure.
% One row in the "caCodesTable" is one C/A code. The row number is the PRN
% number of the C/A code.
%
% caCodesTable = makeCaTable(settings)
%
%   Inputs:
%       PRN             - specified PRN for C/A code
%       settings        - receiver settings
%   Outputs:
%       caCodesTable    - an array of arrays (matrix) containing C/A codes
%                       for all satellite PRN-s
%--------------------------------------------------------------------------

% Find number of samples per spreading code 
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));

% Find time constants 
ts = 1/settings.samplingFreq;   % Sampling period in sec
tc = 1/settings.codeFreqBasis;  % C/A chip period in sec

%  Generate CA code for given PRN 
caCode = generateCAcode(PRN);

%% Digitizing
%  Make index array to read C/A code values 
% The length of the index array depends on the sampling frequency -
% number of samples per millisecond (because one C/A code period is one
% millisecond).
codeValueIndex = ceil((ts * (1:samplesPerCode)) / tc);

% Correct the last index (due to number rounding issues)
codeValueIndex(end) = 1023;

% Make the digitized version of the C/A code 
% The "upsampled" code is made by selecting values form the CA code
% chip array (caCode) for the time instances of each sample.
caCodesTable = caCode(codeValueIndex);
