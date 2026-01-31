function L5ICodesTable = makeL5ITable(PRN,settings)
% Function generates L5I codes for all 32 satellites based on the settings
% provided in the structure "settings". The codes are digitized at the
% sampling frequency specified in the settings structure.
% One row in the "L5ICodesTable" is one L5I code. The row number is the PRN
% number of the L5I code.
%
% L5ICodesTable = makeL5ITable(settings)
%
%   Inputs:
%       PRN             - specified PRN for L5I code
%       settings        - receiver settings
%   Outputs:
%       L5ICodesTable   - an array of arrays (matrix) containing L5I codes
%                       for all satellite PRN-s
%--------------------------------------------------------------------------

% Find number of samples per spreading code
samplesPerCode = round(settings.samplingFreq / ...
    (settings.codeFreqBasis / settings.codeLength));

% Find time constants
ts = 1/settings.samplingFreq;   % Sampling period in sec
tc = 1/settings.codeFreqBasis;  % L5I chip period in sec

% Generate L5I code for given PRN 
L5ICode = generateL5Icode(PRN,settings);

%% Digitizing
%  Make index array to read L5I code values 
codeValueIndex = ceil((ts * (1:samplesPerCode)) / tc);

% Correct the last index (due to number rounding issues) 
codeValueIndex(end) = settings.codeLength;

%  Make the digitized version of the L5I code 
% The "upsampled" code is made by selecting values form the L5I code
% chip array (caCode) for the time instances of each sample.
L5ICodesTable = L5ICode(codeValueIndex);

