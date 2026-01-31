function [pf3,pf2,pf1] = calcLoopCoefCarr(settings,LBW)
%Function finds loop coefficients. The coefficients are used then in PLL-s
%and DLL-s.
%
%[tau1, tau2] = calcLoopCoef(LBW, zeta, k)
%
%   Inputs:
%       settings     - Receiver settings
%       LBW          - Loop noise bandwidth
%
%   Outputs:
%       tau1, tau2   - Loop filter coefficients 

%--------------------------------------------------------------------------
%                             MJTL_SDR v1.0
% Copyright (C) 2026 J. Sun, Z. Tang, J. Wei and J. Zhou.
% Developed for GPS L1/L5 multi-frequency joint tracking by J. Sun, 
% Z. Tang, J. Wei and J. Zhou.
% -------------------------------------------------------------------------

% Summation interval
intTime = settings.intTime;

% loop constant coefficients
a3 = 2;
b3 = 2;

% Solve natural frequency
Wn = 1.2 * LBW;

% solve for [pf3,pf2,pf1]
pf3 = Wn^3 * intTime^2;
pf2 = a3 * Wn^2 * intTime;
pf1 = b3 * Wn;


