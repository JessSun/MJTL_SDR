function [CNodB]= CNoVSM(I,Q,T)
%Calculate CNo using the Variance Summing Method
%
%[CNodB]= CNoVSM(I,Q,T)
%
%   Inputs:
%       I           - Prompt In Phase values of the signal from Tracking
%       Q           - Prompt Quadrature Phase values of the signal from Tracking
%       T          - Accumulation interval in Tracking (in sec)
%   Outputs:
%       CNo         - Estimated C/No for the given values of I and Q

%--------------------------------------------------------------------------
%                             MJTL_SDR v1.0
% Copyright (C) 2026 J. Sun, Z. Tang, J. Wei and J. Zhou.
% Developed for GPS L1/L5 multi-frequency joint tracking by J. Sun, 
% Z. Tang, J. Wei and J. Zhou.
% -------------------------------------------------------------------------

%Calculate Power
Z=I.^2+Q.^2;
%Calculate the mean and variance of the Power
Zm=mean(Z);
Zv=var(Z);
%Calculate the average carrier power
Pav=sqrt(Zm^2-Zv);
%Calculate the variance of the noise
Nv=0.5*(Zm-Pav);
%Calculate C/No
CNodB=10*log10(abs((1/T)*Pav/(2*Nv)));

% if length(I)==100
%     K=5;
%     M=20;
%     for i = 1:K
%         Pwb(i) = sum(   I((i-1)*M+1 : i*M).^2   +   Q((i-1)*M+1 : i*M).^2   );
%     end
%     for i = 1:K
%         Pnb(i) = sum(abs(I((i-1)*M+1 : i*M))).^2 +...
%             sum(abs(Q((i-1)*M+1 : i*M))).^2;
%     end
%     Pnw = Pnb./Pwb;
%     miuP = 1/K*sum(Pnw);
%     CNo = abs((1/T)*(miuP-1)/(M-miuP));
%     CNodB = 10*log10(CNo);
% end