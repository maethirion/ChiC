% This computes the alpha and gamma for The Chicago (Code) Construction,
% or, the ChiC, for an (n,k,d) DSS when the inner code (post-MRD) being
% used is an (n,kIC,d) DSS.
% 
% The alpha and gamma are ultimately plotted for a file size of 1 unit; the
% plot can be assumed to be between alpha per file size vs gamma per file
% size. 
% 
% To understand the use of MRD (maximum rank distance) codes in the
% construction, please refer to "New Codes and Inner Bounds for Exact
% Repair in Distributed Storage Systems" by Goparaju et al. (ISIT 2014),
% and "High-Rate Regenerating Codes Through Layering" by Sasidharan et al. 
% (http://arxiv.org/abs/1301.6157). The reference to Ernvall refers to
% "Exact-Regenerating Codes between MBR and MSR Points" by Ernvall
% (http://arxiv.org/abs/1304.5357).
% 
% File created by: Sreechakra Goparaju (SG)
% File created on: Jan 31, 2014
% File last modified on: Jul 14, 2014
% File last modified by: SG
%
% @inputs:  n = number of nodes in the DSS
%           k = number of "information" nodes in the DSS
%           d = number of helper nodes for a single failure repair
%           kIC = the "working" k of the Inner Code (code after MRD);
%                 for Ernvall, k_inner = k; for Sasidharan et al, k_inner = d.
% @outputs: 
%           tradeoffCurve = the alphaBar vs gammaBar tradeoff curve
% $other vars and notation:
%           alphaBar refers to alpha/M (M is the file size)
%           gammaBar refers to gamma/M
%           MSR and MBR refer to the MSR and the MBR points.


function [tradeoffCurve] = ChicagoConstruction(n,k,d,kIC)

% wMax = floor((d*k-kIC*(k-1))/(d-k+1));
wMax = k;
% Refer to Goparaju et al., ISIT 2014.
% This gives the maximum value for kHat: the "k" for the small code. Any
% other value beyond that leads to a point beyond MSR, and so needn't be
% calculated (IF UNIVERSAL/OPPORTUNISTIC REPAIR IS NOT USED - So, ideally,
% the value of wMax can be larger. But we work with this for now and remove
% the restriction later).
% Remark (Feb 03, 2014; SG): For (61,55,59), removing the restriction 
% didn't change anything in the regime which matters (between MBR and MSR).

tradeoffCurve = zeros(wMax,2);

for i = 1:wMax,
%     i
    kHat = i;
    nHat = i + (n-kIC);
    
    % calculations needed to compute alphaBar
    MAux = 0;
    for j = max([1,k-(n-nHat)]):min([k,nHat]),
        MAux = MAux + (min([kHat,j])*nchoosek(nHat,j)*nchoosek(n-nHat,k-j));
    end
    MAux = MAux/nchoosek(n,k);
    alphaAux = nHat/n;
    
    % alphaBar when kHat = i
    tradeoffCurve(i,1) = alphaAux/MAux;
    
    % calculations needed to compute gammaBar
    gammaAux = 0;
    for j = d-(n-nHat):min([d,nHat-1]),
        gammaAux = gammaAux + ((j/(j-kHat+1))*nchoosek(n-nHat,d-j)*nchoosek(nHat-1,j));
    end
    gammaAux = gammaAux/nchoosek(n-1,d);
    
    % gammaBar when kHat = i 
    tradeoffCurve(i,2) = gammaAux*tradeoffCurve(i,1);
    
end

% % THESE GO IN THE OUTER CODE! (ISITPoints.m)
% 
% alphaBarMSR = 1/k;
% gammaBarMSR = d/(k*(d-k+1));
% 
% alphaBarMBR = d/(d*k-nchoosek(k,2));
% gammaBarMBR = alphaBarMBR;