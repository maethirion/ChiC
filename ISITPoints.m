% This computes the alpha and gamma for the code constructions in "New
% Codes and Inner Bounds for Exact Repair in Distributed Storage Systems"
% co-authored by Sreechakra Goparaju, Salim El Rouayheb, and Robert
% Calderbank (presented at ISIT 2014), and beyond.
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
% There are two main constructions in the paper and the code. Both use an
% existing construction overloaded with the concept of opportunistic
% repair. 
% For Construction 1, refer to "Exact-Regenerating Codes between
% MBR and MSR Points" by Ernvall (http://arxiv.org/abs/1304.5357).
% For Construction 2, refer to "High-Rate Regenerating Codes Through
% Layering" by Sasidharan and Kumar (http://arxiv.org/abs/1301.6157).
% For opportunistic repair, refer to "Distributed Data Storage Systems with
% Opportunistic Repair" by Aggarwal et al (http://arxiv.org/abs/1311.4096)
% and "Asymptotic Interference Alignment for Optimal Repair of MDS codes in
% Distributed Data Storage" by Cadambe et al
% (http://www.mit.edu/~viveck/resources/Research/asymptotic_storage.pdf).
% 
% File created by: Sreechakra Goparaju (SG)
% File created on: Feb 03, 2014
% File last modified on: Jul 14, 2014
% File last modified by: SG
% 
% @inputs:  n = number of nodes in the DSS
%           k = number of "information" nodes in the DSS
%           d = number of helper nodes for a single failure repair
% @outputs: 
% @other vars and notation:
%           alphaBar refers to alpha/M (M is the file size)
%           gammaBar refers to gamma/M
%           MSR and MBR refer to the MSR and the MBR points.

function [] = ISITPoints(n,k,d)

% Calculating the alphaBar and gammaBar for MSR and MBR
alphaBarMSR = 1/k;
gammaBarMSR = d/(k*(d-k+1));

alphaBarMBR = d/(d*k-nchoosek(k,2));
gammaBarMBR = alphaBarMBR;

% Construction 1 (à la Toni)
tradeoffCurve1 = ChicagoConstruction(n,k,d,k);
tC1aux = [alphaBarMBR, gammaBarMBR;
        tradeoffCurve1;
       alphaBarMSR, gammaBarMSR];
tC1 = convhull(tC1aux(:,1),tC1aux(:,2));

% Construction 2 (à la Birenjith)
tradeoffCurve2 = ChicagoConstruction(n,k,d,d);
tC2aux = [alphaBarMBR, gammaBarMBR;
        tradeoffCurve2;
       alphaBarMSR, gammaBarMSR];
tC2 = convhull(tC2aux(:,1),tC2aux(:,2));

% plot(tC1aux(tC1,1),tC1aux(tC1,2),'-ro');
% hold on;
% axis([alphaBarMSR alphaBarMBR gammaBarMBR gammaBarMSR])
% plot(tC2aux(tC2,1),tC2aux(tC2,2),'-bo');

% What lies between Constructions 1 and 2? (NOTE: new constructions after
% ISIT 2014)
figure;
hold all;
axis([alphaBarMSR alphaBarMBR gammaBarMBR gammaBarMSR]);
legendArray = [];
tCAllAux = [];

% NOTE: The Constructions 1 and 2 above are redundant if kIC = k:d is used
% below.
for kIC = k:d,
    kIC
    tradeoffCurve = ChicagoConstruction(n,k,d,kIC);
    tCAux = [alphaBarMBR, gammaBarMBR;
                tradeoffCurve;
             alphaBarMSR, gammaBarMSR];
    tCAllAux = [tCAllAux; tCAux];
    tC = convhull(tCAux(:,1),tCAux(:,2));
    
    % This is to remove the unrequired lines in the plot (on the right side
    % of the space sharing line)
    tC = tC(find(tC==max(tC)):end);
    
    plot(tCAux(tC,1),tCAux(tC,2),'-');
    legendArray = [legendArray;num2str(kIC)];
end
legend(legendArray);

% tCBoth refers to the convex hull obtained by Constructions 1 and 2
tCBothaux = [tC1aux; tC2aux];
tCBoth = convhull(tCBothaux(:,1),tCBothaux(:,2));
plot(tCBothaux(tCBoth,1),tCBothaux(tCBoth,2),'-');

% tCBoth refers to the convex hull obtained by all the constructions in
% this program
tCAll = convhull(tCAllAux(:,1),tCAllAux(:,2));
plot(tCAllAux(tCAll,1),tCAllAux(tCAll,2),'-');

% legend show;