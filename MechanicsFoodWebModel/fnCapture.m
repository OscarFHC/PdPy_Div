%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matlab Code supplementing the paper
% The mechanics of predator-prey interactions: first principles of physics predict predator-prey size ratios
% by Portalier, Fussmann, Loreau, Cherif
%
% Functional Ecology

% November 2018
%
% Matlab version: R2018b
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%

%%%% READ ME %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% The following code defines the fnCapture function used within the Mainfile
% to find a root to predator and prey distances covered during capture sequence 
%
%%% The functions takes two arguments: a single value (x) and a cell array (p)  
% 
% x is the value at which the function is evaluated
%
% p is a cell array (2 cells):
% p{1}: the polynom values for distance covered by the predator
% p{2}: the polynom values for distance covered by the prey
%
%%% The function evaluates the predator and prey distances using polynomial coefficients 
%
%%% The function returns the absolute difference between the two distances (y)
% This distance is evaluated by fmincon function (Mainfile)
%  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Capture sequence
function [y] = fnCapture(x,p)
    P1=p{1};
    P2=p{2};
    Predict1=polyval(P1,x);
    Predict2=polyval(P2,x);
    y=abs(Predict2-Predict1);
end