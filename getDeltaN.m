function [delN1, N1, N2, M] = getDeltaN(N1,arrival_day,p)

% -- [delN1, N1, N2, M] = getDeltaN(N1,arrival_day,p)
%
% This function serves two purposes. First, it can be used
% to calculate the change in population size in one year
% given the parameters supplied. Therefore it can be used
% with a numerical solver to find the population steady
% state. Second, it can return the daily dynamics of the
% focal and alternative-prey species during the predation
% and territory-acquisition period. This is helpful for
% calculating invasion fitness of mutant strategies. The
% function is used in both ways in calc_fit.m.
%
%
% INPUTS:
%
% N1: The population size of the focal species at the start
% of the year.
%
% arrival_day: The prevailing arrival-day strategy.
%
% p: A dictionary of parameter values. See params.m for an example.
%
% OUTPUTS:
%
% delN1: The change in population size over 1 year.
%
% N1: A vector of the focal species' population size over
% the predation and territory-acquisition period.
% 
% N2: A vector of the alternative-prey species' population size 
% over the predation and territory-acquisition period.
% 
% M: Focal population size at the end of the predation and
% territory-acquisition period.
%
%
% See also: calc_fit.m

tmax = p.tmax; P0 = p.P0; Vg = p.Vg; Vb = p.Vb; Rg = p.Rg; Rb = p.Rb; gamma = p.gamma; N2 = p.N2; X = p.X; a = p.a; b = p.b;


for t = 1:tmax;
    if t < arrival_day
        N1(t+1) = N1(t) - P0*N1(t);
        N2(t+1) = max(0,N2(t) - (1-a)*N2(t)*X/(1+(1-a)*b*N2(t)));
    else
        if arrival_day == t
            N1_ad = N1(t);
            P1_ad = a*N1(t)*X/(1+a*b*N1(t)+(1-a)*b*N2(t));
        end
        N1(t+1) = max(0,N1(t) - a*N1(t)*X/(1+a*b*N1(t)+(1-a)*b*N2(t)));
        N2(t+1) = max(0,N2(t) - (1-a)*N2(t)*X/(1+a*b*N1(t)+(1-a)*b*N2(t)));
    end
end
M = N1(tmax+1);

% Parameter-values have been confined to the saturated case
Mg = Vg; 
Mb = Vb;

% Population dynamics equation
N_x = (1-gamma)*(M + Mg*Rg + Mb*Rb);

% Derivative for solver
delN1 = N_x - N1(1);

