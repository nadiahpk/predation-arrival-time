function [fit_e,fit_l] = calc_fit(ad,p)

% -- [fit_e,fit_l] = calc_fit(ad,p)
%
% The purpose of this function is to calculate the fitness of mutants
% invading the system with arrival days one day earlier and one day
% later than the prevailing arrival-day strategy.
%
% If either fitness is greater than 1, then the population can be
% invaded by mutants pursuing the mutant arrival-day strategy.
%
%
% INPUTS:
%
% ad: The prevailing arrival-day strategy
%
% p: A dictionary of parameter values. See params.m for an example.
%
% OUTPUTS:
%
% fit_e: The fitness of a mutant arriving one day earlier than the
% prevailing arrival-day strategy.
% 
% fit_l: The fitness of a mutant arriving one day later than the
% prevailing arrival-day strategy.
% 
%
% See also: getDeltaN.m

% Read in parameter values
tmax = p.tmax; P0 = p.P0; Vg = p.Vg; Vb = p.Vb; Rg = p.Rg; Rb = p.Rb; gamma = p.gamma; N2 = p.N2; X = p.X; a = p.a; b = p.b;

% Find equilibrium for resident population, N1
deltaN = @(N0) getDeltaN(N0,ad,p);
[N1,fval] = fzero(deltaN,1000); 

% Also get the N1 on each day, N2 on each day, and M
[deltaN, N1, N2, M] = getDeltaN(N1,ad,p); 

if N1(ad) < Vb+Vg || M < Vg+Vb
    % If we're not in the saturated parameter space
    % with a viable population
    fit_e = NaN;
    fit_l = NaN;
else
    % = Mutant Survival =

    % Survival for mutant that arrives one day earlier
    if ad == 1
        Ps_e = NaN; % Cannot arrive one day earlier
    else
        Ps_e = (M/N1(1)) * max(0, ( 1 - (a*X)/( 1 + (1-a)*b*N2(ad-1) ) ) / (1-P0));
    end

    % Survival for mutant that arrives one day later
    if ad == tmax
        Ps_l = NaN; % Cannot arrive one day later
    else
        Ps_l = (M/N1(1)) * max(0, (1-P0) / ( 1 - (a*X)/( 1 + a*b*N1(ad) + (1-a)*b*N2(ad) ) ));
    end

    % = Mutant Territory Acquisition =

    % Territory acquisition for mutant that arrives one day earlier
    Pg_e = 1; 
    Pb_e = 0;

    % Territory acquisition for mutant that arrives one day later
    Pg_l = Vg/M - Vg/N1(ad+1);
    Pb_l = Vb/M - Vb/N1(ad+1);

    % Fitness of mutants, one day earlier and one day later
    fit_e = (1-gamma)*Ps_e*(1 + Pg_e*Rg + Pb_e*Rb);
    fit_l = (1-gamma)*Ps_l*(1 + Pg_l*Rg + Pb_l*Rb);
end
