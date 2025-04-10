function A2 = ADstatistic(Z)
% Original, basic (non-modified) AD statistic from p.100 D'Agostino
% Return NaN if Z is invalid for doing the calculations.

    % This statistic measures the squared distance
    % between experimental and theoretical Zs, and, indirectly, between 
    % theoretical and experimental Xs
    Zsorted = sort(Z);
    if (Zsorted(1) < eps) || (1 - Zsorted(end) < eps) % AD test cannot work with samples that produce 0 or 1 Z values
        %warning('Invalid value for Z -either 0 or 1- in AD test; Z1=%.19f, Zn=%.19f',Zsorted(1),Zsorted(end));
        A2 = NaN;
        return;
    end
    n = length(Zsorted);
    is = 1:n;
    sumatoria = sum((is*2-1).*log(Zsorted)+(2*n+1-2*is).*log(1-Zsorted)); % page 101, bottom formula
    A2 = -n - (1/n)* sumatoria;

end