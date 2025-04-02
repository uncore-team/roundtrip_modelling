function [ok,ps] = ReqCheckInSample(r,S)
% Check whether the requirement R is satisfied by sample S.
%
% R -> Requirement (see ReqCreate).
% S -> set of roundtrip times.
%
% OK <- 1 if satisfied, 0 if not.
% PS <- prob of satisfying the req in the sample.

    ok = 0;
    ps = 0;
    if isempty(S)
        return;
    end

    % Estimate the Bernoulli distribution parameter, which equals the prob PS.
    % This Bernoulli distribution is a discrete one in the support {"satisfies
    % the requirement", "does not satisfy the requirement"}.
    nyes = length(find((S >= r.t0) & (S <= r.t1)));
    ntot = length(S);
    ps = nyes / ntot;
    if ps >= r.p
        ok = 1;
    end

end