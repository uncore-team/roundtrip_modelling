function r = LognormalRnd(offs, mu, sigma, m, n)
% See lognormal parms in LognormalFit.m
%
% The resulting sample may have values equal to the offset (offs) due to a
% reason only (since the LN cannot generate naturally such cases): when 
% adding a large offset to a small value; this is due to numerical 
% imprecissions.
% This can happens also in a real scenario when we are measuring delays:
% some of them may come from a LN (continuous support, i.e., infinite
% precision) but become truncated in their precision.

global TOLROUNDTRIPS

    ConstantsInit();

    LognormalCheckParms(offs,mu,sigma);

    % logrnd() can produce datum too close to 0 that, if added to a large
    % offset, produce data equal to the offset.
    r = lognrnd(mu,sigma,m,n);

    % no solution is done here since they either distort the distribution
    % or take too long to work
    % alreadygenerated = 0;
    % while ~alreadygenerated
    % 
    %     % this solution does not distort the distribution but takes longer,
    %     % specially if the sample is large
    %     if isempty(find(r < TOLROUNDTRIPS, 1))
    %         alreadygenerated = 1;
    %     else
    %         r = lognrnd(mu,sigma,m,n);
    %     end
    % 
    %     % % this solution distorts the distribution but it is very fast
    %     % r(r < TOLROUNDTRIPS) = TOLROUNDTRIPS;
    %     % alreadygenerated = 1;
    % 
    %     % % this solution distorts the distribution (since prevents the
    %     % % appearance of otherwise valid -though small- values).
    %     % [rowsmall,colsmall] = find(r < TOLROUNDTRIPS);
    %     % if isempty(rowsmall) % then, colsmall is empty too
    %     %     alreadygenerated = 1;
    %     % else % iid regeneration of too small values
    %     %     nr = length(rowsmall);
    %     %     nc = length(colsmall);
    %     %     r2 = lognrnd(mu,sigma,nr,nc);
    %     %     r(rowsmall,colsmall) = r2;
    %     % end
    % end
    
    r = r + offs;

end
