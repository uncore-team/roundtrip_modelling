function x = ExponentialRnd(alpha,beta,n,m)
% See the distrib. parameters in ExponentialFit.
% The resulting sample may have values equal to the offset (alpha) due to 2
% reasons: they are naturally equal to it; they become equal to it when
% adding a large offset to a small value, due to numerical imprecissions.
% The latter can happens also in a real scenario when we are measuring delays:
% some of them may come from an EXP (continuous support, i.e., infinite
% precision) but become truncated in their precision.

	ExponentialCheckParms(alpha,beta);

    x = exprnd(beta,n,m) + alpha; % matlab requires the mean (without offset) as
    							  % the first parameter of exprnd().

end
