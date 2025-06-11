function x = ExponentialRnd(alpha,beta,n,m)
% See the distrib. parameters in ExponentialFit.

	ExponentialCheckParms(alpha,beta);
    x = exprnd(beta,n,m) + alpha; % matlab requires the mean (without offset) as
    							  % the first parameter of exprnd().

end
