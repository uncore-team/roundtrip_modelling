function ex = ExponentialToExpectation(alpha,beta)
% See the distrib. parameters in ExponentialFit.

	ExponentialCheckParms(alpha,beta);
    ex = alpha + beta;
    
end
