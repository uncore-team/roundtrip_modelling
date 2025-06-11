function ys = ExponentialPdf(xs,alpha,beta)
% See the distrib. parameters in ExponentialFit.
    
	ExponentialCheckParms(alpha,beta);
    ys = exp(-(xs - alpha)/beta)/beta;

end
