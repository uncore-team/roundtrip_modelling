function ExponentialCheckParms(alpha,beta)
% Gives an error if ALPHA,BETA are not valid.
% See the distrib. parameters in ExponentialFit.

	if (alpha < 0) || (beta <= 0)
        error('Invalid parameters for exponential distr.');
    end

end
