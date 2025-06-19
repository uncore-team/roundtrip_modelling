function [reject,stat,thresh] = ExponentialGof(x,alpha,beta,modelnotfromdata)
% Based on D'Agostino p. 141: both parameters unknown. 
% See ExponentialFit()

    ExponentialCheckParms(alpha,beta);

    n = numel(x);
    x = reshape(x,1,n); % force X is a row vector

	% ---- transform sample to a theoretically uniform one (not sorted)

    w = (x - alpha) / beta;
    Z = 1 - exp(-w);

    % ---- calculate EDF statistic
    [kindstat,stat] = StatisticGof(Z);
    if strcmp(kindstat,'A2')

        if modelnotfromdata
            % no correction of stat in this case (see table 4.2 D'Agostino)
    	    thresh = 2.492; % known parameters (n>=5)
                            % table 4.2 D'Agostino; confirmed with MonteCarlo
                            % experiments (test_tabulategofthrs.m)
        else % model comes from data
            % correction: A2 for case 3 (both parameters were deduced from the same sample). 
            % do the following since parameters come from sample (D'Agostino table 4.14)
            stat = stat * (1 + 5.4/n - 11/n^2); 
            % the following is D'Agostino and it does work even when it 
            % assumes only one parameter unknown while we have 2; we have
            % conducted new MonteCarlo experiments to deduce the threshold in
            % this case in test_tabulategofthrs.m and got the same threshold
            % for sample sizes ranging from 20 to 10000.            
            thresh = 1.321; 
        end

    elseif strcmp(kindstat,'W2')

        if modelnotfromdata
            % correction if parms are true (not from sample); table 4.2
            % D'Agostino
            stat = (stat - 0.4/n + 0.6/n^2) * (1 + 1/n);
            thresh = 0.461; % we have confirmed this value with MonteCarlo (test_tabulategofthrs.m)
        else % model comes from data
            % correction if parms come from sample; table 4.14 D'Agostino
            stat = stat * (1 + 2.8/n - 3/n^2);
            thresh = 0.222; % we have confirmed this value with MonteCarlo (test_tabulategofthrs.m)
        end

    else
        error('Unknown gof statistic');
    end

    if (stat > thresh) % equivalently, the p-value is smaller than the significant level
        reject=1; % reject
    else
        reject=0; % cannot reject the null hypothesis of the data coming from an exponential
    end 

end
