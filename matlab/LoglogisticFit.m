function [a, b, c, exitflag] = LoglogisticFit(x)
% Fit through MLE a 3-loglogistic to the data in X
%
% EXITFLAG <- the exitflag of the last optimization call if it is incorrect
% that is, <0 (then, parameters a,b,c are filled with initial guesses). If
% exitflag is >= 0, the model can be considered acceptable, even if it is 0
% (too many iterations)

	minlen=10;
	tole=1e-9;
	optimalg='trust-region-reflective';%'trust-region-reflective'; % 'trust-region-reflective' 'levenberg-marquardt'
	minc=0.05; % 0.05 has been used in the c++ hub because it produces almost gaussians if c below that
               % and because a too small value of minc produces overflows
               % in the calculations of the functions to be optimized
	maxc=1/2-eps; % c<1 for having expectation; c<1/2 for having variance
    debug=0;

    lx=length(x);
	if (lx<minlen)
		error('Cannot fit anything with less than %d values',minlen);
	end
    minx=min(x);
	if (minx<=0)
		error('Cannot fit a loglogistic with minx<=0');
	end
    maxx=max(x);
    mux=mean(x);
    
    % initial guess for a
    ahat=minx-(maxx-minx)/1e4; % just a crude estimate of where the a could be, to start with
    if (ahat<0)
        ahat=0;
    end
    
    % initial guess for b
    bhat = median(x-ahat);
    if (bhat < 1e-6)
        bhat = 1e-6; % to avoid some numerical errors in the fitting pass
    end
%     phat=mle(log(x-ahat),'distribution','logistic'); % have not obtained better estimates in several tests
%     bhat=exp(phat(1)); % in principle, this is a better estimation of b than the sample median, at the cost of another optimization by mle()

	% initial guess for c, given those initial a and b
    options=optimset('Display', 'off', 'Algorithm', optimalg,'TolFun', tole, 'TolX', tole); 
    boundmin=[minc];
    boundmax=[maxc]; 
	c0=minc;
    [chat,RESNORM,RESIDUAL,exitflag,OUTPUT] = lsqnonlin(@(chat) ahat + bhat*beta(1+chat, 1-chat) - mux,c0,...
                                                        boundmin,boundmax,options);
    a=ahat;
    b=bhat;
    c=chat;    
    if (exitflag<0)
        return;
    end

    % do the fitting of the entire loglogistic from that seed
    x0 = [ahat; bhat; chat]; % start at initial guess
    %fprintf('Fitting ll3 with ahat = %.24f, bhat = %.24f, chat = %.24f...\n',x0(1),x0(2),x0(3));
    options=optimset('Display', 'off', 'Algorithm', optimalg, 'Jacobian', 'on', 'TolFun', tole, 'TolX', tole);%, 'Display','iter');
    lb = [eps; eps; minc]; 
    ub = [minx-eps; inf; maxc];        
    [xx,RESNORM,RESIDUAL,exitflag,OUTPUT] = lsqnonlin(@(x0) Loglogistic_fittingfunctions(x0,x,lx), x0,lb,ub,options);
    if (exitflag<0)
        return;
    end
    a = xx(1);
    b = xx(2);
    c = xx(3);
    
    % some post-security checks that display warnings (should not happen)
	smallval=1e-9;
    if (a<=0.0) 
        if (debug==1)
            fprintf('A!');
        end
        a=smallval;
    end
    if (a>=minx) 
        if (debug==1)
            fprintf('M!');
        end
        a=minx-smallval;
    end
    if (b<=0.0) 
        if (debug==1)
            fprintf('B!');
        end
        b=smallval;
    end
    if (c<=0.0) 
        if (debug==1)
            fprintf('C!');
        end
        c=smallval;
    end

end