function [ok, offs, mu, sigma, casesign] = LognormalFit(x)
% Fit through MLE a 3-lognormal to the data in X.
% 
% Taken from our Open-paper (thirdsubmission)
% 
% We estimate the offset and then fit a non-shifted lognormal. A positive 
% random variable X is lognormally distributed if the natural logarithm of 
% X is normally distributed.
%
% pdf(x) = 1/(x * sigma * sqrt(2*pi)) * exp(-( (log(x - offs) - mu)^2/(2*sigma^2) ))
% Expectation: offset + exp(mu + sigma^2 / 2); median: exp(mu); mode: offset + exp(mu - sigma^2); variance: exp(sigma^2) - 1) * exp(2*mu + sigma^2)
%
% Ours, wikipedia and Matlab use the same formulation.
% A lognormal sample cannot have any value <= to its offset. The offset
% estimation procedure deals correctly with such values.
%
% X -> data sample, with min(x) > 0.
%
% OK <- 1 indicating a fit is found, always.
% OFFS <- offset in the data, >= 0 and < min(x).
% MU <- mean of the natural logarithm of the non-shifted data.
% SIGMA <- std of the natural logarithm of the non-shifted data.
% CASESIGN <- for debugging (see estimateOffset() below).

global TOLROUNDTRIPS

    ConstantsInit();

    ok = 0;
    offs = NaN;
    mu = NaN;
    sigma = NaN;

    minx = min(x);
    if minx <= 0
        error('Invalid data sample for lognormal fit');
    end

    [ok,offs,casesign] = estimateOffset(x,0);
    % OK is 3 if MLE has been done; otherwise, a heuristic offset has been
    % provided. In any case, offset < min(x)
    % Given offset, these are the MLE for mu and sigma according to
    % wikipedia (both unbiased):
    %  mu = sum(log(xi - offs)) / n
    %  sigma = sqrt(  sum((log(xi - offs) - mu)^2) / (n - 1) )
    
    % the natural logarithm of the non-shifted data should be normal, so:
    ds = log(x - offs); % reduce to non-shifted normal
    mu = mean(ds); % estimate mu and sigma as the ones of that normal (MLE)
    sigma = std(ds);

    % or, equivalently:
    %[phat,~] = lognfit(x-offs); % phat(1) is mean, phat(2) is std
	%mu = phat(1);
	%sigma = phat(2);

    % figure;
    % [hfreqs,hxs] = hist(ds,50);
    % bar(hxs,hfreqs);
    % hold on;
    % grid;
    % xs = linspace(min(ds) * 0.99,max(ds)*1.01,100000);
    % ys = normpdf(xs,mu,sigma);
    % ys = ys / trapz(xs,ys) * trapz(hxs,hfreqs);
    % plot(xs,ys,'b-');

    ok = 1;

end

function [ok,offset,casesign] = estimateOffset(reg,trace)
% MLE estimation of the lognormal offset of the sample REG; when there is
% no solution in the MLE formula, a heuristic estimation is done.
%
% REG -> sample
% TRACE -> 0 to not trace; 1 to trace console; 2 to trace figs besides
% 
% OK <- 1 if the offset has been estimated heuristically at one of the extremes
%       because the zero-crossing function did not changed its sign at those
%       extremes; 
%       2 if it has been estimated at one of the extremes because
%       the search for it in between has failed; 
%       3 if it has been estimated
%       correctly within the possible range and was unique; 
%       4 in the same
%       situation but has been selected among all the offsets found.
% OFFSET <- offset estimated for the data.
% CASESIGN <- 0 for the case the signs of the zero-crossing function have
%             the same sign in both extremes; 1 if different signs.

global TOLROUNDTRIPS

    ConstantsInit();

    if trace
        fprintf('estimate offset..............\n');
    end
    
    % --------- DATA SANITIZING

    % this reduces the magnitude of the values and thus diminishes the prob. of 
    % numerical inaccuracies without affecting the distribution
    minreg = min(reg);
    % shift the sample by correctionoffset
    if minreg < 0.1
        correctionoffset = 0;
    else
        correctionoffset = min(reg) - 0.1; 
    end
    correctedsample = reg - correctionoffset; 
    if trace
        fprintf('-------- reg[%f,%f] mean(%f) median(%f) mean/median(%f) \n',...
                min(correctedsample),max(correctedsample),mean(correctedsample),median(correctedsample),mean(correctedsample)/median(correctedsample));
    end
    orderedsample = sort(correctedsample);
    
    % ------ ZERO-CROSSING FUNCTION THAT ESTIMATES THE OFFSET WITH MLE

    % This comes from our paper ICCRC2012_AGB_JAFM_ACM, which in turn is
    % based on the modification MMLE-I of Cohen's paper, that is recommended
    % for both small and large samples; we use r = 1.

	r = 1;
    n = length(orderedsample);
    kr = norminv(r/(n+1)); %sqrt(2)*erfinv(2*r/(n+1)-1);
    
    if trace == 2
        gxs = linspace(orderedsample(r)-10,orderedsample(r)-TOLROUNDTRIPS,1000000);
        gys = zeros(1,length(gxs));
        for f = 1:length(gys)
%            gys(f) = MLEfun(orderedsample,gxs(f),n);
            gys(f) = MMLEIfun(orderedsample,gxs(f),r,kr,n);
        end
        hfig = figure;
        plot(gxs,gys,'r.-');
        hold on;
        grid;
        plot(gxs,zeros(1,length(gxs)),'k-','LineWidth',1);
        title(sprintf('Function to find root ofs (3rd central moment of data = %f, x(1) = %f)',moment(correctedsample,3),orderedsample(1)));
        drawnow;
        pause;
        close(hfig);
    end            
    
%    funtofindzeroes = @(g)(MLEfun(orderedsample,g,n)); 
    funtofindzeroes = @(g)(MMLEIfun(orderedsample,g,r,kr,n)); 
    x0 = [TOLROUNDTRIPS, ...
          orderedsample(r) - TOLROUNDTRIPS];
        % Range of search for G
        
    % ---------- EXAMINE ZERO-CROSSING FUNCTION BEFORE USING IT

    f0 = funtofindzeroes(x0(1));
    f1 = funtofindzeroes(x0(end));
    if sign(f0) == sign(f1) % Cannot do the search

        casesign = 0;

        % The crossing can be around data(1), because the data is highly 
        % exponential or around 0 in the case the histogram is very 
        % gaussian.

        offset = heuristicoffsetatextreme(orderedsample,TOLROUNDTRIPS,trace);

        if trace
            fprintf('press a key - signs the same\n');
        end
        if trace == 2
            gxs = linspace(TOLROUNDTRIPS,orderedsample(1),1000000);
            gys = zeros(1,length(gxs));
            for f = 1:length(gys)
                gys(f) = funtofindzeroes(gxs(f));
            end
            hfig = figure;
            subplot(1,2,1);
            plot(gxs,gys,'r.-');        
            hold on;
            grid;
            plot(gxs,zeros(1,length(gxs)),'k-','LineWidth',1);
            plot(orderedsample(1)*[1 1],[min(gys),max(gys)],'k-');
            gxs = linspace(x0(1),x0(end),100000);
            gys = zeros(1,length(gxs));
            for f = 1:length(gys)
                gys(f) = funtofindzeroes(gxs(f));
            end        
            plot(gxs,gys,'mo');        
            title(sprintf('Function to find root ofs (3rd central moment of data = %f, x(1) = %f)',moment(correctedsample,3),orderedsample(1)));
            subplot(1,2,2);
            histogram(orderedsample);
            grid;
            drawnow;
            pause;
            close(hfig);
        end
        ok = 1;
        offset = offset + correctionoffset;

    else % there is a change of sign in function, thus some possible zero crossing

        casesign = 1;

        % Finding zero
        options = optimset('Algorithm', 'levenberg-marquardt'); %'LevenbergMarquardt', 'on');  
        try
            [offsets,fval,exitflag,output] = fzero(funtofindzeroes, x0, options);
            offsets = offsets(isreal(offsets) & (offsets < orderedsample(1)));
            numsols = length(offsets);
            if numsols == 0
               ok = 2;
               offset = heuristicoffset(orderedsample,TOLROUNDTRIPS,trace);
            elseif numsols > 1 % several solutions: cohen chooses the one giving expectation closest to sample mean
               mu = zeros(1,numsols);
               s = zeros(1,numsols);
               for i = 1:lnumsols
                    p = lognfit(orderedsample - offsets(i));
                    mu(i) = p(1);	% mu parameter 2params LN
                    s(i) = p(2);
               end
               esperanzas = offsets + exp(mu+(s.^2)./2);	% expectation of the LLN3 (cohen first page)
               sm = mean(orderedsample); % sample mean
               % busca en esperanza el valor que se parezca mas a mu
               erresps = (esperanzas - sm).^2;
               [~,in] = min(errsps);
               offset = offsets(in);
               ok = 4;
            else
               offset = offsets(1);
               ok = 3;
            end
            if isnan(offset)
                error('nan offset');
            end
            offset = offset + correctionoffset;
        catch errRecord    
            ok = 2;
            offset = heuristicoffset(orderedsample,TOLROUNDTRIPS,trace) + correctionoffset;
            disp(errRecord)
        end
    
        if trace
            fprintf('press a key - signs different ok %d\n',ok);
        end        
        if trace == 2
    
            gxs = linspace(TOLROUNDTRIPS,orderedsample(1),1000000);
            gys = zeros(1,length(gxs));
            for f = 1:length(gys)
                gys(f) = funtofindzeroes(gxs(f));
            end
            hfig = figure;
            plot(gxs,gys,'r.-');
            hold on;
            grid;
            plot(gxs,zeros(1,length(gxs)),'k-','LineWidth',1);
            plot(orderedsample(1)*[1 1],[min(gys),max(gys)],'k-');
            title(sprintf('Function to find root ofs (3rd central moment of data = %f, x(1) = %f)',moment(correctedsample,3),orderedsample(1)));
            subplot(1,2,2);
            histogram(orderedsample);
            grid;
            drawnow;
            pause;
            close(hfig);
        end

    end

    if (offset < 0) || (offset >= minreg)
        error('Invalid offset calculation');
    end
    
end

function offheur = heuristicoffsetatextreme(orderedsample,tole,trace)
% Calculate the offset either as close to 0.0 or close to the minimum of
% the data, by examining the data histogram shape (more gaussian or more
% exponential, respectively).

    [yshist,ehist] = histcounts(orderedsample,'Normalization','pdf');
    xshist = (ehist(2:end)+ehist(1:end-1))/2;
    
    [alpha,beta,expok] = ExponentialFit(orderedsample); % exponential fit
    if expok
        yspdfexp = ExponentialPdf(xshist,alpha,beta);
        mseexp = sum((yspdfexp - yshist).^2);
    else
        mseexp = Inf;
    end

    mu = mean(orderedsample); % normal fit (always exists)
    sigma = std(orderedsample);
    yspdfnorm = normpdf(xshist,mu,sigma);
    msenorm = sum((yspdfnorm - yshist).^2);
    
    if trace
        fprintf('extreme offset; mseexp: %.15f, msenorm: %.15f\n',mseexp,msenorm);
        if mseexp <= msenorm
            fprintf('\tEXP wins\n');
        else
            fprintf('\tNORM wins\n');
        end
    end
    if trace == 2
        fi = figure;
        bar(xshist,yshist);
        hold on;
        grid;
        plot(xshist,yspdfexp,'r.-');
        plot(xshist,yspdfnorm,'b.-');
        fprintf('press a key - ok 0a\n');
        pause;
        close(fi);
    end

    if mseexp <= msenorm % exponential model is better; take its offset
        offheur = orderedsample(1) - tole;
    else % normal model is better; take its offset
        offheur = 0;
    end
end

function ga = MLEfun(orderedsample,gamma,n)
% Original function of Cohen, paper 1951, not used here

    ordsample0 = orderedsample - gamma;
    ordsamplelog = log(ordsample0);
    sumlog = sum(ordsamplelog);
    ga = sum(1./ordsample0) * ...
         ( sumlog - sum(ordsamplelog.^2) + 1/n * sumlog^2 ) - ...
         n * sum(ordsamplelog./ordsample0);

end

function ga = MMLEIfun(orderedsample,gamma,r,kr,n)
% Modified function of Cohen, used here

    lsg = log(orderedsample-gamma);
    sg = sum(lsg);
    ga = log(orderedsample(r)-gamma) - sg/n - kr*sqrt( sum(lsg.*lsg)/n-(sg/n).^2 );
    if ~isreal(ga) % may occur due to the numerical noise
        ga = real(ga);
    end
    
end
