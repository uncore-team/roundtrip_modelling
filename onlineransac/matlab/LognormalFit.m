function [ok, offs, mu, sigma] = LognormalFit(x)
% Fit through MLE a 3-lognormal to the data in X.
% Taken from our Open-paper (thirdsubmission)

global TOLROUNDTRIPS

    ConstantsInit();

    minsepoffset = TOLROUNDTRIPS; % x-offset must be >0
    minx = min(x);
    [ok,offs] = estimateOffset(x,0);
    if ~ok
        offs = NaN;
        mu = NaN;
        sigma = NaN;
        return;
    elseif ok == 2
        offs = minx - minsepoffset;
    else
        if (minx - offs < minsepoffset)
            offs = minx - minsepoffset; 
        end
    end
    if offs > min(x)
        error('invalid lognormal fit');
    end
    ds = log(x-offs);
    mu = mean(ds);
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

end

function [ok,offset] = estimateOffset(regimen,trace)
% Calcula el offset segun el MLE de la Lognormal
% Si no se encuentra cambio de signo al minimizar la funcion para calcular
% el offset, se devuelve ok == 0; si se encuentra pero luego no es capaz de
% encontrar el offset mediante optimizacion, se devuelve ok == 2.
% Esto viene del paper de la lognormal nuestro: ICCRC2012_AGB_JAFM_ACM, que está
% basado en Cohen, en la modificación MMLE-I, que es la que al final recomiendan
% tanto para muestras pequeñas como grandes, usando además r = 1. Cohen usa
% la misma formulacion de la LLN que en la wikipedia y que matlab: mu, sigma (y gamma ==
% offset).
global TOLROUNDTRIPS

    ConstantsInit();

    fprintf('..............\n');
    
    corrsample = min(regimen) - 0.1;
    regimen = regimen - corrsample; % to reduce the magnitude of the values and thus diminish the prob. of numerical noise
    
	r = 1;
    orderedsample = sort(regimen);
    n = length(orderedsample);
    kr = norminv(r/(n+1)); %sqrt(2)*erfinv(2*r/(n+1)-1);
    
    if trace == 1
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
        title(sprintf('Function to find root ofs (3rd central moment of data = %f, x(1) = %f)',moment(regimen,3),orderedsample(1)));
        drawnow;
        pause;
        close(hfig);
    end            
    
%    f = @(g)(MLEfun(orderedsample,g,n)); 
    f = @(g)(MMLEIfun(orderedsample,g,r,kr,n)); 
    x0 = [TOLROUNDTRIPS, ...
          orderedsample(r) - TOLROUNDTRIPS];
        % Range of search for G
        
    f0 = f(x0(1));
    f1 = f(x0(end));

    if sign(f0) == sign(f1) 
   
        % If both are negative, indicating that the crossing should be
        % around gamma == 0, is because the data is highly normal
        
        % If both are positive, indicating that the crossing should be
        % around data(1), is because the data is highly exponentiel
        
        %warning('Cannot do the search unless sign(F(minbound)) ~= sign(F(maxbound))');
        
        fprintf('--------\n');
        
        ok = 0;
        offset = NaN;

        if trace == 2
            gxs = linspace(TOLROUNDTRIPS,orderedsample(1),1000000);
            gys = zeros(1,length(gxs));
            for f = 1:length(gys)
    %            gys(f) = MLEfun(orderedsample,gxs(f),n);
                gys(f) = MMLEIfun(orderedsample,gxs(f),r,kr,n);
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
    %            gys(f) = MLEfun(orderedsample,gxs(f),n);
                gys(f) = MMLEIfun(orderedsample,gxs(f),r,kr,n);
            end        
            plot(gxs,gys,'mo');        
            title(sprintf('Function to find root ofs (3rd central moment of data = %f, x(1) = %f)',moment(regimen,3),orderedsample(1)));
            subplot(1,2,2);
            histogram(regimen,100);
            grid;
            drawnow;
            ok
            fprintf('press a key - ok 0\n');
            pause;
            close(hfig);
        end

        return;
    end    
    
    ok = 1;

    % Finding zero
    options = optimset('Algorithm', 'levenberg-marquardt'); %'LevenbergMarquardt', 'on');  
    try
        [offsets,fval,exitflag,output] = fzero(f, x0, options);
        offsets = offsets(isreal(offsets) & (offsets < orderedsample(1)));
        numsols = length(offsets);
        if numsols == 0
           ok = 2;
           offset = NaN;
        elseif numsols > 1 % several solutions: cohen chooses the one giving expectation closest to sample mean
           mu = zeros(1,numsols);
           s = zeros(1,numsols);
           for i = 1:lnumsols
                p = lognfit(orderedsample-offsets(i));
                mu(i) = p(1);	% mu parameter 2params LN
                s(i) = p(2);
           end
           esperanzas = offsets + exp(mu+(s.^2)./2);	% expectation of the LLN3 (cohen first page)
           sm = mean(orderedsample); % sample mean
           % busca en esperanza el valor que se parezca mas a mu
           erresps = (esperanzas - sm).^2;
           [~,in] = min(errsps);
           offset = offsets(in);
        else
           offset = offsets(1);
        end
        if ~isnan(offset)
            offset = offset + corrsample;
        end
    catch errRecord    
        ok = 2;
        offset = NaN;
        disp(errRecord)
    end

    if (ok ~= 1) && (trace == 2)

        fprintf('=========\n');
        
        gxs = linspace(TOLROUNDTRIPS,orderedsample(1),1000000);
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
        plot(orderedsample(1)*[1 1],[min(gys),max(gys)],'k-');
        title(sprintf('Function to find root ofs (3rd central moment of data = %f, x(1) = %f)',moment(regimen,3),orderedsample(1)));
        subplot(1,2,2);
        histogram(regimen,100);
        grid;
        drawnow;
        ok
        fprintf('press a key - ok 2\n');
        pause;
        close(hfig);
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
