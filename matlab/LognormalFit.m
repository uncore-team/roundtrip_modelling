function [ok, offs, mu, sigma] = LognormalFit(x)
% Fit through MLE a 3-lognormal to the data in X.
% Taken from our Open-paper (thirdsubmission)

    minsepoffset = 1e-9; % x-offset must be >0
    minx = min(x);
    [ok,offs] = estimateOffset(x,0); 
    if ~ok
        offs = NaN;
        mu = NaN;
        sigma = NaN;
        return;
    end
    if (minx - offs < minsepoffset)
        offs = minx - minsepoffset; 
    end
    if offs > min(x)
        error('invalid lognormal fit');
    end
    xs = log(x-offs);
    mu = mean(xs);
    sigma = std(xs);
    % or, equivalently:
    %[phat,~] = lognfit(x-offs); % phat(1) is mean, phat(2) is std
	%mu = phat(1);
	%sigma = phat(2);

end

function [ok,offset] = estimateOffset(regimen,trace)
% Calcula el offset segun el MLE de la Lognormal
% Si no se encuentra cambio de signo al minimizar la funcion para calcular
% el offset, se elegira el minimo como offset
% Esto viene del paper de la lognormal nuestro: ICCRC2012_AGB_JAFM_ACM, que está
% basado en Cohen, en la modificación MMLE-I, que es la que al final recomiendan
% tanto para muestras pequeñas como grandes, usando además r = 1. Cohen usa
% la misma formulacion de la LLN que en la wikipedia y que matlab: mu, sigma (y gamma ==
% offset).

	r = 1;
    orderedsample = sort(regimen);
    n = length(orderedsample);
    kr = norminv(r/(n+1)); %sqrt(2)*erfinv(2*r/(n+1)-1);
    
    if trace
        gxs = linspace(orderedsample(r)-10,orderedsample(r)-1e-9,1000000);
        gys = zeros(1,length(gxs));
        for f = 1:length(gys)
%            gys(f) = MLEfun(orderedsample,gxs(f),n);
            gys(f) = MMLEIfun(orderedsample,gxs(f),r,kr,n);
        end
        figure;
        plot(gxs,gys,'r-','LineWidth',2);
        hold on;
        grid;
        plot(gxs,zeros(1,length(gxs)),'k-','LineWidth',1);
        title(sprintf('Function to find root ofs (3rd central moment of data = %f, x(1) = %f)',moment(regimen,3),orderedsample(1)));
    end

%    f = @(g)(MLEfun(orderedsample,g,n)); 
    f = @(g)(MMLEIfun(orderedsample,g,r,kr,n)); 
    x0 = [0 orderedsample(r)-1e-9]; % range of search for G
    if sign(MMLEIfun(orderedsample,x0(1),r,kr,n)) == sign(MMLEIfun(orderedsample,x0(end),r,kr,n))
        % in many cases, this occurs: the entire gamma curve is below 0
        %warning('Cannot do the search unless sign(F(minbound)) ~= sign(F(maxbound))');
        ok = 0;
        offset = NaN;
        return;
    end

    ok = 1;

    % LevenbergMarquardt
    options = optimset('Algorithm', 'levenberg-marquardt'); %'LevenbergMarquardt', 'on');  
    try
        offsets = fzero(f, x0, options);
        offsets = offsets(offsets < orderedsample(1));
        numsols = length(offsets);
        if numsols == 0
            offset = orderedsample(1) - 1e-9;
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
    catch errRecord    
        ok = 0;
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
    fprintf('%f\n',ga);
end
