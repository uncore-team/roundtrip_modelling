function [b,c,exitflag]=LoglogisticFromEV(E,V)
% Given in E and V the expectation and variance, respectively, of a
% loglogistic distribution with 0 offset, return in B and C the values of
% its parameters that achieve such E and V
% https://en.wikipedia.org/wiki/Log-logistic_distribution
% E[] = b*(c*pi)/sin(c*pi)
% V[] = b^2 * ( (2*pi*c)/sin(2*pi*c) - (pi*c)^2 / (sin(pi*c)^2) )
% c>0, c<1/2, b>0

    % the following algorithm is in
    % calculo_de_b_y_c_de_LL3_a_partir_de_E_y_V.png
    
    constantopt = log(V/(E^2)+1); % log of the function for avoiding large values
    [pic,RESNORM,RESIDUAL,exitflag,OUTPUT,lambda,jacobian] = lsqnonlin(@(pic) log(tan(pic))-log(pic)-constantopt,...
                                                        0.25*pi,...
                                                        [eps],[0.5*pi-eps],...
                                                        optimset('Display', 'off', 'Algorithm','trust-region-reflective','TolFun',1e-9,'TolX',1e-9));
    c=pic/pi;
    b=E*sin(pic)/pic;
    
    [e,v]=LoglogisticToEV(0,b,c);
    dife = abs(e-E);
    difv = abs(v-V);
    if (dife/E>0.1)||(difv/V>0.1)
%        warning('Poor deduction of b,c for loglogistic. E=%.24f, V=%.24f\n\tDeduction: b=%.15f, c=%.15f\n\twith those b,c: E=%.24f, V=%.24f\n\tdifference in E=%.2f%%, difference in V=%.2f%%',E,V,b,c,e,v,dife/E*100,difv/V*100);
        exitflag=[dife/E difv/V];
    else
        exitflag=1; % found a solution
    end
    
    return;

    
    
    % --- older algorithm: it achieves the same but with more computational
    % cost (search in two-variable space) and obtaining worse results for
    % both b,c in the case that c is close to 0.5 (the previous algorithm
    % obtains a good result for c in that case).
    
    lb = [eps; eps]; % lower bounds for b and c
    ub = [Inf; 0.5-eps]; % upper bounds for b and c
    x0 = [0.05; 0.25]; % initial guesses        
    
    x0=[b;c];
    options=optimset('Display', 'off', ...
                     'Algorithm', 'trust-region-reflective',...%'levenberg-marquardt', ...
                     'Jacobian', 'off', ...
                     'TolFun', 1e-24, 'TolX', 1e-24); % needed a very small tolerance to obtain good approximations
    [xx,RESNORM,RESIDUAL,exitflag,OUTPUT] = lsqnonlin(@(x0) LLfun(x0,E,V),x0,lb,ub,options);
    b=xx(1);
    c=xx(2);

end

function F = LLfun(x,E,V)
% Implement the two functions to optimize

    F=zeros(2,1);
    b=x(1);
    c=x(2);
    [ee,vv]=LoglogisticToEV(0,b,c);
    F(1)=ee-E;
    F(2)=vv-V;
%     cpi=c*pi;
%     scpi=sin(cpi);
%     F(1)=b*cpi/scpi-E;
%     F(2)=b^2*(2*cpi/sin(2*cpi)-(cpi^2)/(scpi^2))-V;

end
