
x=load('rtts.txt');
lx=length(x);

tole=1e-9;
optimalg='trust-region-reflective'; % 'trust-region-reflective' 'levenberg-marquardt'

minc=0.05; % 0.05 has been used in the c++ hub because it produces almost gaussians if c below that
           % and because a too small value of minc produces overflows
           % in the calculations of the functions to be optimized
maxc=1/2-eps; % c<1 for having expectation; c<1/2 for having variance

minx=min(x);

lb = [eps; eps; minc]; 
ub = [minx-eps; inf; maxc];

% do the fitting of the entire loglogistic from that seed
x0 = [  1000.0005049499999; 0.1518950500000642; 0.41512074819075279 ]; % start at initial guess

options=optimset('Display', 'off', ...
                'Algorithm', optimalg, ...
                'Jacobian', 'on', ...
                'TolFun', tole, ...
                'TolX', tole ); % 'Display','iter'

[x_,~,~,exitflag,~] = lsqnonlin(@(x0) Loglogistic_fittingfunctions(x0,x,lx), x0,lb,ub,options);

if (exitflag<0)
    error('Optimization has not converge.');
else
	x_
end
