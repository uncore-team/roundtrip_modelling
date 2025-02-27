function r = LoglogisticRnd(a, b, c, m, n)
% generate m*n data from loglogistic (a,b,c)

	alpha=b;
	betha=1/c;
	off=a;

	mu=log(alpha);
	sigma=1/betha;

	p=rand(m,n);
	r=log(p./(1-p)).*sigma+mu;

	r=exp(r)+off;

end