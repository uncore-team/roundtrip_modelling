function r = LoglogisticRnd(a, b, c, m, n)
% Generate m*n data from loglogistic (a,b,c)
%
% The resulting sample may have values equal to the offset (a) due to a
% reason only (since the LL cannot generate naturally such cases): when 
% adding a large offset to a small value; this is due to numerical 
% imprecissions.
% This can happens also in a real scenario when we are measuring delays:
% some of them may come from a LL (continuous support, i.e., infinite
% precision) but become truncated in their precision.

	alpha=b;
	betha=1/c;
	off=a;

	mu=log(alpha);
	sigma=1/betha;

	p=rand(m,n);
	r=log(p./(1-p)).*sigma+mu;

	r=exp(r)+off;

end