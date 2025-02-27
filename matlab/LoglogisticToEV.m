function [E,V]=LoglogisticToEV(a,b,c)
% Return the expectation and variance for the 3-logl, which will include the
% contribution of the offset.

    E = NaN;
    V = NaN;
    
    if (c>=1)
        warning('Loglogistic does not have expectation if c>=1');
    end
    if (c>=0.5)
        warning('Loglogistic does not have variance if c>=1/2');
    end

	gammac=gamma(c);
	gamma1_c=gamma(1-c);
	gamma2c=gamma(2*c);
	gamma1_2c=gamma(1-2*c);
    if (c < 1)
    	E = a+b*c*gammac*gamma1_c; % formulas from the hidrology paper
    end
    if (c < 0.5)
    	V = b^2*(2*c*gamma2c*gamma1_2c-(c*gammac*gamma1_c)^2);
    end
	
	% formulas for wikipedia, less accurate when c is close to 0.5, since in that
	% case 2*c is close to an integer, a situation where the Euler's reflection
	% formula used for transforming gamma() to sin() does not apply.
%	pic = pi*c;
%	pic2 = 2*pic;
%	sinpic = sin(pic);
%	sinpic2 = sin(pic2);
%	E = a+b*pic/sinpic;
%	V = b^2*(pic2/sinpic2-(pic/sinpic)^2);

end
