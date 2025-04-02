function y = LoglogisticLikelihood(x, a,b,c)
% Return the value of the pdf only for a given X

    if (x<=a)
        y=0;
    else
        y = ((x-a)/b)^(-1/c) / ( c*(x-a) * (1 + ((x-a)/b)^(-1/c))^2 );
    end
    
end
