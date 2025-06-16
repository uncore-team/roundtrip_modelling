function y = LoglogisticPdf(x, a,b,c)
% Return the pdf for that support and parameters

    y = ((x-a)/b).^(-1/c) ./ ( c*(x-a) .* (1 + ((x-a)./b).^(-1/c)).^2 ); 
    y(x <= a) = 0;
    
end
