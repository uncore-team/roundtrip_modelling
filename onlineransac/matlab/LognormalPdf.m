function y = LognormalPdf(x, offset,mu,sigma)
% Return the pdf for that support and parameters

    y = lognpdf(x-offset,mu,sigma);
% equivalently:     y= 1./((x-offset)*sigma*sqrt(2*pi)) .* exp ( (- (log(x-offset)-mu).^2)./(2*sigma^2) );
    y(x-offset<=0) = 0;
    
end
