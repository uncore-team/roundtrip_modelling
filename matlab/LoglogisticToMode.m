function m=LoglogisticToMode(a,b,c)
% return the mode of that 3-loglogistic

    if (c>=1)
        m = NaN;
        warning('Loglogistic does not have mode if c>=1');
        return;
    end
    m=b*((1-c)/(1+c))^c + a;

end