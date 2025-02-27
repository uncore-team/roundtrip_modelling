function v = LoglogisticToVariance(a,b,c)
% Return the variance of that loglogistic.

    if c >= 0.5
        v = NaN;
    else
        pic = pi * c;
        v = b * (2 * pic / sin(2 * pic) - (pic/sin(pic))^2);
    end

end