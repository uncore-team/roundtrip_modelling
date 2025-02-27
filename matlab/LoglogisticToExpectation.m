function ex = LoglogisticToExpectation(a,b,c)
% Return the expectation of that loglogistic

    if c >= 1
        ex = NaN;
    else
        ex = a + b * pi * c / sin(pi * c);
    end

end