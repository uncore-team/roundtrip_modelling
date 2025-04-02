function h = drawHisto(hdata,tit,xlab)
% Draw an histogram of HDATA with the given TIT and XLABel, showing it as
% relative frequencies.

    validdata = hdata((~isnan(hdata)) & (~isinf(hdata)));
    if ~isempty(validdata)
        h = histogram(validdata,'Normalization','pdf');
        hold on;
        text(mean(h.BinLimits),max(h.Values)/2,...
             sprintf('%s: mu %.3f, s: %.3f',tit,mean(validdata),std(validdata)));
        xlabel(xlab);
        ylabel('relative freq.');
    else
        h = NaN;
    end

end