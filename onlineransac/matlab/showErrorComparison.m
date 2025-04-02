function showErrorComparison(basertts,indserrs,errsignal1,errsignal2,tittot,tit1,tit2,units)
% Draw a figure with both error signals being compared and also print some
% stuff.

    n = length(errsignal1);
    if n ~= length(errsignal2)
        error('Invalid comparison signals');
    end
    
    mrtts = mean(basertts);

    figure;
    subplot(3,2,1);
    plot((basertts-min(basertts))/(max(basertts)-min(basertts))*(max(errsignal1)-min(errsignal1))+min(errsignal1),'-');
    hold on;
    grid;
    plot(indserrs,errsignal1,'.-');
    plot(indserrs,errsignal2,'.-');
    legend(sprintf('base %s',units),tit1,tit2);
    title(tittot);
    subplot(3,2,2);
    plot(indserrs,errsignal1/mrtts*100,'.-');
    hold on;
    grid;
    plot(indserrs,errsignal2/mrtts*100,'.-');
    title(sprintf('% w.r.t. mean(%s)',units));

    subplot(3,2,3);
    drawHisto(errsignal1,tit1,'');
    hold on;
    drawHisto(errsignal2,tit2,'');
    subplot(3,2,4);
    drawHisto(errsignal1/mrtts*100,sprintf('%s %%',tit1),'');
    hold on;
    drawHisto(errsignal2/mrtts*100,sprintf('%s %%',tit2),'');

    subplot(3,2,5);
    drawHisto(abs(errsignal1),sprintf('abs(%s)',tit1),'');
    hold on;
    drawHisto(abs(errsignal2),sprintf('abs(%s)',tit2),'');
    subplot(3,2,6);
    drawHisto(abs(errsignal1/mrtts*100),sprintf('abs(%s)/mean(%s) %%',tit1,units),'');
    hold on;
    drawHisto(abs(errsignal2/mrtts*100),sprintf('abs(%s)/mean(%s) %%',tit2,units),'');

    nbett = length(find(abs(errsignal1) <= abs(errsignal2)));
    fprintf('[%s] %s is better than %s in %d cases out of %d (%.2f%%)\n',...
            tittot,tit1,tit2,nbett,n,nbett/n*100);

end