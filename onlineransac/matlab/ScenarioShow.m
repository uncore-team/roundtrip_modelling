function hfigesc = ScenarioShow(S,tit,regs,models,units)
% Show common figures about the scenario S.
%
% S -> entire scenario
% TIT -> title
% REGS -> regimes detected (see regimeDetection); empty if none.
% MODELS -> a cell with as many entries as models detected, each entry
%           being a struct with:
%               .index -> index within S when the model was assessed.
%               .model -> as in ModelCreate.m.
%           Empty if no model detected.
% UNITS -> text with the units of the values in S.
%
% HFIGESC <- figure handler for the one that shows the scenario.

    ts = cumsum(S);

    hfigesc = figure;
    subplot(2,1,2);
    plot(S,'.-');
    grid;
    xlabel('index');
    subplot(2,1,1);
    hold on;
    grid;
    plotregs(ts,S,regs);
    plot(ts,S,'.-');
    xlabel(units);
    ylabel('round-trip time (ms)');
    title(sprintf('Scenario [%s] - sequence of %d roundtrip times',tit,length(S)));

    figure;
    drawHisto(S,'Moments',units);

    figure;
    autocorr(S);

    if ~isempty(regs)
        figure;
        ss = regs(:,2)-regs(:,1) + 1;
        histogram(ss,size(regs,1));
        xlabel('Length of regime');
        ylabel('freq');
        title(sprintf('Length of regimes (%.2f%% > 20)',length(find(ss>20)) / size(regs,1) * 100));
    end

    if ~isempty(models)
        mnames = ModelTypes(1);
        nmodels = length(mnames);
        trajs = cell(1,nmodels); % each one with the trajectory of each model
        for f = 1:length(models)
            if ~models{f}.model.defined
                error('Model %d not defined',f);
            end
            indtraj = NaN;
            for g = 1:nmodels
                if strcmp(mnames{g},models{f}.model.type)
                    indtraj = g;
                    break;
                end
            end
            if isnan(indtraj)
                error('Model %s not known',models{f}.model.type);
            end

            mps = ModelToCoeffs(models{f}.model);
            trajs{indtraj} = [trajs{indtraj} ; mps(2:end)];
        end
        for f = 1:nmodels
            if ~isempty(trajs{f})
                figure;
                plot3(trajs{f}(:,1),trajs{f}(:,2),trajs{f}(:,3),'.-');
                hold on;
                plot3(trajs{f}(1,1),trajs{f}(1,2),trajs{f}(1,3),'go');
                plot3(trajs{f}(end,1),trajs{f}(end,2),trajs{f}(end,3),'r*');
                xlabel('parm-1');
                ylabel('parm-2');
                zlabel('parm-3');
                grid;
                title(sprintf('Trajectory of models %s',mnames{f}));
            end
        end
    end

end

function plotregs(ts,S,regs)

    if isempty(S) || isempty(ts)
        return;
    end

    N = length(S);
    [nr,~] = size(regs);

    if nr == 0 % the entire data sequence is a non-regime
        drawOneRegime(0,ts,S,1,N,[]);
    else
        if regs(1,1) > 1 % heading non-regime
            drawOneRegime(0,ts,S,1,regs(1,1)-1,[]);
        end
        for f = 1:nr
            indini = regs(f,1);
            indfini = regs(f,2);
            if (f > 1) && (indini > regs(f-1,2)+1) % gap before
                drawOneRegime(0,ts,S,regs(f-1,2)+1,indini-1,[]);
            end
            model = ModelFromCoeffs(regs(f,3:end));
            drawOneRegime(f,ts,S,indini,indfini,model);
        end
        if regs(end,2) < N % trailing non-regime
            drawOneRegime(0,ts,S,regs(end,2)+1,N,[]);
        end
    end

end

function drawOneRegime(indreg,ts,S,indini,indfini,model)

    medy = mean(S);
    stdy = std(S);
    h = 12 * stdy;
    minry = medy - h/2;
    heightr = h;
    resolpdf = 100;

    datareg = S(indini:indfini);

    tini = ts(indini);
    tfini = ts(indfini);
    if ~isempty(model)
        co = [77,100,71]/100;
    else
        co = [100,77,71]/100;
    end
    rectangle('Position',[tini,minry,tfini-tini,heightr],...
              'FaceColor',co,...
              'EdgeColor', 'none');
    if indreg > 0
        txt = sprintf('#%d',indreg);
    else
        txt = 'noreg';
    end
    text(tini,minry,sprintf('%s (%d-%d: %d rtts; %f s))',txt,indini,indfini,indfini-indini+1,(tfini-tini)/1000));

    minds = min(datareg);
    maxds = max(datareg);
    [nhist,edgeshist] = histcounts(datareg);

    xdatareg = linspace(minds,maxds,resolpdf);
    xdrawings = linspace(tini,tfini,resolpdf);

    midxhist = (edgeshist(2:end) - edgeshist(1:end-1)) / 2 + ...
               edgeshist(1:end-1);
    xplothist = (midxhist - minds) / (maxds - minds) * ...
                    (xdrawings(end) - xdrawings(1)) + ...
                xdrawings(1);
    yplothist = minry + heightr + ...
                nhist / max(nhist) * heightr / 2;
    
    modelbuilthere = 0;
    if ~isempty(model)
        ypdf = ModelPdf(model,S,xdatareg);
        yplotpdf = minry + heightr + ypdf / max(ypdf) * heightr / 2;
        xmode = ModelToMode(model,S);
    elseif length(datareg) >= 20
        % provide a LL3 calculated here just as reference
        [a, b, c, exitflag] = LoglogisticFit(datareg);
        if exitflag > 0
            modelbuilthere = 1;
            ypdf = LoglogisticPdf(xdatareg,a,b,c);
            yplotpdf = minry + heightr + ypdf / max(ypdf) * heightr / 2;
            xmode = LoglogisticToMode(a,b,c);            
            [reject,stat,~] = LoglogisticGoF(datareg,a,b,c,0);
        else
            xmode = mode(datareg);
        end
    else
        xmode = mode(datareg);
    end
    xmode = (xmode - xdatareg(1)) / (xdatareg(end) - xdatareg(1)) * (xdrawings(end) - xdrawings(1)) + xdrawings(1);

    bar(xplothist,yplothist,'BaseValue',minry + heightr); %,'FaceAlpha',0.5);
    if ~isempty(model)
        txt = sprintf('#%d [%s]',indreg,model.type);
    else
        if modelbuilthere
            if ~reject
                txt = 'yes';
            else
                txt = 'no';
            end
        else
            txt = '?';
        end
    end
    text(xmode,max(yplothist),txt);
    if ~isempty(model) || modelbuilthere
        if modelbuilthere
            co = 'r-';
            if ~reject
                txtparms = sprintf('%s: %.2f, %.2f, %.2f (accept)',txt,a,b,c);
            else
                txtparms = sprintf('%s: %.2f, %.2f, %.2f (reject; stat %f)',txt,a,b,c,stat);
            end
        else
            co = sprintf('%c-',ModelColor(model));
            mps = ModelToCoeffs(model);
            t = '';
            for f = 1:length(mps)
                if ~isnan(mps(f))
                    if f > 1
                        t = sprintf('%s,',t);
                    end
                    t = sprintf('%s%.2f',t,mps(f));
                end
            end
            txtparms = sprintf('%s: %s',txt,t);
        end
        plot(xdrawings,yplotpdf,co,'LineWidth',2);
        text(xmode,(max(yplotpdf)-min(yplotpdf))/2+min(yplotpdf),txtparms);
    end

end