function [ok,pm] = ReqCheckInModel(r,m,scenario,withfigs)
% Check whether the requirement R is satisfied by model M.
%
% R -> Requirement (see ReqCreate).
% M -> Model (see ModelCreate).
% SCENARIO -> complete scenario (required for some models).
% WITHFIGS -> 1 to illustrate the process.
%
% OK <- 1 if satisfied, 0 if not.
% PM <- prob of satisfying the req in the model.

    if ~m.defined
        error('Undefined model');
    end

    vs = ModelCdf(m,scenario,[r.t0,r.t1]);
    cdf0 = vs(1);
    cdf1 = vs(2);
    if withfigs
        figure;
        xs = linspace(r.t0 - r.t0/8,r.t1 + r.t0/8,10000);
        ys = ModelPdf(m,scenario,xs);
        plot(xs,ys,'r-','LineWidth',2);
        grid;
        hold on;
        plot(r.t0 * [1 1],[0,max(ys)],'g--');
        plot(r.t1 * [1 1],[0,max(ys)],'m--');
        text(r.t0,max(ys),sprintf('cdf0 = %f',cdf0));
        text(r.t1,max(ys),sprintf('cdf1 = %f',cdf1));
        legend('pdf','t0','t1');
        xlabel('Rtt');
        ylabel('pdf');
        title(sprintf('Req. satisfaction. Prob = %f',cdf1-cdf0));
    end
    pm = cdf1 - cdf0;
    if (pm >= r.p)
        ok = 1;
    else
        ok = 0;
    end

end