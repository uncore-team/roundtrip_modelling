% This script checks the behaviour of GOF functions with synthetic data.

clear all
close all

logpath = './tests';

dists = { 'EXP2', 'LN3', 'LL3' };
ndists = length(dists);

% sizes = [20, 50, 150, 200, 250, 500, 10000];
% sizes = 20:50:1030;
sizes = 1000:500:10000;
nsizes = length(sizes);

ntests = 1000; % number of tests with different samples
fixedtrue = 0; % 0- use random models each time; 1- use always the true model

flog = create_logfile(sprintf('%s/logfile.csv', logpath));
print_data(flog, 'dist;size;kmalpha;kmbeta;umalpha;umbeta;unfit;time');

v_km_alphas = zeros(ndists,nsizes);
v_um_alphas = zeros(ndists,nsizes);
v_unfitrates = zeros(ndists,nsizes);

t1 = tic;
for d = 1:ndists
    dist = dists{d};
    [numparms, namparms] = ModelParmsDef(dist);

    parfor s = 1:nsizes
        sz = sizes(s);
        tmodel = NaN;

        f = create_logfile(sprintf('%s/data_dist=%s_size=%d.csv', logpath, dist, sz));
        print_data(f,...
            'iter;tc(1);tc(2);tc(3);treject;tstat;fc(1);fc(2);fc(3);freject;fstat;time');

        numunfit = 0;
        km_nrejects = 0;
        um_nrejects = 0;

        if fixedtrue
            tmodel = ModelCreateRnd(dist, 'typrnd'); % generate model
        end

        t2 = tic;
        for n = 1:ntests
            t3 = tic;

            if ~fixedtrue
                tmodel = ModelCreateRnd(dist, 'typrnd'); % generate model
            end

            sample = ModelRnd(tmodel, 1, sz); % generate sample
    
            [treject, tstat, ~] = ModelGof(tmodel, sample, 1);
            km_nrejects = km_nrejects + treject;
            tcoeffs = ModelToCoeffs(tmodel);
    
            fmodel = ModelFit(sample, 1, length(sample), dist);
            if fmodel.defined
                [freject, fstat, ~] = ModelGof(fmodel, sample, 0); % gof with parms coming from the sample
                um_nrejects = um_nrejects + freject;
                fcoeffs = ModelToCoeffs(fmodel);
            else
                freject = 2;
                fstat = NaN;
                fcoeffs = zeros(size(tcoeffs));
                numunfit = numunfit + 1;
            end
            t = toc(t3);
    
            buffer = sprintf('%d;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f\n', ...
                n, ...
                tcoeffs(2), tcoeffs(3), tcoeffs(4), treject, tstat, ...
                fcoeffs(2), fcoeffs(3), fcoeffs(4), freject, fstat, ...
                t);
            print_data(f, buffer);
        end
        t = toc(t2);
        fclose(f);
    
%         fprintf('\n\nALPHA (TYPE I ERROR) ESTIMATES:\n');
        fprintf('[dist=%4s][sz=%5d][time=%f]\n', dist, sz, t);
    
        unfitrate = numunfit/ntests;
        km_alpha = km_nrejects/ntests;
        um_alpha = um_nrejects/(ntests - numunfit);
    
%         fprintf('Assuming parameters known:\n');
%         fprintf('\tEst.alpha (Type I error): %f\n', km_alpha);
%         fprintf('\tCorrect detection: %f\n', km_beta);
%         fprintf('\n');
%         fprintf('Assuming parameters unknown (undefined: %d; %.2f%%):\n', numunfit, unfitrate);
%         fprintf('\tEst.alpha (Type I error): %f\n',um_alpha);
%         fprintf('\tCorrect detection: %f\n',um_beta);
    
        buffer = sprintf('%s;%d;%f;%f;%f;%f;%f;%f\n',...
            dist, sz, km_alpha, 1-km_alpha, um_alpha, 1-um_alpha, unfitrate, t);
        print_data(flog, buffer);
    
        v_km_alphas(d, s) = km_alpha;
        v_um_alphas(d, s) = um_alpha;
        v_unfitrates(d, s) = unfitrate;

    end

end
toc(t1)
fclose(flog);

figure;
for d=1:ndists

    subplot(3,1,d)
    plot(sizes, v_km_alphas(d,:), 'b*-',...
         sizes, 1-v_km_alphas(d,:), 'bx-',...
         sizes, v_um_alphas(d,:), 'r*-',...
         sizes, 1-v_um_alphas(d,:), 'rx-',...
         sizes, v_unfitrates(d,:), 'ko-',...
         [sizes(0) sizes(end)], [.05 .05], 'k-.',...
         [sizes(0) sizes(end)], [.95 .95], 'k-.');
    subtitle(['dist=' dists{d}]);
    xlabel('sizes'), ylabel('alpha / beta');
    legend('alpha (known model)', ...
           'beta (known model)', ...
           'alpha (unknown model)', ...
           'beta (unknown model)', ...
           'rate of unfits', ...
           '.05', '.95');
    drawnow;

end

return;

function file = create_logfile(filename)
    file = fopen(filename, 'w');
end

function print_data(file, string)
    fprintf(file, '%s\n',string);
end
