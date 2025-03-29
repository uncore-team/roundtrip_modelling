% This script checks the behaviour of GOF functions with synthetic data.

clear all
close all

dists = { 'EXP2', 'LN3', 'LL3' };

sizes = [20, 50, 150, 200, 250];

nsizes = length(sizes);
ndists = length(dists);
ntests = 100;

fig = figure;

for s = 1:nsizes
    sz = sizes(s);

    % buffers
    data = zeros(ndists,2);
    data2 = zeros(ndists,2);

    for d = 1:ndists
        dist = dists{d};

        % f  = fopen(sprintf('./tests/%s_sz%d_true.txt', dist, sz), 'w');
        % print_header(f);
        % f2 = fopen(sprintf('./tests/%s_sz%d_false.txt', dist, sz), 'w');
        % print_header(f2)

        for n = 1:ntests

            if strcmp(dist, 'EXP2')
                params  = generate_params([1e-4 76e3], [1e-5 10]); % beta=1/mean => mean in [eps, 5 times 2 seconds]
                params2 = generate_params([1e-4 76e3], [1e-5 10]);

                samples = ExponentialRnd(params(1), params(2), 1, sz);
                [reject1, stat1, thresh1] = ExponentialGof(samples, params(1), params(2));
                [reject2, stat2, thresh2] = ExponentialGof(samples, params2(1), params2(2));

            elseif strcmp(dist, 'LN3')
                params  = generate_params([1e-4 76e3], [-10 10], [1e-3 8]);
                params2 = generate_params([1e-4 76e3], [-10 10], [1e-3 8]);

                samples = LognormalRnd(params(1), params(2), params(3), 1, sz);
                [reject1, stat1, thresh1] = LognormalGof(samples, params(1), params(2), params(3));
                [reject2, stat2, thresh2] = LognormalGof(samples, params2(1), params2(2), params2(3));

            elseif strcmp(dist, 'LL3')
                params  = generate_params([1e-4 76e3], [1e-4 32e3], [0.05 0.4]);
                params2 = generate_params([1e-4 76e3], [1e-4 32e3], [0.05 0.4]);

                samples = LoglogisticRnd(params(1), params(2), params(3), 1, sz);
                [reject1, stat1, thresh1] = LoglogisticGoF(samples, params(1), params(2), params(3));
                [reject2, stat2, thresh2] = LoglogisticGoF(samples, params2(1), params2(2), params2(3));

            else
                error('wrong distribution!.');
            end

            % print_data(f, params, params, reject1, stat1, thresh1);
            % print_data(f2, params, params2, reject2, stat2, thresh2);

            data(d,1) = data(d,1) + (1 - reject1)*100/ntests;
            data(d,2) = data(d,2) + 100 - (1 - reject1)*100/ntests;
            data2(d,1) = data2(d,1) + reject2*100/ntests;
            data2(d,2) = data2(d,2) + 100 - reject2*100/ntests;

        end

        % fclose(f);
        % fclose(f2);
    end

    subplot(2,nsizes,s)
    bar(data,'stacked');
    subtitle(['GoF with true params (size=' num2str(sz) ')']);
    xlabel('Distribution'); ylabel('GoF');
    legend('Success=reject0', 'Failure=reject1');
    xticklabels(dists);
    ylim([0,100]);

    subplot(2,nsizes,nsizes+s)
    bar(data2,'stacked');
    subtitle(['GoF with false params (size=' num2str(sz) ')']);
    xlabel('Distribution'); ylabel('GoF');
    legend('Success=reject1', 'Failure=reject0');
    xticklabels(dists);
    ylim([0,100]);

    % data
end

return;

function num = rnd(xmin, xmax)
    num = xmin + (xmax - xmin) * rand();
end 

function params = generate_params(alim, blim, clim)

    if nargin == 0
        error('wrong params number!');
    end
    if nargin >= 1
        params(1) = rnd(alim(1), alim(2));
    end
    if nargin >= 2
        params(2) = rnd(blim(1), blim(2));
    end
    if nargin == 3
        params(3) = rnd(clim(1), clim(2));
    end
end 

function print_header(f)
    fprintf(f, '%% p1(1)\t\tp1(2)\t\tp1(3)\t\tp2(1)\t\tp2(2)\t\tp2(3)\t\treject\t\tstat\t\tthresh\n');
end

function print_data(f, params, params2, reject, stat, thresh)
    if numel(params) == 2, params(3) = 0; end
    if numel(params2) == 2, params2(3) = 0; end
    fprintf(f, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', ...
        params(1),  params(2),  params(3), ...
        params2(1), params2(2), params2(3), ...
        reject, stat, thresh);
end


% dists = { create_model('EXP2', 2, @ExponentialRnd, @ExponentialGof)
%           create_model('LN3', 3, @LognormalRnd, @LognormalGof)
%           create_model('LL3', 3, @LoglogisticRnd, @LoglogisticGoF) };
% 
% sizes = [10, 50, 100, 500, 1000, 5000];
% 
% nsizes = length(sizes);
% ndists = length(dists);
% ntests = 100;
% myeps = 1e-9;
% 
% fig = figure;
% 
% for s = 1:nsizes
%     sz = sizes(s);
% 
%     % buffers
%     data = zeros(ndists,2);
%     data2 = zeros(ndists,2);
%     labels = cell(1,ndists);
% 
%     for d = 1:ndists
%         dist = dists{d};
%         labels{d} = dist.name;
% 
%         f  = fopen(sprintf('./tests/%s_sz%d_true.txt', dist.name, sz), 'w');
%         print_header(f);
%         f2 = fopen(sprintf('./tests/%s_sz%d_false.txt', dist.name, sz), 'w');
%         print_header(f2)
% 
%         for n = 1:ntests
% 
%             % params used to generate & check the samples
%             params = [ rnd(0.01, 2)          % 0.01 <= location <= 2
%                        rnd(myeps, 10-myeps)      % 0 < scale? < 10
%                        rnd(myeps, 0.5-myeps) ];  % 0 < shape? < 0.5
% 
%             % params used to check the samples (only)
%             params2 = [ rnd(0.01, 2)          % 0.01 <= location <= 2
%                         rnd(myeps, 10-myeps)      % 0 < scale? < 10
%                         rnd(myeps, 0.5-myeps) ];  % 0 < shape? < 0.5
% 
%             if dist.nparams == 2
%                 samples = dist.rnd(params(1), params(2), 1, sz);
%                 [reject1, stat1, thresh1] = dist.gof(samples, params(1), params(2));
%                 [reject2, stat2, thresh2] = dist.gof(samples, params2(1), params2(2));
%                 params(3) = 0; params2(3) = 0; % for printing purposes only
% 
%             elseif dist.nparams == 3
%                 samples = dist.rnd(params(1), params(2), params(3), 1, sz);
%                 [reject1, stat1, thresh1] = dist.gof(samples, params(1), params(2), params(3));
%                 [reject2, stat2, thresh2] = dist.gof(samples, params2(1), params2(2), params2(3));
% 
%             else
%                 error('number of params wrong!.');
%             end
% 
%             print_data(f, params, params, reject1, stat1, thresh1);
%             print_data(f2, params, params2, reject2, stat2, thresh2);
% 
%             data(d,1) = data(d,1) + (1 - reject1)*100/ntests;
%             data(d,2) = data(d,2) + 100 - (1 - reject1)*100/ntests;
%             data2(d,1) = data2(d,1) + reject2*100/ntests;
%             data2(d,2) = data2(d,2) + 100 - reject2*100/ntests;
% 
%         end
% 
%         fclose(f);
%         fclose(f2);
%     end
% 
%     subplot(2,nsizes,s)
%     bar(data,'stacked');
%     subtitle(['GoF with true params (size=' num2str(sz) ')']);
%     xlabel('Distribution'); ylabel('GoF');
%     legend('Success=reject0', 'Failure=reject1');
%     xticklabels(labels);
%     ylim([0,100]);
% 
%     subplot(2,nsizes,nsizes+s)
%     bar(data2,'stacked');
%     subtitle(['GoF with false params (size=' num2str(sz) ')']);
%     xlabel('Distribution'); ylabel('GoF');
%     legend('Success=reject1', 'Failure=reject0');
%     xticklabels(labels);
%     ylim([0,100]);
% 
%     % data
% end
% 
% return;
% 
% function num = rnd(xmin, xmax)
%     num = xmin + (xmax - xmin) * rand();
% end 
% 
% function model = create_model(name, nparams, frnd, fgof)
%     model.name = name;
%     model.nparams = nparams;
%     model.rnd = frnd;
%     model.gof = fgof;
% end
% 
% function print_header(f)
%     fprintf(f, '%% p1(1)\t\tp1(2)\t\tp1(3)\t\tp2(1)\t\tp2(2)\t\tp2(3)\t\treject\t\tstat\t\tthresh\n');
% end
% 
% function print_data(f, params, params2, reject, stat, thresh)
%     % if numel(params) == 2, params(3) = 0; end
%     % if numel(params2) == 2, params2(3) = 0; end
%     fprintf(f, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', ...
%         params(1),  params(2),  params(3), ...
%         params2(1), params2(2), params2(3), ...
%         reject, stat, thresh);
% end
