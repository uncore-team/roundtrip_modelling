% This script checks the behaviour of GOF functions with synthetic data.

clear all
close all

dists = { create_model('EXP2', 2, @ExponentialRnd, @ExponentialGof)
          create_model('LN3', 3, @LognormalRnd, @LognormalGof)
          create_model('LL3', 3, @LoglogisticRnd, @LoglogisticGoF) };

sizes = [10, 50, 100, 500, 1000, 5000];

nsizes = length(sizes);
ndists = length(dists);
ntests = 100;

fig = figure;

for s = 1:nsizes
    sz = sizes(s);

    % buffers
    data = zeros(ndists,2);
    data2 = zeros(ndists,2);
    labels = cell(1,ndists);

    for d = 1:ndists
        dist = dists{d};
        labels{d} = dist.name;

        % filename = sprintf('./tests/%s_sz%d_true.txt', dist.name, sz);
        % f = fopen(filename, 'w');
        % fprintf(f, '%% p1(1)\t\tp1(2)\t\tp1(3)\t\tp2(1)\t\tp2(2)\t\tp2(3)\t\treject\t\tstat\t\tthresh\n');
        % 
        % filename = sprintf('./tests/%s_sz%d_false.txt', dist.name, sz);
        % f2 = fopen(filename, 'w');
        % fprintf(f2, '%% p1(1)\t\tp1(2)\t\tp1(3)\t\tp2(1)\t\tp2(2)\t\tp2(3)\t\treject\t\tstat\t\tthresh\n');

        for n = 1:ntests
            params = rand(1, dist.nparams)*10; % params used to generate & check the samples
            params2 = rand(1, dist.nparams)*10; % params used to calculate Gof (false)

            if dist.nparams == 2
                samples = dist.rnd(params(1), params(2), 1, sz);
                [reject1, stat1, thresh1] = dist.gof(samples, params(1), params(2));
                [reject2, stat2, thresh2] = dist.gof(samples, params2(1), params2(2));

            elseif dist.nparams == 3
                samples = dist.rnd(params(1), params(2), params(3), 1, sz);
                [reject1, stat1, thresh1] = dist.gof(samples, params(1), params(2), params(3));
                [reject2, stat2, thresh2] = dist.gof(samples, params2(1), params2(2), params2(3));

            else
                error('number of params wrong!.');
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
    xticklabels(labels);
    ylim([0,100]);

    subplot(2,nsizes,nsizes+s)
    bar(data,'stacked');
    subtitle(['GoF with false params (size=' num2str(sz) ')']);
    xlabel('Distribution'); ylabel('GoF');
    legend('Success=reject1', 'Failure=reject0');
    xticklabels(labels);
    ylim([0,100]);

    % data
end

return;

% function print_data(f, params, params2, reject, stat, thresh)
% 
%     if numel(params) == 2, params(3) = 0; end
%     if numel(params2) == 2, params2(3) = 0; end
% 
%     fprintf(f, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', ...
%         params(1),  params(2),  params(3), ...
%         params2(1), params2(2), params2(3), ...
%         reject, stat, thresh);
% 
% end

function model = create_model(name, nparams, frnd, fgof)

    model.name = name;
    model.nparams = nparams;
    model.rnd = frnd;
    model.gof = fgof;

end