% This script checks the behaviour of GOF functions with synthetic data.

sizes = [10, 50, 100, 500, 1000, 5000];

exp2.name = 'EXP2';
exp2.nparams = 2;
exp2.rnd = @ExponentialRnd;
exp2.gof = @ExponentialGof;

ln3.name = 'LN3';
ln3.nparams = 3;
ln3.rnd = @LognormalRnd;
ln3.gof = @LognormalGof;

ll3.name = 'LL3';
ll3.nparams = 3;
ll3.rnd = @LoglogisticRnd;
ll3.gof = @LoglogisticGoF;

dists = {exp2, ln3, ll3};

ntests = 500;
for d = 1:length(dists)

    dist = dists{d};
    for s = sizes

        filename = sprintf('./tests/%s_sz%d_true.txt', dist.name, s);
        f = fopen(filename, 'w');
        fprintf(f, '%% p1(1)\t\tp1(2)\t\tp1(3)\t\tp2(1)\t\tp2(2)\t\tp2(3)\t\treject\t\tstat\t\tthresh\n');

        filename = sprintf('./tests/%s_sz%d_false.txt', dist.name, s);
        f2 = fopen(filename, 'w');
        fprintf(f2, '%% p1(1)\t\tp1(2)\t\tp1(3)\t\tp2(1)\t\tp2(2)\t\tp2(3)\t\treject\t\tstat\t\tthresh\n');

        if dist.nparams == 2

            for n = 1:ntests
                params = rand(1, dist.nparams)*10; % params used to generate & check the samples
                samples = dist.rnd(params(1), params(2), 1, s);
                [reject, stat, thresh] = dist.gof(samples, params(1), params(2));
                print_data(f, params, params, reject, stat, thresh);

                params2 = rand(1, dist.nparams)*10; % params used to calculate Gof (false)
                [reject, stat, thresh] = dist.gof(samples, params2(1), params2(2));
                print_data(f2, params, params2, reject, stat, thresh);
            end

        elseif dist.nparams == 3

            for n = 1:ntests
                params = rand(1, dist.nparams)*10; % params used to generate & check the samples
                samples = dist.rnd(params(1), params(2), params(3), 1, s);

                [reject, stat, thresh] = dist.gof(samples, params(1), params(2), params(3));
                print_data(f, params, params, reject, stat, thresh);

                params2 = rand(1, dist.nparams)*10; % params used to calculate Gof (false)
                [reject, stat, thresh] = dist.gof(samples, params2(1), params2(2), params2(3));
                print_data(f2, params, params2, reject, stat, thresh);
            end

        else
            error('number of params wrong!.');
        end

        fclose(f);
        fclose(f2);
    end
end

return

function print_data(f, params, params2, reject, stat, thresh)

    if numel(params) == 2, params(3) = 0; end
    if numel(params2) == 2, params2(3) = 0; end

    fprintf(f, '%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n', ...
        params(1),  params(2),  params(3), ...
        params2(1), params2(2), params2(3), ...
        reject, stat, thresh);

end