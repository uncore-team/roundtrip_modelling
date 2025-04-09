% This script checks the behaviour of GOF functions with synthetic data.

clear all
close all

logpath = './tests';

dists = { 'EXP2', 'LN3', 'LL3' };
sizes = [20, 50, 150, 200, 250, 10000];

nsizes = length(sizes);
ndists = length(dists);

ntests = 100; % number of tests with different samples

f = fopen(sprintf('%s/data.csv', logpath), 'w');
print_header(f);

tic
status = 0;
for s = 1:nsizes
    sz = sizes(s);

    for prevmodel = 0:1 % 0: unknown parameters / 1: known parameters
    for d = 1:ndists
        dist = dists{d};

        cont = 0;
        while cont < ntests
            cont = cont + 1;

            if strcmp(dist, 'EXP2')
                params1 = generate_params([1e-3 76e3], [1e-5 10]);

                cont2 = 0;
                while cont2 < ntests
                    cont2 = cont2 + 1;

                    params2 = zeros(size(params1));
                    samples = ExponentialRnd(params1(1), params1(2), 1, sz);
    
                    if prevmodel
                        [reject, stat, thresh] = ExponentialGof(samples, params1(1), params1(2), 1);
                    else
                        [params2(1), params2(2)] = ExponentialFit(samples);
                        % if ~ok
                        %     print_data(f, -1, dist, prevmodel, params1, params2, -1, -1, -1);
                        %     continue;
                        % end
                        [reject, stat, thresh] = ExponentialGof(samples, params2(1), params2(2), 0);
                    end

                    print_data(f, cont, cont2, dist, prevmodel, params1, params2, reject, stat, thresh);
                end

            elseif strcmp(dist, 'LN3')
                params1 = generate_params([1e-4 76e3], [-10 10], [1e-3 8]);

                cont2 = 0;
                while cont2 < ntests
                    cont2 = cont2 + 1;

                    params2 = zeros(size(params1));
                    samples = LognormalRnd(params1(1), params1(2), params1(3), 1, sz);
    
                    if prevmodel
                        [reject, stat, thresh] = LognormalGof(samples, params1(1), params1(2), params1(3), 1);
                    else
                        [ok, params2(1), params2(2), params2(3)] = LognormalFit(samples);
                        if ~ok
                            print_data(f, -1, -1, dist, prevmodel, params1, params2, -1, -1, -1);
                            continue;
                        end
                        [reject, stat, thresh] = LognormalGof(samples, params2(1), params2(2), params2(3), 0);
                    end

                    print_data(f, cont, cont2, dist, prevmodel, params1, params2, reject, stat, thresh);
                end
    
            elseif strcmp(dist, 'LL3')
                params1 = generate_params([1e-4 76e3], [1e-4 32e3], [0.05 0.4]);

                cont2 = 0;
                while cont2 < ntests
                    cont2 = cont2 + 1;

                    params2 = zeros(size(params1));
                    samples = LoglogisticRnd(params1(1), params1(2), params1(3), 1, sz);
    
                    if prevmodel
                        [reject, stat, thresh] = LoglogisticGoF(samples, params1(1), params1(2), params1(3), 1);
                    else
                        [params2(1), params2(2), params2(3), ko] = LoglogisticFit(samples);
                        if ko
                            print_data(f, -1, -1, dist, prevmodel, params1, params2, -1, -1, -1);
                            continue;
                        end
                        [reject, stat, thresh] = LoglogisticGoF(samples, params2(1), params2(2), params2(3), 0);
                    end

                    print_data(f, cont, cont2, dist, prevmodel, params1, params2, reject, stat, thresh);
                end
            else
                error('wrong distribution!.');
            end
        end

        status = status + 100/nsizes/ndists/2;
        fprintf('[%4s][sz=%5d]> %3.1f\n', dist, sz, status);
    end
    end
end
toc
fclose(f);
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
    fprintf(f, 'iter1;iter2;dist;pmodel;p1(1);p1(2);p1(3);p2(1);p2(2);p2(3);reject;stat;thresh\n');
end

function print_data(f, iter1, iter2, dist, pmodel, params1, params2, reject, stat, thresh)

    if numel(params1) == 2, params1(3) = 0; end
    if numel(params2) == 2, params2(3) = 0; end

    fprintf(f, '%d;%d;%s;%d;%f;%f;%f;%f;%f;%f;%f;%f;%f\n', ...
        iter1, iter2, dist, pmodel, ...
        params1(1), params1(2), params1(3), ...
        params2(1), params2(2), params2(3), ...
        reject, stat, thresh);
end
