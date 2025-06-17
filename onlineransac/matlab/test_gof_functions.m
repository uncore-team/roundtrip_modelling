% This script checks the behaviour of GOF functions with synthetic data.

clear all
close all

logpath = './tests';

dists = { 'EXP2', 'LN3', 'LL3' };
%dists = { 'LL3' };
sizes = [20, 50, 150, 200, 250, 10000];

nsizes = length(sizes);
ndists = length(dists);

ntests = 100; % number of tests with different samples

ncombs = 2*nsizes*ndists;
vsize = zeros(2*nsizes*ndists,1);
vdist = cell(2*nsizes*ndists,1);
vkmod = zeros(2*nsizes*ndists,1);

cont = 0;
for s = 1:nsizes
    for d = 1:ndists
        for kmod = 0:1
            cont = cont + 1;
            vsize(cont) = sizes(s);
            vdist(cont) = dists{d};
            vkmod(cont) = kmod;
        end
    end
end

flog = create_logfile(sprintf('%s/logfile.csv', logpath));
print_header(flog, 'size;dist;kmodel;time');

t1 = tic;
% status = 0;
parfor c = 1:ncombs
    sz = vsize(c);
    dist = vdist{c};
    kmod = vkmod(c);

    f = create_logfile(sprintf('%s/data_size=%d_dist=%s_kmodel=%d.csv', logpath, sz, dist, kmod));
%     print_header(f, 'iter1;iter2;dist;pmodel;p1(1);p1(2);p1(3);p2(1);p2(2);p2(3);reject;stat;thresh;time');
    print_header(f, 'iter1;iter2;p1(1);p1(2);p1(3);p2(1);p2(2);p2(3);reject;stat;thresh;time');

    t2 = tic;
    cont = 0;
    while cont < ntests
        cont = cont + 1;

        if strcmp(dist, 'EXP2')
            params1 = generate_params([1e-3 76e3], [1e-5 10]);

            cont2 = 0;
            cont3 = 0;
            while (cont2 < ntests) && (cont3 < 2*ntests)

                params2 = zeros(size(params1));
                samples = ExponentialRnd(params1(1), params1(2), 1, sz);                

                t3 = tic;
                cont3 = cont3 + 1;
                if kmod
                	% CARE! Here the sample can be distorted due to numerical imprecisions
                    [reject, stat, thresh] = ExponentialGof(samples, params1(1), params1(2), kmod);
                else
                    [params2(1), params2(2)] = ExponentialFit(samples);
                    % if ~ok
                    %     print_data(f, cont, -1, dist, knownmodel, params1, params2, -1, -1, -1);
                    %     continue;
                    % end
                    [reject, stat, thresh] = ExponentialGof(samples, params2(1), params2(2), kmod);
                end
                t = toc(t3);
                cont2 = cont2 + 1;
                print_data(f, cont, cont2, dist, kmod, params1, params2, reject, stat, thresh, t);
            end

        elseif strcmp(dist, 'LN3')
            params1 = generate_params([1e-4 76e3], [-10 10], [1e-3 8]);

            cont2 = 0;
            cont3 = 0;
            while (cont2 < ntests) && (cont3 < 2*ntests)

                params2 = zeros(size(params1));
                samples = LognormalRnd(params1(1), params1(2), params1(3), 1, sz);

                t3 = tic;
                cont3 = cont3 + 1;
                if kmod
                	% CARE! Here the sample can be distorted due to numerical imprecisions up to the point
                	% of having values == offset, which is not allowed theoretically (see LognormalRnd)
                    [reject, stat, thresh] = LognormalGof(samples, params1(1), params1(2), params1(3), kmod);
                else
                    [ok, params2(1), params2(2), params2(3)] = LognormalFit(samples);
                    if ~ok
                        t = toc(t3);
                        print_data(f, cont, cont2, dist, kmod, params1, params2, -1, -1, -1, t);
                        continue;
                    end
                    [reject, stat, thresh] = LognormalGof(samples, params2(1), params2(2), params2(3), kmod);
                end
                t = toc(t3);
                cont2 = cont2 + 1;
                print_data(f, cont, cont2, dist, kmod, params1, params2, reject, stat, thresh, t);

            end

        elseif strcmp(dist, 'LL3')
            params1 = generate_params([1e-4 76e3], [1e-4 32e3], [0.05 0.4]);

            cont2 = 0;
            cont3 = 0;
            while (cont2 < ntests) && (cont3 < 2*ntests)

                params2 = zeros(size(params1));
                samples = LoglogisticRnd(params1(1), params1(2), params1(3), 1, sz);

                t3 = tic;
                cont3 = cont3 + 1;
                if kmod
                	% CARE! Here the sample can be distorted due to numerical imprecisions up to the point
                	% of having values == offset, which is not allowed theoretically (see LoglogisticRnd)
                    [reject, stat, thresh] = LoglogisticGoF(samples, params1(1), params1(2), params1(3), kmod);
                else
                    [params2(1), params2(2), params2(3), exitflag] = LoglogisticFit(samples);
                    if exitflag < 0
                        t = toc(t3);
                        print_data(f, cont, cont2, dist, kmod, params1, params2, -1, -1, -1, t);
                        continue;
                    end
                    [reject, stat, thresh] = LoglogisticGoF(samples, params2(1), params2(2), params2(3), kmod);
                end
                t = toc(t3);
                cont2 = cont2 + 1;
                print_data(f, cont, cont2, dist, kmod, params1, params2, reject, stat, thresh, t);

            end
        else
            error('wrong distribution!.');
        end
    end
    t=toc(t2);
    fclose(f);

    % status = status + 100/nsizes/ndists/2;
    fprintf('[sz=%5d][dist=%4s][km=%d][time=%f]> %d/%d\n', sz, dist, kmod, t, c, ncombs);
    print_log(flog, sz, dist, kmod, t);

end
toc(t1)
fclose(flog);

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

function file = create_logfile(filename)
    file = fopen(filename, 'w');
end

function print_header(file, string)
    fprintf(file, '%s\n',string);
end

function print_log(file, size, dist, knownmodel, t)
    fprintf(file, '%d;%s;%d;%f\n', size, dist, knownmodel,t);
end

function print_data(file, iter1, iter2, dist, pmodel, params1, params2, reject, stat, thresh, t)
    if numel(params1) == 2, params1(3) = 0; end
    if numel(params2) == 2, params2(3) = 0; end
    % fprintf(file, '%d;%d;%s;%d;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f\n', ...
    %     iter1, iter2, dist, pmodel, ...
    %     params1(1), params1(2), params1(3), ...
    %     params2(1), params2(2), params2(3), ...
    %     reject, stat, thresh, t);
    fprintf(file, '%d;%d;%f;%f;%f;%f;%f;%f;%f;%f;%f;%f\n', ...
        iter1, iter2, ...
        params1(1), params1(2), params1(3), ...
        params2(1), params2(2), params2(3), ...
        reject, stat, thresh, t);
end


