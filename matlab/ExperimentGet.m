function [fdir,fexp,data] = ExperimentGet(c,n,minind,maxind,transfintospeed,robotspeed,trace)
% Given an experiment identification, return the path and filename for it.
%
% C -> class; a string: 'realoct2023','realpapersensors','sim'
% N -> index, from 1, in the class.
% MININD, MAXIND -> range of data to extract from the scenario; all the
%                   scenario if MININD <= 1, MAXIND >= length(scenario).
% TRANSFINTOSPEED -> 1 to transform the data from time (millisecs) into
%                    meters travelled by a hypothetical robot; 0 to not
%                    transform.
% ROBOTSPEED -> hypothetical robot speed in m/s if TRANSFINTOSPEED == 1;
%               ignored otherwise.
% TRACE -> 1 to see some trace.
%
% FDIR <- folder path, not ended in '/'
% FEXP <- filename with extension
% DATA <- rtts read from the experiment, as a column vector

    if n <= 0
        error('Invalid N');
    end
    if maxind < minind
        error('Invalid range');
    end

    if strcmp(c,'realoct2023')

        fdir = '/DATOS_RSYNC/Investigacion/Papers/web/19 estudio de metodos de modelado de delays para revista con vicente y ana desde proyecto garthim 2023 (pendiente)/paper de modeladores solo/sw/datos experimentos/tomados del roundtripserver';
        fexps = { 'dens1024_local_pcdespacho_02y03deoct2023_1715pm_ff_031023.m' , ...
                  'dens1024_pinos_cylonwifiyfibra_05y06deoct2023_1223am_ch_061023.m', ...
                  'dens1024_smartphone_4Gyfibra_07y08oct2023_1601pm_ff_071023.m', ...
                  'dens65536_local_pcdespacho_09y10deoct2023_1658pm_ff_101023.m', ...
                  'dens65536_pinos_cylonwifiyfibra_20y21deoct2023_1815pm_ff_211023.m',...       
                  'dens65536_smartphone_4Gyfibra_11y12oct2023_1806pm_ch_121023.m', ...
                  'dens786432_local_pcdespacho_16y17deoct2023_1425pm_ff_171023.m', ...
                  'dens786432_pinos_pcdespacho_21y22deoct2023_1827pm_ff_221023.m', ...
                  'dens786432_smartphone_4Gyfibra_19y20oct2023_1946pm_ch_201023.m', ...
                  'dens65536_pinos_laptopwifrepfib_29y30deoct2023_1204am_ff_301023.m', ...
                  'den1024_pinos_laptopwifirepfib_04y05denov2023_847am_ff_051123.m',...
                  'dens786432_pinos_laptopwifirepfib_05y06denov2023_1042_ff_061123.m'};

        if n > length(fexps)
            error('N out of range');
        end
        fexp = fexps{n};
        fn = sprintf('%s/%s',fdir,fexp);
        if trace
            fprintf('Loading experiment [%s]...\r\n',fexp);
        end
        run(fn);
        data = historial(:,2); % only the roundtrip times
    
    elseif strcmp(c,'realpapersensors')

        fdir = '/DATOS_RSYNC/Investigacion/Papers/web/19 estudio de metodos de modelado de delays para revista con vicente y ana desde proyecto garthim 2023 (pendiente)/paper de modeladores solo/sw/datos experimentos/paper sensors open ana gago';
        if n > 18
            error('N out of range');
        end
        fexp = sprintf('sc%d.txt',n);
        fn = sprintf('%s/%s',fdir,fexp);
        if trace
            fprintf('Loading experiment [%s]...\r\n',fexp);
        end
        data = load(fn);

    elseif strcmp(c,'sim')

        fdir = '<sim>';
        fexp = sprintf('%d',n);
        if trace
            fprintf('Generating experiment [%s]...\r\n',fexp);
        end
        rng(54);
        switch n
            case 1
                data = LoglogisticRnd(1000,2,0.3,2000,1);                
            case 2
                data = [ LoglogisticRnd(1000,2,0.25,2000,1) ; ...
                         LoglogisticRnd(1010,10,0.1,500,1) ; ...
                         LoglogisticRnd(950,20,0.4,500,1) ];
            case 3
                data = LognormalRnd(1000,5,0.1,1,2000);
            case 4
                data = ExponentialRnd(1000,5,1,2000);
            otherwise
                error('N out of range');
        end
        
    else
        error('Unknown experiment class');
    end

    data = data(:);
    l = length(data);
    if (minind > 1) || (maxind < l)
        data = data(max(1,minind):min(l,maxind));
    end
    if transfintospeed
        data = data / 1000 * robotspeed; % from milliseconds to meters travelled in each step
    end

end
