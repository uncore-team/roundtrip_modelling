function [fs,factornames,factorlevels,totlevs] = ExperimentFactorDefs(trace)
% Return all the factors used for experiments, and the possible levels 
% of each one.
%
% TRACE -> 1 to show in console factors and levels.
%
% FS <- a cell with as many elements as factors. Each one a struct:
%           .name <- factor name
%           .levels <- a cell with the names of the levels of the factors,
%                      from 1 to the last one.
% FACTORNAMES <- a cell with the names of all factors.
% FACTORLEVELS <- a row vector with the number of levels of each factor.
% TOTLEVS <- total number of levels when considering all factors.
    
    ca = ExperimentCatalog(0);

    fs = {};   

    % NOTE: the numbering within FS of the factors must change if
    % ExperimentFactors changes.

    % Factors related to the scenario (those that have some monotonic
    % meaning are ordered in increasing order of their names)

    fs{1} = struct('name','distance',...
                   'levels',{get_uniques_in_cell(get_structfields_as_cell(ca,{'dist'},0))});
    fs{1}.levels = sort(fs{1}.levels);
    fs{2} = struct('name','server-sw',...
                   'levels',{get_uniques_in_cell(get_structfields_as_cell(ca,{'serversw'},0))});
    fs{3} = struct('name','client-sw',...
                   'levels',{get_uniques_in_cell(get_structfields_as_cell(ca,{'clientsw'},0))});
    fs{4} = struct('name','network',...
                   'levels',{get_uniques_in_cell(get_structfields_as_cell(ca,{'netw'},0))});
    fs{4}.levels = sort(fs{4}.levels);
    fs{5} = struct('name','density',...
                   'levels',{get_uniques_in_cell(get_structfields_as_cell(ca,{'dens'},1))});
    denslevs = str2double(fs{5}.levels);
    [~, sorted_denslevs] = sort(denslevs);
    fs{5}.levels = fs{5}.levels(sorted_denslevs);

    % Factors related to the algorithm

    fs{6} = struct('name','algorithm',...
                   'levels',{{'onlineransac','bernoulli'}});

    % Factors related to the onlineransac algorithm

    % the following names must follow the corresponding numeric order of parms in AlgOnlineransac
    % 'noransac' must be as it is and always be the first
    fs{7} = struct('name','onlineransac-modelpreserv',...
                   'levels',{{'noransac','no-preserv','preserv-rebuild','rebuild-preserv'}});
    fs{8} = struct('name','onlineransac-sampleslide',...
                   'levels',{{'noransac','no-slide','slide'}});
    fs{9} = struct('name','onlineransac-datapreserv',...
                    'levels',{{'noransac','no-preserv','preserv'}});
    fs{10}= struct('name','onlineransac-modelforms',...
                   'levels',{{'noransac',...
                              'LL3','EXP2','LN3',...
                              'LL3-EXP2','LL3-LN3','EXP2-LL3','EXP2-LN3','LN3-LL3','LN3-EXP2',...
                              'LL3-EXP2-LN3','LL3-LN3-EXP2','EXP2-LL3-LN3','EXP2-LN3-LL3','LN3-LL3-EXP2','LN3-EXP2-LL3'}});

    % factor names and number of levels

    numfactors = length(fs);
    factornames = cell(1,numfactors);
    factorlevels = zeros(1,numfactors);
    totlevs = 1;
    if trace
        fprintf('FACTORS (%d) AND LEVELS:\n',numfactors);
    end
    for f = 1:numfactors
        factornames{f} = fs{f}.name;
        factorlevels(f) = length(fs{f}.levels);
        totlevs = totlevs * length(fs{f}.levels);
        if trace
            fprintf('\t%s (%d levels)\n',fs{f}.name,length(fs{f}.levels));
        end
    end
    if trace
        fprintf('TOTAL No. OF LEVELS: %d\n',totlevs);
    end

end