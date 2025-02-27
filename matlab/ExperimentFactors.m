function [levsstr,levsnum] = ExperimentFactors(expclass,expindex,...
                                               algorithm,algparms)
% Given an experiment already launched where some performance measures have
% been collected, return the factor levels that define that experiment.
%
% EXPCLASS -> class of the experiment (see ExperimentCatalog).
% EXPINDEX -> index within the class (see ExperimentCatalog).
% ALGORITHM -> modelling algorithm used. One name in the catalog of
%              algorithms (see AlgCatalog).
% ALGPARMS -> a struct with the parameters used for that algorithm.
%
% LEVSSTR <- a cell with as many elements as factors; each one will be the 
%            level of that factor in this experiment.
% LEVSNUM <- a row vector with as many elements as factors; each one will
%            be the level index within the levels definition of that factor
%            for this experiment.

    ca = ExperimentCatalog(0);
    [fds,~,~,~] = ExperimentFactorDefs(0);
    
    % NOTE: ExperimentFactorDefs must change if these indexes in FDS
    % change.

    indcat = ExperimentCatalogFind(ca,expclass,expindex);

    % factors that depend on the scenario

    indfactordist = find_string_in_cell(fds{1}.levels, ca{indcat}.dist);
    indfactorserver = find_string_in_cell(fds{2}.levels, ca{indcat}.serversw);
    indfactorclient = find_string_in_cell(fds{3}.levels, ca{indcat}.clientsw);
    indfactornet = find_string_in_cell(fds{4}.levels, ca{indcat}.netw);
    indfactordens = find_string_in_cell(fds{5}.levels, ca{indcat}.dens);

    % factors that depend on the algorithm

    indfactoralg = find_string_in_cell(fds{6}.levels, algorithm);

    % factors that depend on the algorithm parameters

    if strcmp(algorithm,'onlineransac')
        switch algparms.modelpreserving
            case 0
                indfactorransmodelpres = find_string_in_cell(fds{7}.levels, 'no-preserv');
            case 1
                indfactorransmodelpres = find_string_in_cell(fds{7}.levels, 'preserv-rebuild');
            case 2
                indfactorransmodelpres = find_string_in_cell(fds{7}.levels, 'rebuild-preserv');
            otherwise
                error('Invalid model preserving');
        end
        switch algparms.samplesliding
            case 0
                indfactorranssamplesld = find_string_in_cell(fds{8}.levels, 'no-slide');
            case 1
                indfactorranssamplesld = find_string_in_cell(fds{8}.levels, 'slide');
            otherwise
                error('Invalid sample sliding');
        end
        switch algparms.datapreserving
            case 0
                indfactorransdatapres = find_string_in_cell(fds{9}.levels, 'no-preserv');
            case 1
                indfactorransdatapres = find_string_in_cell(fds{9}.levels, 'preserv');
            otherwise
                error('Invalid data preserving');
        end
        tystxt = '';
        for f = 1:length(algparms.mtypes)
            if f == 1
                tystxt = algparms.mtypes{f};
            else
                tystxt = sprintf('%s-%s',tystxt,algparms.mtypes{f});
            end
        end
        indfactorransmodelform = find_string_in_cell(fds{10}.levels, tystxt);
    elseif strcmp(algorithm,'bernoulli')
        indfactorransmodelpres = find_string_in_cell(fds{7}.levels, 'noransac');
        indfactorranssamplesld = find_string_in_cell(fds{8}.levels, 'noransac');
        indfactorransdatapres = find_string_in_cell(fds{9}.levels, 'noransac');
        indfactorransmodelform = find_string_in_cell(fds{10}.levels, 'noransac');
    else
        error('Unknown algorithm');
    end

    % collecting everything

    levsstr = {fds{1}.levels{indfactordist},...
               fds{2}.levels{indfactorserver},...
               fds{3}.levels{indfactorclient},...
               fds{4}.levels{indfactornet},...
               fds{5}.levels{indfactordens},...
               fds{6}.levels{indfactoralg},...
               fds{7}.levels{indfactorransmodelpres},...
               fds{8}.levels{indfactorranssamplesld},...
               fds{9}.levels{indfactorransdatapres},...
               fds{10}.levels{indfactorransmodelform}};
    levsnum = [indfactordist,...
               indfactorserver,...
               indfactorclient,...
               indfactornet,...
               indfactordens,...
               indfactoralg,...
               indfactorransmodelpres,...
               indfactorranssamplesld,...
               indfactorransdatapres,...
               indfactorransmodelform];

end