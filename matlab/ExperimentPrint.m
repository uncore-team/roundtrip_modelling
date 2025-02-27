function ExperimentPrint(ca,ind,withfig)
% Given an experiment, print in console its main content.
%
% CA -> catalog of experiments produced by ExperimentCatalog.
% IND -> index in the catalog of the selected experiment.
% WITHFIG -> 1 to plot the data of the 

    if ~iscell(ca)
        error('CA does not seem to be an experiment catalog');
    end
    if (ind < 1) || (ind > length(ca)) || (ind ~= floor(ind))
        error('Invalid index in experiment catalog');
    end

    edef = ca{ind};
    fprintf('Experiment #%d (%s):\n',ind,edef.name);
    fprintf('\tClass: %s; index: %d; dens: %d\n',edef.class,edef.index,edef.dens);
    fprintf('\tdist: %s; server: %s; client: %s; net: %s\n',edef.dist,edef.serversw,edef.clientsw,edef.netw);

end