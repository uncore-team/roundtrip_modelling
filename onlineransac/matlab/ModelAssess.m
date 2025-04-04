function [m,g] = ModelAssess(S,ind0,ind1,models)
% If MODELS is a cell with some model types, select the first one that 
% fits ok to the data in S and is assessed with GoF; if models
% is a model, just assesses it with GoF.
%
% S -> sample to fit/gof.
% IND0, IND1 -> first and last indexes of the data within the complete 
%               scenario.
% MODELS -> either a cell with as many entries as model types, each one a 
%           text, or a model already defined (see ModelCreate).
%
% M <- model fitted/gofed to the data or a copy of MODELS if
%      MODELS is a model. Empty in any case the assessment fails.
% G <- empty if assessment fails; otherwise a struct with info about the 
%      validity of the gof:
%       .stat <- statistic of the GoF test.
%       .thresh <- threshold used for testing the statistic value.

    m = [];
    g = [];
    if iscell(models) % --- a list of possible models to find for these data
        nm = length(models);
        for f = 1:nm
            mo = ModelFit(S,ind0,ind1,models{f});
            if mo.defined
                [reject,stat,thresh] = ModelGof(mo,S,0); % model comes from data
                if ~reject
                    g = struct('stat',stat,'thresh',thresh);
                    m = mo;
                    return;
                end
            end
        end
    elseif isstruct(models) % --- just one model, previously deduced, not from these data
        if ~models.defined
            error('Undefined model to assess');
        end
        [reject,stat,thresh] = ModelGof(models,S,1); % model does not come from data
        if ~reject
            g = struct('stat',stat,'thresh',thresh);
            m = models;
            return;
        end        
    else
        error('Invalid type of MODELS parms');
    end

end