function [exitbranch,newstate] = AlgOnlineransac(parms,state,rt)
% Do one step of the on-line RANSAC algorithm for modelling round-trip
% times.
%
% PARMS -> a struct with the parameters of the algorithm:
%           .s -> minimum number of round-trip times to start modelling.
%           .N -> maximum number of round-trip times to be modelled. 
%                 If NaN, no maximum.
%           .mtypes -> a cell with as many strings as model types we wish
%                      to use for modelling, in the same order we wish to
%                      try them at modelling. See ModelCreate().
%           .modelpreserving -> either 0 (no preserving), 1 (first
%                               preserve, then rebuild), or 2 (first
%                               rebuild, then preserve).
%           .samplesliding -> either 0 (forgetting) or 1 (sliding).
%           .datapreserving -> either 0 (forgetting) or 1 (preserving).
% STATE -> a struct with the current state of the algorithm:
%           .sample -> current round-trip sample used for modelling.
%           .ind0 -> index within the global scenario of the first 
%                    round-trip time in the sample used
%                    so far.
%           .ind1 -> index of the last round-trip time in the sample used
%                    so far, within the complete scenario.
%           .model -> the model currently active (see ModelCreate) or empty
%                     if none. The model has been assessed for the SAMPLE.
%           .consensus -> info about the validity of the consensus set (see
%                         ModelGof) or empty if none.
%          Initially, this struct must contain an empty sample, model and
%          consensus, and .ind* == NaN.
% RT -> the current round-trip time, to be processed in this call.
%
% EXITBRANCH <- reason of exit:
%               '1' <- no model yet
%               '2' <- model updated
%               '3' <- model preserved
% NEWSTATE <- updated STATE after processing the round-trip time RT.

    % --- sanity checks

    if (parms.s <= 1)
        error('Invalid params.s');
    end
    if (length(parms.mtypes) < 1)
        error('No model specified');
    end

    % --- preparation of data

    S = [state.sample, rt];
    if isnan(state.ind0)
        ind0 = 1;
        ind1 = 1;
    else
        ind0 = state.ind0;
        ind1 = state.ind0 + length(state.sample);
    end
    ns = length(S);
    newstate = state;
    
    % --- core of the algorithm

    if (ns < parms.s)

        exitbranch = 1;

    elseif (ns == parms.s)

        [m,g] = ModelAssess(S,ind0,ind1,parms.mtypes); % create new model
        if (~isempty(m)) && (~isempty(g))
            exitbranch = 2;
        else
            exitbranch = 1;
        end

    else
        
        switch parms.modelpreserving
            case 0

                [m,g] = ModelAssess(S,ind0,ind1,parms.mtypes); % create new model
                if (~isempty(m)) && (~isempty(g))
                    exitbranch = 2;
                else
                    if ~parms.samplesliding
                        S = rt;
                        ind0 = ind1;
                    end
                    exitbranch = 1;
                end

            case 1

                preserved = 0;
                if ~isempty(state.model)
                    % can fail the model coefficients if they are
                    % incompatible with the new sample: not the same as the
                    % sample on which the model was fitted
                    [~,g] = ModelAssess(S,ind0,ind1,state.model); % assess existing model
                    if ~isempty(g) 
                        newstate.consensus = g; % newstate.model is preserved
                        preserved = 1;
                        exitbranch = 3;
                    end
                end
                if ~preserved
                    [m,g] = ModelAssess(S,ind0,ind1,parms.mtypes); % create new model
                    if (~isempty(m)) && (~isempty(g))
                        exitbranch = 2;
                    else
                        if ~parms.samplesliding
                            S = rt;
                            ind0 = ind1;
                        end
                        exitbranch = 1;
                    end
                end

            case 2

                oldmodel = state.model;
                [m,g] = ModelAssess(S,ind0,ind1,parms.mtypes); % create new model
                if (~isempty(m)) && (~isempty(g))
                    exitbranch = 2;
                else
                    if ~isempty(oldmodel) % oldmodel == state.model
                        [~,g] = ModelAssess(S,ind0,ind1,state.model); % assess existing model
                        if ~isempty(g) 
                            newstate.consensus = g;
                            exitbranch = 3;
                        else
                            if ~parms.samplesliding
                                S = rt;
                                ind0 = ind1;
                            end
                            exitbranch = 1;
                        end
                    else
                        exitbranch = 1;
                    end
                end

            otherwise 
                error('Invalid model-preserving');
        end

    end

    % --- preparing new state

    switch exitbranch
        case 1 % no model yet
            ns = length(S); % may have varied due to sliding; ind0,ind1 updated already
            if (~parms.datapreserving) && (ns >= parms.s)  % parms.s > 2 here
                % preserve only the newest s-1 values
                newstate.sample = S(ns - (parms.s - 1) + 1 : ns);
                newstate.ind0 = ind1 - (parms.s - 1) + 1; % ind1 points to the last one in S right before this
                newstate.ind1 = ind1;
            else
                newstate.sample = S;
                newstate.ind0 = ind0;
                newstate.ind1 = ind1;
            end
            newstate.model = [];
            newstate.consensus = [];
        case 2 % model + consensus created now
            newstate.sample = S;
            newstate.ind0 = ind0;
            newstate.ind1 = ind1;
            newstate.model = m; % m, g calculated in the corresponding case
            newstate.consensus = g;
        case 3 % model preserved
            newstate.sample = S;
            newstate.ind0 = ind0;
            newstate.ind1 = ind1;
            % model and consensus already updated in the corresponding case
        otherwise
            error('Invalid exit branch');
    end

    % --- bounding gathered sample

    ns = length(newstate.sample);
    if (~isnan(parms.N)) && (ns > parms.N)
        newstate.sample = newstate.sample(ns - parms.N + 1 : ns);
        newstate.ind0 = newstate.ind1 - parms.N + 1;
        % if the bounded sample violates current consensus, that will be
        % detected in the next iteration by the core algorithm
    end
    if newstate.ind1 - newstate.ind0 + 1 ~= length(newstate.sample)
        error('Invalid indexes of sample');
    end

end

