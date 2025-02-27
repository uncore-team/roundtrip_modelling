function [exitbranch,newstate] = AlgBernoulli(parms,state,rt)
% Take a new value from an scenario that is being modelled with Bernoulli.
% This is the analogue to onlineransac but using a simpler Bernoulli
% modelling process.
%
% PARMS -> a struct with the parameters of the algorithm:
%           .s -> minimum number of round-trip times to start modelling.
%           .win -> max number of round-trip times to have at any
%                   iteration.
% STATE -> a struct with the current state of the algorithm:
%           .model -> the model currently active (see ModelCreate) or empty
%                     if none. The model is of 'BERN' type.
%           .ind0 -> index within the global scenario of the first 
%                    round-trip time in the sample used
%                    so far.
%           .ind1 -> index of the last round-trip time in the sample used
%                    so far, within the complete scenario.
%          Initially, this struct must contain an empty model, and 
%          .ind* == NaN.
% RT -> the current round-trip time, to be processed in this call.
%
% EXITBRANCH <- reason of exit:
%               '1' <- no model yet
%               '2' <- model updated
% NEWSTATE <- updated STATE after processing the round-trip time RT.

    % --- sanity checks

    if (parms.s <= 1) || (parms.win < parms.s)
        error('Invalid params');
    end

    % --- preparation of data

    if isnan(state.ind0)
        ind0 = 1;
        ind1 = 1;
    else
        ind0 = state.ind0;
        ind1 = state.ind1 + 1;
    end
    if ind1 - ind0 + 1 > parms.win
        ind0 = ind0 + 1;
    end
    
    % --- core of the algorithm

    newstate = state;    
    newstate.ind0 = ind0;
    newstate.ind1 = ind1;
    newstate.model = ModelCreate('BERN');
    if ind1 - ind0 + 1 >= parms.s
        newstate.model.coeffs.ind0 = ind0;
        newstate.model.coeffs.ind1 = ind1;
        newstate.model.defined = 1;
        exitbranch = 2;
    else
        newstate.model = [];
        exitbranch = 1;
    end

end