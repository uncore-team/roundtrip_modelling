function regs = regimeDetection(S,onlineransacparms,withtrace)
% Perform detection of regimes in a given scenario.
%
% S -> sequence of roundtrip times.
% ONLINERANSACPARMS -> parameters for the online RANSAC algorithm (see
%                      AlgOnlineransac.m).
%
% REGS <- matrix with the detected regimes: one row per regime with these
%         columns: index of first roundtrip time, index of last roundtrip
%         time, model_as_coeffs (see ModelToCoeffs) for the regime
%         Notice that regimes may overlap, for instance if the parameters
%         include sample sliding.

    regs = [];
    n = length(S);
    if n > 0

        if withtrace
            fprintf(sprintf('Detecting regimes in %d data:\r\n',n));
        end
        ransacstate = struct('sample',[],'ind0',NaN,'ind1',NaN,...
                             'model',[],'consensus',[]);
        oldexitbranch = 1; % no model yet
        oldperc = 0;
        tic;
        for f = 1:n

            perc = f / n * 100;
            if (round(perc) > round(oldperc))
                oldperc = perc;
                if withtrace
                    fprintf('%.2f%%, f=%d, ',perc,f);
                    toc
                end
            end

            [exitbranch,ransacstate] = AlgOnlineransac(onlineransacparms,...
                                                    ransacstate,...
                                                    S(f));
            if (oldexitbranch == 1) && (exitbranch ~=1) % got a regime
                if f <= length(ransacstate.sample) - 1
                    error('Invalid starting index for a regime');
                end
                initreg = f - length(ransacstate.sample) + 1; % started in the past
                if initreg >= f
                    error('overlaping regimes');
                end
                model = ransacstate.model;
            elseif (oldexitbranch ~= 1) && (exitbranch == 1) % lost regime (current rtt is not in it)
                mps = ModelToCoeffs(model);
                regs = [regs ; initreg, f-1, mps];
            end

            oldexitbranch = exitbranch;
        end

        % complete trailing regime, if any
        if (oldexitbranch ~=1) % got a trailing regime
            mps = ModelToCoeffs(model);
            regs = [regs ; initreg, n, mps];
        end

    end

end