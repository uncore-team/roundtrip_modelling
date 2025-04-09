function ConstantsInit()
% Define global constants.
% Call this function in any other that needs the constants.

    global TOLROUNDTRIPS % minimum significative increment in a roundtrip time

    if isempty(TOLROUNDTRIPS)
        TOLROUNDTRIPS = 1e-9;
    end

end