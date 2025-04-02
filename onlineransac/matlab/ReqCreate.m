function r = ReqCreate(prob,t0,t1)
% Create a soft real-time requirement.
%
% PROB -> probability of having roundtrip times within [t0,t1], in [0,1].
% T0,T1 -> bounds of the requirement, in the same time unit as the 
%          roundtrip samples and models.
%
% R <- resulting struct:
%       .p <- prob
%       .t0,.t1 <- bounds

    if (prob < 0) || (prob > 1) || (t0 < 0) || (t1 < t0)
        error('Invalid requirement');
    end
    r = struct('p',prob,'t0',t0,'t1',t1);

end