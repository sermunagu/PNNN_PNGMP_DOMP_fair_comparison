function v = nmse_db(ref, est)
% nmse_db - Compute time-domain NMSE in dB for complex signals.
%
% This helper reports normalized prediction error after outputs have been
% reconstructed in the complex domain.
%
% Inputs:
%   ref - Reference complex signal.
%   est - Estimated complex signal.
%
% Outputs:
%   v - NMSE value in dB.

v = nmseComplexDb(ref, est);
end
