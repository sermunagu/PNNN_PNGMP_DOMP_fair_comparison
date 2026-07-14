function [x_in, y_out] = selectXYByMapping(x, y, mappingMode)
% selectXYByMapping - Select modeled-block input and output signals.
%
% This helper applies the local PNNN X/Y convention before dataset
% construction. X and Y are local modeled-block variables; mappingMode must
% not be interpreted automatically as PA-forward physical direction.
%
% Inputs:
%   x, y - Signals loaded from the measurement file.
%   mappingMode - Local modeled-block mapping mode.
%
% Outputs:
%   x_in - Input of the block being modeled.
%   y_out - Output of the block being modeled.

x = x(:);
y = y(:);
switch string(mappingMode)
    case "xy_forward"
        x_in = x;
        y_out = y;
    case "yx_inverse"
        x_in = y;
        y_out = x;
    otherwise
        error("mappingMode debe ser 'xy_forward' o 'yx_inverse'.");
end
end
