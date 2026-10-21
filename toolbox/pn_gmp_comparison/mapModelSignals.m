function [x_in, y_out] = mapModelSignals(x, y, mappingMode)
% Select local modeled-block input/output, not the physical PA direction.

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
