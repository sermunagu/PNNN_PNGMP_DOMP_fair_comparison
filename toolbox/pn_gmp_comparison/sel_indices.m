function indices = sel_indices(xdpd, ydpd, perc)
%sel_indices Selects a segment of the signal for modeling purposes.
%
%   indices = sel_indices(xdpd, ydpd, perc) returns the index range of a 
%   segment from the output signal ydpd, centered around its maximum 
%   amplitude, for use in modeling (e.g., Volterra-based or DPD modeling).
%
%   INPUTS:
%       xdpd  - Input signal (not used in this function, but included for
%               interface consistency or future extension).
%       ydpd  - Output signal from which the modeling segment is extracted.
%       perc  - Fraction (in [0,1]) indicating the relative length of the 
%               segment to extract (e.g., 0.015 for 1.5% of the signal).
%
%   OUTPUT:
%       indices - Index range [indmodinf : indmodsup] corresponding to the
%                 selected segment of ydpd.

[~, indy] = max(abs(ydpd));
indmodinf = floor(length(ydpd)*perc*floor((indy/length(ydpd))/perc))+1;
indmodsup = ceil(length(ydpd)*perc*ceil((indy/length(ydpd))/perc))-1;

if(indmodsup>length(ydpd))
    indmodsup = length(ydpd);
    indmodinf = length(ydpd)-floor(perc*length(ydpd));
end

indices = indmodinf:indmodsup;