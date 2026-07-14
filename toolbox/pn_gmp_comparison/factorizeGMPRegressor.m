function descriptor = factorizeGMPRegressor(regressor, index)
% factorizeGMPRegressor - Describe one GMP column from X/Xconj/Xenv metadata.
% Canonical carrier-envelope terms are separated from the two legacy
% auxiliaries so explicit PN-IQ features never rely on population ordering.

if nargin < 2 || isempty(index)
    index = NaN;
end
if ~isobject(regressor) || ~isprop(regressor, 'X') || ...
        ~isprop(regressor, 'Xconj') || ~isprop(regressor, 'Xenv')
    error('factorizeGMPRegressor:InvalidRegressor', ...
        'regressor must expose X, Xconj, and Xenv properties.');
end

x_terms = double(regressor.X(:).');
x_conjugate_terms = double(regressor.Xconj(:).');
envelope_terms = double(regressor.Xenv(:).');
[envelope_lags, envelope_powers] = groupTerms(envelope_terms);

descriptor = struct();
descriptor.index = double(index);
descriptor.originalExpression = string(regressor.print());
descriptor.family = "unsupported/noncanonical term";
descriptor.canonicalGMP = false;
descriptor.auxiliaryType = "none";
descriptor.carrierLag = NaN;
descriptor.carrierType = "none";
descriptor.envelopeLags = envelope_lags;
descriptor.envelopePowers = envelope_powers;
descriptor.envelopeDescription = describeEnvelope( ...
    envelope_lags, envelope_powers);
descriptor.structurallyRealAfterPhaseNormalization = false;
descriptor.IColumnStructurallyZero = false;
descriptor.QColumnStructurallyZero = false;
descriptor.iSignature = "";
descriptor.qSignature = "";
descriptor.magnitudeSignatures = strings(0, 1);
descriptor.notes = "Exact normalized I/Q fallback required.";

if numel(x_terms) == 1 && isempty(x_conjugate_terms)
    descriptor.family = "canonical GMP carrier-envelope term";
    descriptor.canonicalGMP = true;
    descriptor.carrierLag = x_terms(1);
    descriptor.carrierType = "direct complex carrier";
    descriptor.structurallyRealAfterPhaseNormalization = ...
        x_terms(1) == 0;
    descriptor.QColumnStructurallyZero = x_terms(1) == 0;
    base_signature = carrierEnvelopeSignature( ...
        x_terms(1), envelope_lags, envelope_powers);
    descriptor.iSignature = "IQ:I:" + base_signature;
    descriptor.qSignature = "IQ:Q:" + base_signature;
    [magnitude_lags, magnitude_powers] = addMagnitudeFactor( ...
        envelope_lags, envelope_powers, x_terms(1));
    descriptor.magnitudeSignatures = "MAG:" + ...
        envelopeSignature(magnitude_lags, magnitude_powers);
    descriptor.notes = "Explicit carrier I/Q and envelope product.";
elseif isempty(x_terms) && isequal(x_conjugate_terms, 0) && ...
        isempty(envelope_terms)
    descriptor.family = "conjugate auxiliary term";
    descriptor.auxiliaryType = "conjugate";
    descriptor.carrierLag = 0;
    descriptor.carrierType = "conjugate auxiliary";
    descriptor.iSignature = "AUX_CONJ:I";
    descriptor.qSignature = "AUX_CONJ:Q";
    descriptor.magnitudeSignatures = ["AUX_CONJ:I"; "AUX_CONJ:Q"];
    descriptor.notes = ["Noncanonical auxiliary retained as exact " + ...
        "phase-normalized I/Q fallback features."];
elseif isempty(x_terms) && isempty(x_conjugate_terms) && ...
        isequal(envelope_terms, 0)
    descriptor.family = "envelope-only auxiliary term";
    descriptor.auxiliaryType = "envelope";
    descriptor.carrierType = "none (envelope-only)";
    descriptor.iSignature = "AUX_ENV:I";
    descriptor.qSignature = "AUX_ENV:Q";
    descriptor.magnitudeSignatures = ["AUX_ENV:I"; "AUX_ENV:Q"];
    descriptor.notes = ["Noncanonical auxiliary retained as exact " + ...
        "phase-normalized I/Q fallback features."];
else
    suffix = compose('%d', index);
    descriptor.iSignature = "UNSUPPORTED_" + suffix + ":I";
    descriptor.qSignature = "UNSUPPORTED_" + suffix + ":Q";
    descriptor.magnitudeSignatures = [ ...
        "UNSUPPORTED_" + suffix + ":I"; ...
        "UNSUPPORTED_" + suffix + ":Q"];
end
end

function [lags, powers] = groupTerms(terms)
if isempty(terms)
    lags = zeros(1, 0);
    powers = zeros(1, 0);
    return;
end
lags = unique(terms, 'sorted');
powers = zeros(size(lags));
for k = 1:numel(lags)
    powers(k) = nnz(terms == lags(k));
end
end

function [lags, powers] = addMagnitudeFactor(lags, powers, carrier_lag)
position = find(lags == carrier_lag, 1);
if isempty(position)
    lags(end+1) = carrier_lag;
    powers(end+1) = 1;
else
    powers(position) = powers(position) + 1;
end
[lags, order] = sort(lags);
powers = powers(order);
end

function signature = carrierEnvelopeSignature(carrier_lag, lags, powers)
signature = "c" + compose('%+d', carrier_lag) + ":" + ...
    envelopeSignature(lags, powers);
end

function signature = envelopeSignature(lags, powers)
if isempty(lags)
    signature = "e1";
    return;
end
parts = strings(numel(lags), 1);
for k = 1:numel(lags)
    parts(k) = "a" + compose('%+d', lags(k)) + ...
        "^" + compose('%d', powers(k));
end
signature = strjoin(parts, '*');
end

function description = describeEnvelope(lags, powers)
if isempty(lags)
    description = "1";
    return;
end
parts = strings(numel(lags), 1);
for k = 1:numel(lags)
    parts(k) = "|x(n" + compose('%+d', -lags(k)) + ...
        ")|^" + compose('%d', powers(k));
end
description = strjoin(parts, ' ');
end
