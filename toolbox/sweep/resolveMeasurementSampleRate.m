function sampleRateHz = resolveMeasurementSampleRate(measurement)
% resolveMeasurementSampleRate - Resolve the two supported measurement rates.

hasFs = isfield(measurement, 'fs');
hasInfoFs = isfield(measurement, 'info_signal') && ...
    isstruct(measurement.info_signal) && ...
    isfield(measurement.info_signal, 'fsovs');

if hasFs
    sampleRateHz = double(measurement.fs);
    if hasInfoFs
        infoSampleRate = double(measurement.info_signal.fsovs);
        tolerance = 1e-9*max([1, abs(sampleRateHz), abs(infoSampleRate)]);
        if abs(infoSampleRate - sampleRateHz) > tolerance
            error('run_selected_comparison:SampleRateMismatch', ...
                'Measurement fs and info_signal.fsovs disagree.');
        end
    end
elseif hasInfoFs
    sampleRateHz = double(measurement.info_signal.fsovs);
else
    error('run_selected_comparison:MissingSampleRate', ...
        'Measurement must contain fs or info_signal.fsovs.');
end
end
