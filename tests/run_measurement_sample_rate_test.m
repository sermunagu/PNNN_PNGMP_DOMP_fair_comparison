% Verify the two supported measurement sample-rate sources and mismatch error.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'toolbox', 'sweep'));

expectedRate = 491.52e6;
sampleRateHz = resolveMeasurementSampleRate(struct('fs', expectedRate));
assert(sampleRateHz == expectedRate);

sampleRateHz = resolveMeasurementSampleRate( ...
    struct('info_signal', struct('fsovs', expectedRate)));
assert(sampleRateHz == expectedRate);

sampleRateHz = resolveMeasurementSampleRate( ...
    struct('fs', expectedRate, ...
    'info_signal', struct('fsovs', expectedRate)));
assert(sampleRateHz == expectedRate);

mismatchRaised = false;
try
    resolveMeasurementSampleRate(struct('fs', expectedRate, ...
        'info_signal', struct('fsovs', expectedRate + 1)));
catch exception
    mismatchRaised = exception.identifier == ...
        "run_selected_comparison:SampleRateMismatch";
end
assert(mismatchRaised);

fprintf('MEASUREMENT SAMPLE RATE TEST: PASS\n');
