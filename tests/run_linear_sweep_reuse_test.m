% Guard the one-path, shared-matrix structure of the linear sweep.
% This source-level fixture prevents target-loop DOMP or regressor rebuilding.
% It reads MATLAB code only and does not load data or fit any model.

clearvars;
projectRoot = fileparts(fileparts(mfilename('fullpath')));
linearFile = fullfile(projectRoot, 'toolbox', 'sweep', ...
    'runLinearComplexitySweep.m');
runnerFile = fullfile(projectRoot, 'run_parameter_sweep.m');
linearSource = fileread(linearFile);
runnerSource = string(fileread(runnerFile));

loopStart = strfind(linearSource, 'for index = 1:numel(targets)');
loopEnd = strfind(linearSource, 'sweep.complexTable = complexTable;');
assert(isscalar(loopStart) && isscalar(loopEnd) && loopStart < loopEnd);
targetLoop = string(linearSource(loopStart:loopEnd - 1));
for forbidden = ["buildGMPRegressorRows", "buildPNDomain", ...
        "buildPhaseNormalizedIQRegressors", "selectDOMPSupport", ...
        "selectSharedIQFeatures"]
    assert(~contains(targetLoop, forbidden));
end

assert(count(string(linearSource), "selectDOMPSupport(") == 2);
assert(count(string(linearSource), "selectSharedIQFeatures(") == 2);
assert(count(string(linearSource), "predictComplexGrid(") == 2);
assert(count(string(linearSource), "predictPNGrid(") == 2);
assert(count(runnerSource, "runLinearComplexitySweep(") == 1);

fprintf('LINEAR SWEEP REUSE TEST: PASS\n');
