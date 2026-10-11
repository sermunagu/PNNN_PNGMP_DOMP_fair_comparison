function run_dpd_online_fresh_session_test(folder)
% Run in a new MATLAB process after run_dpd_offline_online_test(folder).
root = fileparts(fileparts(mfilename('fullpath')));
addpath(root);
fixture = load(fullfile(folder, 'fresh_session_fixture.mat'));
profile clear; profile on;
[outputFile, dpd] = runDPDOnline(fixture.modelA, fixture.inputFile);
profile off; info = profile('info');
for name = ["selectDOMPSupport", "prepare_pnnn_dense_source", "pruneAndFineTunePNNN", "trainnet"]
    assert(~any(strcmp({info.FunctionTable.FunctionName}, name)));
end
assert(isequal(dpd, fixture.expected));
saved = load(outputFile);
assert(isequal(saved.dpd, dpd));
fprintf('DPD fresh MATLAB session: PASS.\n');
end
