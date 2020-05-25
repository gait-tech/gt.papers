% Run experiment for all NeuRA data
% NS2: use yaw offset from Vicon
%
% Example:
%       writetable(results, 'neura-sparse01/explore/results20190812.csv')
%
%% configuration
usebuffer = false; % use buffer in dir/mat if available

%% folder initialization
dir = 'data/neura-sparse01';
outDir = sprintf('%s/output', dir);
if ~exist(outDir, 'dir'), mkdir(outDir); end
matDir = sprintf('%s/mat', dir);
if ~exist(matDir, 'dir'), mkdir(matDir); end
stepParamDir = '+papers/+ckf2019/+stepParam';
if ~exist(stepParamDir, 'dir'), mkdir(stepParamDir); end
addpath('mod-lib');

%% load data list
dataList = readtable('+papers/+ckf2019/data-list.csv');
dataN = size(dataList, 1);
DATARANGE = 1:dataN;
% DATARANGE = 15;

results = table();

ns = "NS2-";
setups = { struct('est', 'ckf', ...
               'accData', 'w__v', 'oriData', 'w__v', 'accDataNoise', 0, ...
               'initSrc', 'w__v', 'stepDetection', 'av03', ...
               'applyPred', 1, 'applyMeas', 76, 'applyCstr', 355, 'P', 0.5, ...
               'sigma2QAcc', 1e1^2, 'c3d', true), ...
          struct('est', 'ckf', ...
               'accData', 'w__s', 'oriData', 'w__s', 'accDataNoise', 0, ...
               'initSrc', 'w__v', 'stepDetection', 'av03', ...
               'applyPred', 1, 'applyMeas', 76, 'applyCstr', 355, 'P', 0.5, ...
               'sigma2QAcc', 1e1^2, 'c3d', true), ...
          struct('est', 'vsxsens', 'c3d', true) };

for i = 1:length(setups)
    setups{i}.label = getLabel(setups{i}, ns);
end

for i = DATARANGE
    n = table2struct(dataList(i, :));
    name = sprintf("%s%s-%s", ns, n.subj, n.act);
    fprintf("Data %3d/%3d: %s\n", i, dataN, name);
    
    % load data
    data = mocapdb.loadNeuRaTrialData(dir, n.subj, n.act, usebuffer, ns);
    % preprocess the data
    [imuData, gfrAcc, step] = papers.ckf2019.preprocessing( ...
                    data.dataS, data.dataV, ...
                    data.calibYawFix, data.revStepDetect, ...
                    sprintf("%s/%s-preproc.mat", dir, name));
    % run experiment for each configuration listed in setups
    actBody = data.dataV.togrBody(1:data.dataV.nSamples, ...
        {'name', name, 'oriUnit', 'deg', 'lnSymbol', '-', 'ptSymbol', '*', ...
         'fs', data.dataV.fs, 'xyzColor', {'m', 'y', 'c'}});
    
    r = papers.ckf2019.runExperiment(imuData, gfrAcc, step, ...
        actBody, data.dataX, setups, outDir, stepParamDir);
    results = [results; struct2table(r)];
end

% Append new results
dataPath = sprintf("%s/results.mat", outDir);
if exist(dataPath, 'file')
    newResults = results;
    load(dataPath);
    [C, ia, ib] = intersect(results(:,{'name', 'label'}), newResults(:,{'name', 'label'}));
    results(ia,:) = [];
    results = [results; newResults];
end
save(sprintf("%s/results.mat", outDir), 'results')

function label = getLabel(setup, prefix)    
    if strcmp(setup.est, 'vsxsens')
        label = sprintf("%sviconvsxsens", prefix);
    else
        label = sprintf("%s%s+A%sO%sI%s+S%s+P%03d+M%03d+C%03d", prefix, setup.est, ...
                setup.accData, setup.oriData, setup.initSrc, ...
                setup.stepDetection, setup.applyPred, ...
                setup.applyMeas, setup.applyCstr);
    end
end