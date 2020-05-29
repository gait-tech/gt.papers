% Run experiment for all NeuRA data
% NS2: use yaw offset from Vicon
%
% Example:
%       writetable(results, 'neura-sparse01/explore/results20190812.csv')
%
%% configuration
usebuffer = false; % use buffer in dir/mat if available
append_results = false; % append results if existing in buffer
export_c3d = false;

%% folder initialization
dir = 'data/neura-sparse01';
outDir = sprintf('%s/output', dir);
if ~exist(outDir, 'dir'), mkdir(outDir); end
matDir = sprintf('%s/mat', dir);
stepParamDir = '+papers/+lgcekf2019/+stepParam';
meas3DistDir = '';

addpath('mod-lib');
addpath('mod-lib/liese3lib');

%% load data list
dataList = readtable(sprintf('+papers/+lgcekf2019/data-list.csv', dir));
dataN = size(dataList, 1);
DATARANGE = 1:dataN;
% DATARANGE = 31;

results = table();

ns = "NS2+";
setups = {
    struct('est', 'vsxsens')
    };

for sI = [ struct('algo', 'ckf', 'pI', 1, 'mI', 76, 'cI', 355), ...
           struct('algo', 'lieekf', 'pI', 101, 'mI', 125, 'cI', 7) ] 
    for sdI = {'av03'}
        for initI = {'w__v'} % 
            setups{end+1} = struct('est', sI.algo, ...
                       'accData', 'w__s', 'oriData', 'w__s', 'accDataNoise', 0, ...
                       'initSrc', initI, 'stepDetection', sdI, ...
                       'applyPred', sI.pI, 'applyMeas', sI.mI, ...
                       'applyCstr', sI.cI, 'P', 0.5, ...
                       'sigma2QAcc', 1e1^2, ...
                       'c3d', export_c3d );
            setups{end+1} = struct('est', sI.algo, ...
                       'accData', initI, 'oriData', initI, 'accDataNoise', 0, ...
                       'initSrc', initI, 'stepDetection', sdI, ...
                       'applyPred', sI.pI, 'applyMeas', sI.mI, ...
                       'applyCstr', sI.cI, 'P', 0.5, ...
                       'sigma2QAcc', 1e1^2, ...
                       'c3d', export_c3d);
        end
    end
end


for i = 1:length(setups)
    setups{i}.label = getLabel(setups{i}, ns);
end

for i = DATARANGE
    n = table2struct(dataList(i, :));
    name = sprintf("%s%s-%s", ns, n.subj, n.act);
    fprintf("Data %3d/%3d: %s\n", i, dataN, name);
    
    % load data
    data = mocapdb.loadNeuRaTrialData(dir, n.subj, n.act, usebuffer);
    % preprocess the data
    [imuData, gfrAcc, step] = papers.lgcekf2019.preprocessing( ...
                    data.dataS, data.dataV, ...
                    data.calibYawFix, data.revStepDetect, ...
                    sprintf("%s/%s-preproc.mat", outDir, name));
    % run experiment for each configuration listed in setups
    actBody = data.dataV.togrBody(1:data.dataV.nSamples, ...
        {'name', name, 'oriUnit', 'deg', 'lnSymbol', '-', 'ptSymbol', '*', ...
         'fs', data.dataV.fs, 'xyzColor', {'m', 'y', 'c'}});
    
    r = papers.lgcekf2019.runExperiment(imuData, gfrAcc, step, ...
        actBody, data.dataX, setups, outDir, stepParamDir, meas3DistDir);
    results = [results; struct2table(r)];
end

% Append new results
if append_results
    dataPath = sprintf("%s/results.mat", outDir);
    if exist(dataPath, 'file')
        newResults = results;
        load(dataPath);
        [C, ia, ib] = intersect(results(:,{'name', 'label'}), newResults(:,{'name', 'label'}));
        results(ia,:) = [];
        results = [results; newResults];
    end
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