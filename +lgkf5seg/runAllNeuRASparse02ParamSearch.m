% Run experiment for all NeuRA data
% NS2: use yaw offset from Vicon
%
% Example:
%       writetable(results, 'neura-sparse01/explore/results20190812.csv')
%
%% configuration
usebuffer = false; % use buffer in dir/mat if available
append_results = true; % append results if existing in buffer
export_c3d = true;

%% folder initialization
dir = 'data/neura-sparse01';
outDir = sprintf('%s/output', dir);
if ~exist(outDir, 'dir'), mkdir(outDir); end
matDir = sprintf('%s/mat', dir);
stepParamDir = '+papers/+lgkf5seg/+stepParam';
meas3DistDir = '+papers/+ckfdist2020/+dist_sim';

addpath('mod-lib');
addpath('mod-lib/liese3lib');

%% load data list
dataList = readtable(sprintf('+papers/+lgkf5seg/data-list-train.csv', dir));
dataN = size(dataList, 1);
% DATARANGE = 1:dataN;
DATARANGE = 28;

results = table();

ns = "NS2-";
setups = {
%     struct('est', 'vsxsens')
    };

for sI = [ struct('algo', 'lieekfdist', 'pI', 00125, 'mI', 4125, 'cI', 7), ...
           struct('algo', 'lieekfdist', 'pI', 10125, 'mI', 4125, 'cI', 7), ...
           struct('algo', 'lieekfdist', 'pI', 20125, 'mI', 4125, 'cI', 7), ...
           struct('algo', 'lieekfdist', 'pI', 30125, 'mI', 4125, 'cI', 7), ...
           struct('algo', 'lieekfdist', 'pI', 40125, 'mI', 4125, 'cI', 7), ...
           struct('algo', 'lieekfdist', 'pI', 00121, 'mI', 4125, 'cI', 7)
           ] 
    for j = 2
        predI = sI.pI + j*1e3;
        setups{end+1} = struct('est', sI.algo, ...
                   'accData', 'w__s', 'oriData', 'w__s', 'accDataNoise', 0, ...
                   'initSrc', 'w__v', 'stepDetection', 'av03', ...
                   'applyPred', predI, 'applyMeas', sI.mI, ...
                   'applyCstr', sI.cI, 'P', 0.5, ...
                   'sigma2QAcc', 1e2, ...
                   'sigma2RLim', 10.0^mod(idivide(int32(predI), 100, 'floor'), 10), ...
                   'sigma2RXYPosPVLSRS', 10.0^mod(idivide(int32(predI), 10, 'floor'), 10), ...
                   'sigma2ROri', 10.0^(double(mod(idivide(int32(predI), 1000, 'floor'),10))-2), ...
                   'sigma2QGyr', 10.0^double(mod(idivide(int32(predI), 10000, 'floor'),10)), ...
                   'c3d', export_c3d );
    end
end

for sI = [ struct('algo', 'lieekfdist', 'pI', 00121, 'mI', 4345, 'cI', 7) , ... 
           struct('algo', 'lieekfdist', 'pI', 00125, 'mI', 4345, 'cI', 7) , ...
           struct('algo', 'lieekfdist', 'pI', 10125, 'mI', 4345, 'cI', 7) , ...
           struct('algo', 'lieekfdist', 'pI', 20125, 'mI', 4345, 'cI', 7) , ...
           struct('algo', 'lieekfdist', 'pI', 30125, 'mI', 4345, 'cI', 7) , ...
           ] 
    for j = 2
        predI = sI.pI + j*1e3;
        for jLegDist = 6
            for jAnkleDist = 5
                for s3dSigma = 0.1 % [0.01:0.01:0.09 0.15:0.05:0.2] % [0.01:0.01:0.1 0.15:0.05:0.2]
                    for s3dI = 1
                        s3dI2 = s3dI + jLegDist*1e1 + jAnkleDist*1e2;
                        setups{end+1} = struct('est', sI.algo, ...
                                   'accData', 'w__s', 'oriData', 'w__s', 'accDataNoise', 0, ...
                                   'initSrc', 'w__v', 'stepDetection', 'av03', ...
                                   'applyPred', predI, 'applyMeas', sI.mI, ...
                                   'applyCstr', sI.cI, 'P', 0.5, 'sigma2QAcc', 1e2, ...
                                   'sigma2RLim', 10.0^mod(idivide(int32(predI), 100, 'floor'), 10), ...
                                   'sigma2RXYPosPVLSRS', 10.0^mod(idivide(int32(predI), 10, 'floor'), 10), ...
                                   'sigma2ROri', 10.0^(double(mod(idivide(int32(predI), 1000, 'floor'),10))-2), ...
                                   'sigma2QGyr', 10.0^double(mod(idivide(int32(predI), 10000, 'floor'),10)), ...
                                   'sigma2DistLeg', 10.0^(double(mod(idivide(int32(s3dI2), 1e1, 'floor'),10))-5), ...
                                   'sigma2DistLARA', 10.0^(double(mod(idivide(int32(s3dI2), 1e2, 'floor'),10))-5), ...
                                   'c3d', export_c3d, ...
                                   'sim3Dist', s3dI2, 'sim3DistSigma', s3dSigma);
                    end
                end
            end
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
    [imuData, gfrAcc, step] = papers.lgkf5seg.preprocessing( ...
                    data.dataS, data.dataV, ...
                    data.calibYawFix, data.revStepDetect, ...
                    sprintf("%s/%s-preproc.mat", outDir, name));
    % run experiment for each configuration listed in setups
    actBody = data.dataV.togrBody(1:data.dataV.nSamples, ...
        {'name', name, 'oriUnit', 'deg', 'lnSymbol', '-', 'ptSymbol', '*', ...
         'fs', data.dataV.fs, 'xyzColor', {'m', 'y', 'c'}});
    
    r = papers.lgkf5seg.runExperiment(imuData, gfrAcc, step, ...
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
        label = sprintf("%s%s+A%sO%sI%s+S%s+P%05d+M%05d+C%03d", prefix, setup.est, ...
                setup.accData, setup.oriData, setup.initSrc, ...
                setup.stepDetection, setup.applyPred, ...
                setup.applyMeas, setup.applyCstr);
    end
end