% Run experiment for all NeuRA data
%
% Example:
%       writetable(results, 'neura-sparse01/explore/results20190812.csv')
%
%% configuration
usebuffer = false; % use buffer in dir/mat if available
append_results = true; % append results if existing in buffer
export_c3d = true;

%% folder initialization
paperDir = '+papers/+lgkf7seg';
stepParamDir = sprintf('%s/+stepParam', paperDir);
% meas3DistDir = ''; % sprintf('%s/+dist_sim', paperDir);
meas3DistDir = '+papers/+ckfdist2020/+dist_sim';
addpath('mod-lib');
addpath('mod-lib/liese3lib');

outDir = 'data/neura-sparse01/output';

%% load data list            
dataList = [];
% neura data list
% for i = ["NS2", "NS5", "NS21", "NS51", "NS52", "NS53"]
for i = ["NS21", "NS51"]
    dataListNeura = readtable(sprintf('%s/neura-data-list.csv', paperDir));
    dataListNeura.db(:) = i;
    dataList = [dataList; dataListNeura];
end

% tcd data list
dataListTCD = readtable(sprintf('%s/tcd-data-list.csv', paperDir));
dataListTCD.startFrame(:) = -1;
dataListTCD.endFrame(:) = -1;
dataListTCD.db(:) = "tcd";
dataList = [dataList; dataListTCD];

dataN = size(dataList, 1);
DATARANGE = 1:dataN;
% DATARANGE0 = [22, 23, 44, 67, 140, 141, 150];
% DATARANGE = [DATARANGE0, DATARANGE0+161];
% DATARANGE = [31];

results = table();

setupsN = {
    struct('est', 'vsxsens', 'c3d', export_c3d), ...
    struct('est', 'vsxsens7s', 'c3d', export_c3d)
    };
setupsT = {
    struct('est', 'vsxsens', 'c3d', export_c3d), ...
    struct('est', 'vsxsens7s', 'c3d', export_c3d)
    };

for sI = [ struct('algo', 'ckfdist', 'pI', 1, 'mI',  76, 'cI', 355), ... %CKF-3I
           struct('algo', 'lieekfdist', 'pI', 32125, 'mI', 4125, 'cI', 7), ... %L5S-3I
           struct('algo', 'lgckf7sb', 'pI', 1, 'mI', 433, 'cI', 3), ... % L7S-3I
           struct('algo', 'lgckf7stwoimu6b', 'pI', 1, 'mI', 9432, 'cI', 3) ... % L7S-2I
            ]
    buf = struct('est', sI.algo, ...
               'accData', 'w__s', 'oriData', 'w__s', 'accDataNoise', 0, ...
               'initSrc', 'w__v', 'stepDetection', 'av03', ...
               'applyPred', sI.pI, 'applyMeas', sI.mI, ...
               'applyCstr', sI.cI, 'P', 0.5, ...
               'sigma2QAcc', 1e1^2, ...
               'c3d', export_c3d, 'debug', true);
    setupsN{end+1} = buf;
    buf.initSrc = 'w__d';
    setupsT{end+1} = buf;
%     buf.stepDetection = 'av02';
%     setupsT{end+1} = buf;
end

for i = 1:length(setupsN)
    setupsN{i}.label = getLabel(setupsN{i}, '');
end
for i = 1:length(setupsT)
    setupsT{i}.label = getLabel(setupsT{i}, '');
end

for i = DATARANGE
    n = table2struct(dataList(i, :));
    name = sprintf("%s-%s-%s", n.db, n.subj, n.act);
    if ~exist(outDir, 'dir'), mkdir(outDir); end
    fprintf("Data %3d/%3d: %s\n", i, dataN, name);
    if strncmp(n.act, 'Trial-Static', 12)
        fprintf("Skip\n");
        continue
    end
    
    if strncmp(n.db, 'NS', 2)
        % load data
        dataDir = 'data/neura-sparse01';
        data = mocapdb.loadNeuRaTrialData(dataDir, n.subj, n.act, usebuffer);
        buf = char(n.db);
        if buf(end) == 'a'
            data.dataS = data.dataS.removeAccBias(data.bias);
            nsmode = buf(4:end-1);
        else
            nsmode = buf(4:end);
        end
        if isempty(nsmode), nsmode = 0;
        else, nsmode = str2double(nsmode); end
            
        if strncmp(n.db, 'NS5', 3)
            dataV = data.dataVKFM;
        else
            dataV = data.dataV;
        end
        if bitand(nsmode, 2)
            calibSW2W = data.calibYawFixKFM;
        else
            calibSW2W = data.calibYawFix;
        end
        if ~bitand(nsmode, 1)
            calibSW2W.Pelvis.ori = [1 0 0 0];
        end
        calibS2BD = [];
        setups = setupsN;
    elseif strcmp(n.db, 'tcd')
        dataDir = 'data/totalcapture';
        data = mocapdb.loadTCDTrialData(dataDir, n.subj, n.act, usebuffer);
        dataV = data.dataV;
        calibSW2W = data.calibW2V.adjustFrame(data.calibV2W, [1,0,0,0], true);
        calibS2BD = data.calibS2B;
        setups = setupsT;
    else
        error('runAllSparse db %s not supported.', n.db);
    end
    
    % run experiment for each configuration listed in setups
    actBody = dataV.togrBody(1:dataV.nSamples, ...
        {'name', name, 'oriUnit', 'deg', 'lnSymbol', '-', 'ptSymbol', '*', ...
         'fs', dataV.fs, 'xyzColor', {'m', 'y', 'c'}});
    B_p_TP = struct('LFT', [0 0 0.5*actBody.calcLFootLength() 1], ...
                    'RFT', [0 0 0.5*actBody.calcRFootLength() 1] );
    for j = 1:length(setups)
        setups{j}.LFA_p_LF = B_p_TP.LFT; 
        setups{j}.RFA_p_RF = B_p_TP.RFT;
    end
    
    % preprocess the data
    [imuData, gfrAcc, step] = papers.lgkf7seg.preprocessing( ...
                    data.dataS, dataV, B_p_TP, ...
                    calibSW2W, calibS2BD, data.revStepDetect, ...
                    sprintf("%s/%s-preproc.mat", outDir, name));
                    
    r = papers.lgkf7seg.runExperiment(imuData, gfrAcc, step, ...
        actBody, data.dataX, setups, outDir, stepParamDir, meas3DistDir);
    r = struct2table(r);
    r.db(:) = n.db;
    results = [results; r];
end

% Append new results
dataPath = sprintf("%s/results.mat", outDir);
if append_results && exist(dataPath, 'file')
    newResults = results;
    load(dataPath);
    [C, ia, ib] = intersect(results(:,{'name', 'label'}), newResults(:,{'name', 'label'}));
    results(ia,:) = [];
    results = [results; newResults];
end
save(dataPath, 'results')

function label = getLabel(setup, prefix)    
    if strcmp(setup.est, 'vsxsens')
        label = sprintf("%sviconvsxsens", prefix);
    elseif strcmp(setup.est, 'vsxsens7s')
        label = sprintf("%sviconvsxsens7s", prefix);
    else
        label = sprintf("%s%s+A%sO%sI%s+S%s+P%03d+M%03d+C%03d", prefix, setup.est, ...
                setup.accData, setup.oriData, setup.initSrc, ...
                setup.stepDetection, setup.applyPred, ...
                setup.applyMeas, setup.applyCstr);
    end
end