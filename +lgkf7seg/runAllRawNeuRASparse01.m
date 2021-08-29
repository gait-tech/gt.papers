% ======================================================================
%> Run experiment for all NeuRA data
%> NS1: use yaw offset from ROM
%> NS2: use yaw offset from Vicon
%> writetable(results, 'C:/Users/z5151460/OneDrive - UNSW/Thesis - Sparse Mocap/Goal03-LGCEKF7seg/results20200303.csv')
% ======================================================================
dir = 'data/neura-sparse01';
expDir = sprintf('%s/output', dir);
mkdir(expDir);
mkdir(sprintf('%s/mat', dir));
addpath('mod-lib');
addpath('mod-lib/liese3lib');

paperDir = '+papers/+lgkf7seg';
DEGRANGE = (0:0.1:359) - 180;
dataList = readtable(sprintf('%s/neura-data-list-all.csv', paperDir));
dataN = size(dataList, 1);
% DATARANGE = 1:dataN;
% DATARANGE = [15, 34];
DATARANGE = 34;

meas3DistLoc = sprintf('%s/+dist_sim', paperDir);
mkdir(meas3DistLoc);
stepParamDir = sprintf('%s/+stepParam', paperDir);
mkdir(stepParamDir);

options = struct('Pelvis', '00B40B91', ...
    'L_UpLeg', '00B40C45', 'R_UpLeg', '00B40C3C', ...
    'L_LowLeg', '00B40C44', 'R_LowLeg', '00B40C47', ...
    'L_Foot', '00B40C55', 'R_Foot', '00B40C48');

results = table();

nsList = {'NS2'};
for nsI = 1:length(nsList)
    ns = nsList{nsI};
    
    setups = {
    };
    
    for sI = [ struct('algo', 'ckfdist', 'pI', 1, 'mI', 76, 'cI', 355), ...
               struct('algo', 'lieekfv1', 'pI', 101, 'mI', 125, 'cI', 7), ...
               struct('algo', 'lgckf7sv2b', 'pI', 1, 'mI', 133, 'cI', 3), ...
               struct('algo', 'lgckf7sv2b', 'pI', 1, 'mI', 533, 'cI', 3), ...
               struct('algo', 'lgckf7sv2twoimu1b', 'pI', 1, 'mI', 533, 'cI', 3) ]
        for sdI = {'av03'}
            for initI = {'w__v'} % , 'w__m'} % {'w__x'}
%                 setups{end+1} = struct('est', sI.algo, ...
%                            'accData', initI, 'oriData', initI, 'accDataNoise', 0, ...
%                            'initSrc', initI, 'stepDetection', sdI, ...
%                            'applyPred', sI.pI, 'applyMeas', sI.mI, ...
%                            'applyCstr', sI.cI, 'P', 0.5, ...
%                            'sigmaQAcc', 1e1, 'sigmaQAngVel', 1e2, ...
%                            'c3d', false);
                       
                accIOptions = {'w__s'};
                smthZPosMode = '';
                
                for accI = accIOptions
                    setups{end+1} = struct('est', sI.algo, ...
                               'accData', accI, 'oriData', 'w__s', 'accDataNoise', 0, ...
                               'initSrc', initI, 'stepDetection', sdI, ...
                               'applyPred', sI.pI, 'applyMeas', sI.mI, ...
                               'applyCstr', sI.cI, 'P', 0.5, ...
                               'smthZPosMode', smthZPosMode, ...
                               'sigmaQAcc', 1e1, 'sigmaQAngVel', 1e2, ...
                               'c3d', false );
                end
            end
        end
    end

    for i = 1:length(setups)
        setups{i}.label = getLabel(ns, setups{i});
    end

    dataN = size(dataList, 1);

    for i = DATARANGE
%     for i = []
        n = table2struct(dataList(i, :));
        
%         uwbDistSigma = 0.0;
%         if strcmp(nxxxxxs(1:3), 'NS2') && (size(ns, 2) > 3)
%             uwbDistSigma = str2double(ns(5:end)); 
%         end
%         
%         if ~(strcmp(n.act(7:end-2), "SpeedSkater") || ...
%            strcmp(n.act(7:end-2), "Static") || ...
%            strcmp(n.act(7:end-2), "HighKneeJog")) && ...
%            uwbDistSigma <= (0.03+1e-4)
% %             setups{1}.sigmaUwbLeg = 1e-1;
%             setups{1}.sigmaUwbLeg = 1e0;
%         else
%             setups{1}.sigmaUwbLeg = 1e0;
%             continue;
%         end

        name = sprintf("%s-%s-%s", ns, n.subj, n.act);
        dataPath = sprintf('%s/mat/%s-%s-%s.mat', dir, ns(1:3), n.subj, n.act);
%         if exist(dataPath, 'file')
%             load(dataPath, 'data');
%         else
            data = struct('name', name, ...
                'fnameV', sprintf('%s/rawvicon/%s-%s.mat', dir, n.subj, n.act), ...
                'fnameX', sprintf('%s/xsens/%s-%s.bvh', dir, n.subj, n.act), ...
                'fnameS', sprintf('%s/rawimu/%s-%s', dir, n.subj, n.act), ...
                'calibFnameSensorYawFixWorldFrame', ...
                sprintf('%s/calib/%s-%s-Calib-SensorYawFixWorldFrame.txt', ...
                        dir, n.subj, n.act), ...
                'calibFnameSensorW2V', ...
                sprintf('%s/calib/%s-Calib-SensorW2V.txt', dir, n.subj), ...
                'fnameRevStepDetect', ...
                sprintf('%s/rawstep-detect/%s-%s-revStepDetect.csv', ...
                        dir, n.subj, n.act));

            data.dataV = mocapdb.ViconBody.loadViconMat(data.fnameV);           
            if exist(data.fnameX, 'file')
                data.dataX = mocapdb.BVHBody.loadXsensBVHFile(data.fnameX, "mm");
            else
                data.dataX = [];
            end
            
            if strcmp(n.subj, 'S03')
                nSamples = min(data.dataV.nSamples, data.dataX.nSamples);
                data.dataV = data.dataV.getSubset(1:nSamples);
                data.dataX = data.dataX.getSubset(1:nSamples);
                
                qXsensV2W = rotm2quat([0 0 1; 1 0 0; 0 1 0]);
                dataXtmp = data.dataX.toWorldFrame(qXsensV2W);
                data.dataV.qLTH = dataXtmp.qLeftUpLeg;
                data.dataV.qLSK = dataXtmp.qLeftLeg;
%                 sIdx = max(1, n.startFrame);
%                 qOffset = quatmultiply(dataXtmp.qLeftLeg(sIdx,:), ...
%                                 quatconj(data.dataV.qLSK(sIdx,:)));
%                 data.dataV.qLSK = quatmultiply(qOffset, data.dataV.qLSK);
%                 data.dataV.qLTH = quatmultiply(axang2quat([0 0 1 deg2rad(-50)]), data.dataV.qLTH);
%                 data.dataV.qLSK = quatmultiply(axang2quat([0 0 1 deg2rad(-50)]), data.dataV.qLSK);
                dLTH = vecnorm(data.dataV.LFEO-data.dataV.LFEP, 2, 2);
                dLSK = vecnorm(data.dataV.LFEO-data.dataV.LTIO, 2, 2);
                LTHz = quat2rotm(data.dataV.qLTH); 
                LTHz = squeeze(LTHz(:,3,:))';
                LSKz = quat2rotm(data.dataV.qLSK); 
                LSKz = squeeze(LSKz(:,3,:))';
                data.dataV.LFEO = data.dataV.LFEP - dLTH.*LTHz;
                data.dataV.LTIO = data.dataV.LFEO - dLSK.*LSKz;
            end
            
            data.dataS = mocapdb.XsensBody.loadMTExport(data.fnameS, options);
            data.dataS.fs = 100;
            data.calibV2W = rotm2quat(mocapdb.loadPendulumCompassMat( ...
                 sprintf('%s/calib/%s-Calib-V2W-Pendulum.mat', dir, n.subj), ...
                 sprintf('%s/calib/%s-Calib-V2W-Compass.mat', dir, n.subj))' );
            if strcmp(ns(1:3), 'NS1')
                if exist(data.calibFnameSensorYawFixWorldFrame, 'file')
                    data.calibYawFix = mocapdb.XsensBody.loadCalibCSV(data.calibFnameSensorYawFixWorldFrame);
                else
                    % using ROM calibration
                    dataS = mocapdb.XsensBody.loadMTExport(sprintf('%s/imu/%s-Trial-Walk-1', dir, n.subj), options);
                    dataS.fs = 100;           
                    W__dataV = mocapdb.ViconBody.loadViconMat(sprintf('%s/vicon/%s-Trial-Walk-1.mat', dir, n.subj));
                    W__dataV.changePosUnit('m', true);
                    W__dataV = W__dataV.toWorldFrame(data.calibV2W);

                    sIdx = max(W__dataV.getStartIndex()+1, 100);

                    viconCalibSB = dataS.calcCalibSB(W__dataV.togrBody(sIdx+1:sIdx+1, {}), sIdx(1));
                    data.calibYawFix = dataS.calcCalibAnkleSensorW2PelvisWFromROM(viconCalibSB, DEGRANGE);
                    data.calibYawFix.saveCalibCSV(data.calibFnameSensorYawFixWorldFrame);
                end
            elseif strcmp(ns(1:3), 'NS2')
                if exist(data.calibFnameSensorYawFixWorldFrame, 'file')
                    data.calibYawFix = mocapdb.XsensBody.loadCalibCSV(data.calibFnameSensorYawFixWorldFrame);
                else
                    W__dataV = data.dataV.toWorldFrame(data.calibV2W);
                    W__dataV.changePosUnit('m', true);
                    data.calibYawFix = data.dataS.calcCalibAnkleSensorW2PelvisWFromVicon(W__dataV);
                    data.calibYawFix.saveCalibCSV(data.calibFnameSensorYawFixWorldFrame);
                end
            elseif strcmp(ns(1:3), 'NS3')
                W__dataX = data.dataX;
                data.calibYawFix = data.dataS.calcCalibAnkleSensorW2PelvisWFromVicon(W__dataX);
            else
                data.calibYawFix = struct();
                for j = ["Pelvis", "L_LowLeg", "R_LowLeg", "L_Foot", "R_Foot"]
                    data.calibYawFix.(j).ori = [1 0 0 0];
                end
            end

            if exist(data.calibFnameSensorW2V, 'file')
                data.calibW2V = mocapdb.XsensBody.loadCalibCSV(data.calibFnameSensorW2V);
            else
                  % using trackingmount calibration
    %             data.calibW2V = mocapdb.XsensBody.loadCalibSensorW2V( ...
    %                  sprintf('%s/calib/%s-Calib-SensorW2V.mat', dir, n.subj), ...
    %                  sprintf('%s/calib/%s-Calib-SensorW2V', dir, n.subj), ...
    %                  options, 100);
                  % using ROM calibration
                dataS = mocapdb.XsensBody.loadMTExport(sprintf('%s/imu/%s-Trial-Walk-1', dir, n.subj), options);
                dataS.fs = 100;           
                V__dataV = mocapdb.ViconBody.loadViconMat(sprintf('%s/vicon/%s-Trial-Walk-1.mat', dir, n.subj));
                V__dataV.changePosUnit('m', true);

                sIdx = max(V__dataV.getStartIndex()+1, 100);

                viconCalibSB = dataS.calcCalibSB(V__dataV.togrBody(sIdx+1:sIdx+1, {}), sIdx(1));
                data.calibW2V = dataS.calcCalibAnkleSensorW2PelvisWFromROM(viconCalibSB, DEGRANGE);
                data.calibW2V.saveCalibCSV(data.calibFnameSensorW2V);
            end

%             if false
%                 % calculate acc pelvis bias from sitting trial.
%                 % pelvis bias results to crouch divergence so not using it for now.
%                 % not 100% tested but it works very well for the static trial
%                 dataS = mocapdb.XsensBody.loadMTExport(sprintf('%s/imu/%s-Trial-TUG-1', dir, n.subj), options);
%                 dataS.fs = 100;               
%                 dataV = mocapdb.ViconBody.loadViconMat(sprintf('%s/vicon/%s-Trial-TUG-1.mat', dir, n.subj));
%                 dataX = mocapdb.BVHBody.loadXsensBVHFile(sprintf('%s/xsens/%s-Trial-TUG-1.bvh', dir, n.subj), "mm");
%                 data.bias = experiment.calcNeuRAPelvisAccBias01(dataS, dataV, ...
%                                         data.calibV2W, data.calibW2V, dataX, ...
%                                         100, 1000);
%             elseif false
%                 % calculate acc pelvis bias from vicon data. Use static trial.
%                 % pelvis bias results to crouch divergence so not using it for now.
%                 % not 100% tested but it works very well for the static trial
%                 dataS = mocapdb.XsensBody.loadMTExport(sprintf('%s/imu/%s-Trial-Static-1', dir, n.subj), options);
%                 dataS.fs = 100;               
%                 dataV = mocapdb.ViconBody.loadViconMat(sprintf('%s/vicon/%s-Trial-Static-1.mat', dir, n.subj));
%                 dataX = mocapdb.BVHBody.loadXsensBVHFile(sprintf('%s/xsens/%s-Trial-Static-1.bvh', dir, n.subj), "mm");
%                 data.bias = experiment.calcNeuRAPelvisAccBias02(dataS, dataV, ...
%                                         data.calibV2W, data.calibW2V, dataX, ...
%                                         1, -1);
%                 data.bias
%             else
                data.bias = struct('w__v', zeros(1, 3), 'v__v', zeros(1, 3), ...
                          'w__x', zeros(1, 3));
%             end
            
            data.revStepDetect = readtable(data.fnameRevStepDetect);
            
            save(dataPath, 'data');
%         end
        
        data.name = name;
        
        fprintf("Data %3d/%3d: %s\n", i, dataN, data.name);
        r = papers.lgkf7seg.runRawNeuRASparse01Experiment(data.dataS, ...
                data.dataV, data.calibV2W, data.calibYawFix, data.calibW2V, [], ...
                data.dataX, data.revStepDetect, meas3DistLoc, ...
                data.name, setups, expDir, stepParamDir, ...
                n.startFrame, n.endFrame, data.bias);
        results = [results; struct2table(r)];
    end
end

%% file list vicon vs xsens comparison
rIdx = size(results, 1) + 1;
results = table2struct(results);

for i = DATARANGE
% for i = []
    n = table2struct(dataList(i, :));
    
    name = sprintf("%s-%s-%s", ns, n.subj, n.act);
    load(sprintf("%s/%s-debug.mat", expDir, name));
    
    sIdx = max(allIdx.w__v(1), allIdx.w__x(1));
    eIdx = min(allIdx.w__v(end), allIdx.w__x(end));
    nSamples = eIdx - sIdx + 1;
    idx = sIdx:eIdx;
    
    viconIdx0 = find(allIdx.w__v==sIdx,1):find(allIdx.w__v==eIdx,1);
    xsensIdx0 = find(allIdx.w__x==sIdx,1):find(allIdx.w__x==eIdx,1);
    
    csActBody = W__viconBody.getSubset(viconIdx0);
    estBody = W__xsensBody.getSubset(xsensIdx0);
    
    idx1EndIdx = csActBody.nSamples;
    found = find(any(isnan(csActBody.MIDPEL), 2), 1);
    if ~isempty(found), idx1EndIdx = min(idx1EndIdx, found-1); end
	found = find(any(isnan(csActBody.LTIO), 2), 1);
    if ~isempty(found), idx1EndIdx = min(idx1EndIdx, found-1); end
	found = find(any(isnan(csActBody.RTIO), 2), 1);
    if ~isempty(found), idx1EndIdx = min(idx1EndIdx, found-1); end
                
    if idx1EndIdx < csActBody.nSamples
        csActBody = csActBody.getSubset(1:idx1EndIdx);
        estBody = estBody.getSubset(1:idx1EndIdx);
        idx = idx(1:idx1EndIdx);
        disp('NaN found in estimator result. Only evaluating until the last non NaN value');
    end
                
    csActBodyRel = csActBody.changeRefFrame('MIDPEL');
    estBodyRel = estBody.changeRefFrame('MIDPEL');
    estBody2 = estBodyRel.toWorldFrame(csActBody.MIDPEL, estBody.qRPV);
    csActBody2 = csActBodyRel.toWorldFrame(csActBody.MIDPEL, csActBody.qRPV);
    
    results0a = estBody.diffRMSEandMean(csActBody);
    results0 = estBody2.diffRMSEandMean(csActBody2);
    
    revStepDetect = readtable(sprintf('%s/rawstep-detect/%s-%s-revStepDetect.csv', ...
                    dir, n.subj, n.act));
    nSamples = csActBody.nSamples;
    intervals = struct('MIDPEL', logical(revStepDetect.stepL(idx)), ...
                       'LTIO', logical(revStepDetect.stepL(idx)), ...
                       'RTIO', logical(revStepDetect.stepR(idx)) );
    results0 = estBody.calcTTDandStepParams(csActBody, intervals, results0);
    if stepParamDir
        estBody.dumpStepParams(csActBody, intervals, ...
            sprintf('%s/%s-%s-viconvsxsens-stepParam', stepParamDir, name, ns));
    end
    
    for j=["MIDPEL", "LTIO", "RTIO"]
        for k=["RMSE", "Std", "Mean"]
            results0.(sprintf("%sW%s",j,k)) = results0a.(sprintf("%s%s",j,k));
        end
    end
            
    results0.dPosW = results0a.dPos;
    results0.name = name;
    results0.label = sprintf("%s+viconvsxsens", ns);
    results0.runtime = 0;
    results0.sim3Dist = 0;
    results0.sim3DistSigma = 0;
        
    results(rIdx) = results0;
    rIdx = rIdx + 1;
    fprintf("%s/%s-%s\n", expDir, name, results0.label);
    
%     targetname = sprintf('%s/%s-viconvsxsens', outDir, name);
%     estBody.exportc3d(sprintf('%s.c3d', targetname), struct(), csActBody);
end

results = struct2table(results);

% Append new results
dataPath = sprintf("%s/results.mat", expDir);
if exist(dataPath, 'file')
    newResults = results;
    load(dataPath);
    [C, ia, ib] = intersect(results(:,{'name', 'label'}), newResults(:,{'name', 'label'}));
    results(ia,:) = [];
    results = [results; newResults];
end
save(sprintf("%s/results.mat", expDir), 'results')

function label = getLabel(ns, setup)
    if setup.accData == 'v'
        if setup.accDataNoise == 0 
            aD = 'v';
        else
            aD = strrep(sprintf('v%.1f', setup.accDataNoise), '.', '');
        end
    else
        aD = setup.accData;
    end
    if strcmp(setup.est, 'ekfv3')
        label = sprintf("%s+A%sO%sI%s+S%s+M%02d+C%03d", ns, aD, ...
            setup.oriData, setup.initSrc, ...
            setup.stepDetection, setup.applyMeas, setup.applyCstr);
    else
        label = sprintf("%s+%s+A%sO%sI%s+S%s+P%03d+M%03d+C%03d", ns, setup.est, ...
            aD, setup.oriData, setup.initSrc, ...
            setup.stepDetection, setup.applyPred, ...
            setup.applyMeas, setup.applyCstr);
    end
end
