% ======================================================================
%> Run experiment for all NeuRA data
%> writetable(results, 'C:/Users/z5151460/OneDrive - UNSW/Thesis - Sparse Mocap/Goal01b-CKFDist/Analysis-ckfdistv1-20200109/results20200109.csv')
% ======================================================================
dir = 'data/neura-sparse01';
expDir = sprintf('%s/output', dir);
mkdir(expDir);
mkdir(sprintf('%s/mat', dir));

addpath('mod-lib');

DEGRANGE = (0:0.1:359) - 180;
dataList = readtable(sprintf('+papers/+ckfdist2020/data-list.csv', dir));
dataN = size(dataList, 1);
DATARANGE = 1:dataN;
% DATARANGE = 34;

meas3DistLoc = '+papers/+ckfdist2020/+dist_sim';
mkdir(meas3DistLoc);

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

    for sI = [ struct('algo', 'ckfdist', 'pI', 1, 'mI',  70, 'cI', 355) , ...
               struct('algo', 'ckfdist', 'pI', 1, 'mI',  76, 'cI', 355) , ...
               struct('algo', 'ckfdist', 'pI', 1, 'mI',  170, 'cI', 355) , ...
               struct('algo', 'ckfdist', 'pI', 1, 'mI',  270, 'cI', 355) , ...
               struct('algo', 'ckfdist', 'pI', 1, 'mI',  130, 'cI', 355) , ...
               struct('algo', 'ckfdist', 'pI', 1, 'mI',  230, 'cI', 355) , ...
               ] 
        for sdI = {'av03'}
            for initI = {'w__v'}
                s3dI = bitand(mod(idivide(int32(sI.mI), 100, 'floor'), 10), 3) > 0;
                setups{end+1} = struct('est', sI.algo, ...
                           'accData', 'w__s', 'oriData', 'w__s', 'accDataNoise', 0, ...
                           'initSrc', initI, 'stepDetection', sdI, ...
                           'applyPred', sI.pI, 'applyMeas', sI.mI, ...
                           'applyCstr', sI.cI, 'P', 0.5, ...
                           'sigma2QAcc', 1e2, 'c3d', false, ...
                           'sim3Dist', s3dI, 'sim3DistSigma', 0);
                setups{end+1} = struct('est', sI.algo, ...
                           'accData', initI, 'oriData', initI, 'accDataNoise', 0, ...
                           'initSrc', initI, 'stepDetection', sdI, ...
                           'applyPred', sI.pI, 'applyMeas', sI.mI, ...
                           'applyCstr', sI.cI, 'P', 0.5, ...
                           'sigma2QAcc', 1e2, 'c3d', false, ...
                           'sim3Dist', s3dI, 'sim3DistSigma', 0);
            end
        end
    end
    for sI = [ ... % struct('algo', 'ckfdist', 'pI', 1, 'mI', 130, 'cI', 355) , ...
               ... % struct('algo', 'ckfdist', 'pI', 1, 'mI', 230, 'cI', 355) , ... 
               struct('algo', 'ckfdist', 'pI', 1, 'mI', 170, 'cI', 355) , ...
               struct('algo', 'ckfdist', 'pI', 1, 'mI', 270, 'cI', 355) , ... 
               ] 
        for sdI = {'av03'}
            for initI = {'w__v'}
                for s3dSigma = [0.01:0.01:0.1 0.15:0.05:0.2]
                    for s3dI = 1 % 1:5
                        setups{end+1} = struct('est', sI.algo, ...
                                   'accData', 'w__s', 'oriData', 'w__s', 'accDataNoise', 0, ...
                                   'initSrc', initI, 'stepDetection', sdI, ...
                                   'applyPred', sI.pI, 'applyMeas', sI.mI, ...
                                   'applyCstr', sI.cI, 'P', 0.5, ...
                                   'sigma2QAcc', 1e2, 'c3d', true, ...
                                   'sim3Dist', s3dI, 'sim3DistSigma', s3dSigma);
%                         setups{end+1} = struct('est', sI.algo, ...
%                                    'accData', initI, 'oriData', initI, 'accDataNoise', 0, ...
%                                    'initSrc', initI, 'stepDetection', sdI, ...
%                                    'applyPred', sI.pI, 'applyMeas', sI.mI, ...
%                                    'applyCstr', sI.cI, 'P', 0.5, ...
%                                    'sigma2QAcc', 1e2, 'c3d', false, ...
%                                    'sim3Dist', s3dI, 'sim3DistSigma', s3dSigma);
                    end
                end
            end
        end
    end

    for i = 1:length(setups)
        setups{i}.label = getLabel(ns, setups{i});
    end

    dataN = size(dataList, 1);

    for i = DATARANGE
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
                sprintf('%s/calib/%s-Calib-SensorYawFixWorldFrame.txt', ...
                        dir, n.subj), ...
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
                W__dataV = data.dataV.toWorldFrame(data.calibV2W);
                W__dataV.changePosUnit('m', true);
                data.calibYawFix = data.dataS.calcCalibAnkleSensorW2PelvisWFromVicon(W__dataV);
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
        r = papers.ckfdist2020.runRawNeuRASparse01Experiment(data.dataS, ...
                data.dataV, data.calibV2W, data.calibYawFix, data.calibW2V, ...
                data.dataX, data.revStepDetect, meas3DistLoc, ...
                data.name, setups, expDir, n.startFrame, n.endFrame, data.bias);
        results = [results; struct2table(r)];
    end
end

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