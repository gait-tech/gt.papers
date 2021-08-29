function results = runRawNeuRASparse01Experiment(dataS, dataV, ...
                        calibV2W, calibYawFix, calibW2V, ...
                        dataX, revStepDetect, meas3Dist, ...
                        name, setups, savedir, startFrame, endFrame, bias)
    % Run experiment on an instance of Vicon and Xsens dataset taken at NeuRA institute
    % 
    % Setup parameters:
    % - label: data instance name (e.g. s1-acting1 or s2-walking1)
    % - est: filter type to be used.
    %      - ekfv3: pelib.est.kf_3_kmus_v3
    %      - pfv1: pelib.est.pf_3_kmus_v1
    %      - mpfv1: pelib.est.mpf_3_kmus_v1
    %      - lieekfv1: pelib.est.lieekf_3_kmus_v1
    %      - lieekfv2: pelib.est.lieekf_3_kmus_v2
    %      - lgckf7sv1: pelib.est.lgcekf7seg_3imus_v1
    %      - ckf: pelib.est.ckf_3imus
    %      - ckfdist: papers.ckfdist2020.ckfdist_3imus
    % - accData: acceleration data to be used
    %      - w__v: vicon (world frame)
    %      - w__s: sparse (world frame)
    %      - v__v: vicon (vicon frame)
    %      - v__s: sparse (vicon frame)
    %      - w__x: xsens
    % - oriData: orientation data to be used
    %      - w__v: vicon (world frame)
    %      - w__s: sparse (world frame)
    %      - v__v: vicon (vicon frame)
    %      - v__s: sparse (vicon frame)
    %      - w__x: xsens
    % - stepDetection: step detection algorithm to be used
    %      - false: turn off
    %      - av01: fixed acceleration variance on tibia accData (var = 1)
    %      - av02: fixed acceleration variance on vicon tibia accData (var = 1)
    %      - av03: use reviewed step detect file. Dimension similar to dataS
    % - initSrc: source of sensor to body orientation and position init
    %      - w__v: vicon (world frame) (default)
    %      - v__v: vicon (vicon frame)
    %      - w__x: xsens
    % - applyMeas: measurement configuration number
    % - applyCstr: constraint configuration number
    % - sigma2QAcc: Q acceleration sigma (variance)
    % - P: initial P matrix
    % - sim3Dist: 0 off, else simulation id
    % - sim3DistSigma: sigma of simulation
    %
    % :param dataS: loaded mocapdb.XsensBody
    % :param dataV: loaded mocapdb.ViconBody
    % :param calibV2W: quaternion (1 x 4) transforming vicon frame to world frame
    % :param calibYawFix: mocapdb.XsensBody fix ankle sensor yaw orientation offset
    % :param calibW2V: mocapdb.XsensBody transforming each sensor's world frame 
    %                  to vicon frame. Includes yaw realignment calibration.
    % :param dataX: loaded mocapdb.BVHBody 
    % :param revStepDetect: manually reviewed table if stepL and stepR was detected
    % :param meas3Dist: 3 point distance measurement. if str or chr, load appropriate simulation file
    % :param name: name of the experiment
    % :param setups: list of experiment parameters (struct) to be run. see details above
    % :param savedir: filepath to save .mat output/debug files (optional)
    % :param startFrame: frame number at which the algorithm will start
    % :param endFrame: frame number at which the algorithm will end
    % :param bias: pelvis accelerometer bias in sensor frame
    % 
    % .. Author: - Luke Wicent Sy (GSBME, Modified 2019 Oct 31)
    
    %% Inputs and Input Check
    validateattributes(dataS, {'mocapdb.XsensBody'}, {});
    validateattributes(dataV, {'mocapdb.ViconBody', 'numeric'}, {});
    validateattributes(calibV2W, {'numeric'}, {});
    validateattributes(calibW2V, {'mocapdb.XsensBody', 'numeric'}, {});
    validateattributes(dataX, {'mocapdb.BVHBody', 'numeric'}, {});
    
    if nargin <= 7, meas3Dist = []; end
    if nargin <= 10, savedir = ''; end
    if nargin <= 11 || startFrame < 0, startFrame = 100; end
    if nargin <= 12 || endFrame < 0, endFrame = inf; end
    if nargin <= 13
        bias = struct('w__v', zeros(1, 3), 'v__v', zeros(1, 3), ...
                      'w__x', zeros(1, 3));
    end
    
    %% Initialization   
    % Initialize other variables
    fs = dataS.fs;
    setupDefault = struct('label', 'ekfv3', 'est', 'ekfv3', ...
        'accData', 'v', 'accDataNoise', 0.0, 'oriData', 'v', ...
        'initSrc', 'v', 'stepDetection', 'av01', ...
        'stepDetectWindow', 0.25, ...
        'stepDetectThreshold', struct('MP', 0, 'LA', 1, 'RA', 1, 'LF', 1, 'RF', 1), ...
        'applyPred', 0, 'applyMeas', 0, 'applyCstr', 0, 'accPelvBias', false, ...
        'sigmaUwbLeg', 1e0, 'c3d', false, ...
        'sigmaQAcc', 0.5, 'P', 100, ...
        'sim3Dist', 0, 'sim3DistSigma', 0.0);
    isMeas3DistDir = ischar(meas3Dist) | isstring(meas3Dist);
    if isMeas3DistDir, meas3DistDir = meas3Dist; end
    
    setupN = length(setups);
    setupDefaultFN = fieldnames(setupDefault);
        
    for i=1:setupN
        for j=1:length(setupDefaultFN)
            if ~isfield(setups{i}, setupDefaultFN{j})
                setups{i}.(setupDefaultFN{j}) = setupDefault.(setupDefaultFN{j});
            end
        end
    end

%     sim3DistSigma = 0.1; %m standard deviation
    
    allIdx = {};
    qOri = {};
    gfrAcc = {};
    
    wbodyOri = {};
    bodyAcc = {};
    x0 = {};
    
    %% Helper variables
    calibYawFixList = ["L_LowLeg", "R_LowLeg"];
    bodyKeyList = struct();
    bodyKeyList.all = [ ...
       struct('ks', 'PELV', 'kp', 'MP', 'x', 'Pelvis', 'xs', 'qHips', 'vs', 'qRPV', 'vp', 'MIDPEL'), ...
       struct('ks', 'LTIB', 'kp', 'LA', 'x', 'L_LowLeg', 'xs', 'qLeftLeg', 'vs', 'qLSK', 'vp', 'LTIO'), ...
       struct('ks', 'RTIB', 'kp', 'RA', 'x', 'R_LowLeg', 'xs', 'qRightLeg', 'vs', 'qRSK', 'vp', 'RTIO')
    ];
    for i = ["ks", "kp", "x", "vs", "vp"]
        bodyKeyList.(i) = arrayfun(@(x) x.(i), bodyKeyList.all, 'UniformOutput', false);
    end
    
    %% Preprocessing in world frame
    if ~isempty(dataV) && ~isempty(calibV2W)
        nSamples = min(dataV.nSamples, dataS.nSamples);
        W__dataV = dataV.getSubset(1:nSamples).toWorldFrame(calibV2W);
        W__dataV.changePosUnit('m', true);
        W__dataS = dataS.getSubset(1:nSamples);
        W__dataS.Pelvis.acc = W__dataS.Pelvis.acc - bias.w__v;
        % apply yaw offset to orientation
        for i=calibYawFixList
            if ~isempty(calibYawFix.(i))
                W__dataS.(i).ori = quatmultiply(calibYawFix.(i).ori, W__dataS.(i).ori);
            end
        end
        
        sIdx = max(W__dataV.getStartIndex()+1, startFrame);
        eIdx = min(length(W__dataV.PELV(:,1)) - 1, endFrame);
        idx = sIdx:eIdx; idx0 = 1:(eIdx-sIdx+1);
        allIdx.w__v = idx;
        
        viconCalibSB = W__dataS.calcCalibSB(W__dataV.togrBody(sIdx+1:sIdx+1, {}), sIdx(1));       
        %% orientation and angular velocity
        % Angular velocity of body in body frame as obtained from vicon input
        W__viconBody = W__dataV.togrBody(1:nSamples, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                         'xyzColor', {'m', 'y', 'c'}});
        % angvel = W__viconBody.calcSegAngVel({'qRPV', 'qLSK', 'qRSK'}, 'B');
        angvel = W__viconBody.calcSegAngVel(bodyKeyList.vs, 'B');
        for i=bodyKeyList.all
            % Orientation of body in world frame as obtained from sparse sensor
            qBufferEst0 = quatmultiply(W__dataS.(i.x).ori, quatconj(viconCalibSB.(i.x).ori));
            qOri.w__sv.(i.ks) = qBufferEst0(sIdx:eIdx, :);
            
            % https://math.stackexchange.com/questions/2282938/converting-from-quaternion-to-angular-velocity-then-back-to-quaternion
            % Angular velocity of body in body frame as obtained from sparse sensor
            wbodyOri.w__sv.(i.ks) = quatrotate(quatconj(viconCalibSB.(i.x).ori), W__dataS.(i.x).gyr);
            wbodyOri.w__sv.(i.ks) = wbodyOri.w__sv.(i.ks)(sIdx:eIdx, :);
            
            % Orientation of body in world frame as obtained from vicon input
            qOri.w__v.(i.ks) = W__dataV.(i.vs)(sIdx+1:eIdx+1, :);
            
            % Angular velocity of body in body frame as obtained from vicon input
            wbodyOri.w__v.(i.ks) = angvel.(i.vs)(sIdx+1:eIdx+1,:);
        end
                                      
        %% position, velocity, acceleration
        vel = W__viconBody.calcJointVel(bodyKeyList.vp); % {'MIDPEL', 'LTIO', 'RTIO'}
        acc = W__viconBody.calcJointAcc(bodyKeyList.vp); % {'MIDPEL', 'LTIO', 'RTIO'}

        x0.w__v = [W__viconBody.MIDPEL(sIdx,:) vel.MIDPEL(sIdx,:) zeros(1,4) ...
                   W__viconBody.LTIO(sIdx,:) vel.LTIO(sIdx,:) zeros(1,4) ...
                   W__viconBody.RTIO(sIdx,:) vel.RTIO(sIdx,:) zeros(1,4)]';     
        
        gfrAcc.w__sv = {};
        % gfrAcc from filtered sparse
        fc = 10;
        [lpf_b, lpf_a] = butter(6, fc/(fs/2));

        vsigma = unique([cellfun(@(x) x.accDataNoise, setups), 0]);
        randnN = size(acc.MIDPEL, 1);
        for i = 1:length(vsigma)
            vLabel = getVLabel('w__v', vsigma(i));
            gfrAcc.(vLabel) = {};
            for j = bodyKeyList.all
                gfrAcc.(vLabel).(j.kp) = acc.(j.vp) + randn(randnN,3).*vsigma(i);
                gfrAcc.(vLabel).(j.kp) = gfrAcc.(vLabel).(j.kp)(sIdx:eIdx,:);
            end
        end
        
        for i = bodyKeyList.all
            % gfrAcc from sparse
            gfrAcc.w__sv.(i.kp) = quatrotate(quatconj(W__dataS.(i.x).ori), ...
                                    W__dataS.(i.x).acc) - [0 0 9.81];
            gfrAcc.w__sv.(i.kp) = gfrAcc.w__sv.(i.kp)(sIdx:eIdx,:);

            % gfrAcc from filtered sparse
            gfrAcc.w__sfv.(i.kp) = filter(lpf_b, lpf_a, gfrAcc.w__sv.(i.kp));

            %% body acceleration
            bodyAcc.w__sv.(i.kp) = quatrotate(quatconj(viconCalibSB.(i.x).ori), W__dataS.(i.x).acc);
            bodyAcc.w__sv.(i.kp) = bodyAcc.w__sv.(i.kp)(sIdx:eIdx, :);
            
            bodyAcc.w__v.(i.kp) = quatrotate(qOri.w__v.(i.ks), gfrAcc.w__v.(i.kp) + [0 0 9.81]);
        end
            
        %% measurement of 3 distance (pelvis to ankles, ankle to ankle)
        %  Simulate distance measurement by generating pairwise combinations
        for sI=1:setupN
            cs = setups{sI};
            if (cs.sim3Dist > 0) && isMeas3DistDir && strcmp(cs.initSrc, 'w__v')
                meas3DistFname = sprintf("%s/%s-I%s-%.3f-%03d.csv", ...
                                         meas3DistDir, name, cs.initSrc,  ...
                                         cs.sim3DistSigma, cs.sim3Dist);
                if ~exist(meas3DistFname, 'file')
                    m3dBuf = table;
                    m3dBuf.PV_LA = vecnorm(W__viconBody.MIDPEL-W__viconBody.LTIO, 2, 2) ...
                        + normrnd(0, cs.sim3DistSigma, [nSamples, 1]);
                    m3dBuf.PV_RA = vecnorm(W__viconBody.MIDPEL-W__viconBody.RTIO, 2, 2) ...
                        + normrnd(0, cs.sim3DistSigma, [nSamples, 1]);
                    m3dBuf.LA_RA = vecnorm(W__viconBody.RTIO-W__viconBody.LTIO, 2, 2) ...
                        + normrnd(0, cs.sim3DistSigma, [nSamples, 1]);
                    m3dBuf.LLeg = vecnorm(W__viconBody.LFEP-W__viconBody.LTIO, 2, 2) ...
                        + normrnd(0, cs.sim3DistSigma, [nSamples, 1]);
                    m3dBuf.RLeg = vecnorm(W__viconBody.RFEP-W__viconBody.RTIO, 2, 2) ...
                        + normrnd(0, cs.sim3DistSigma, [nSamples, 1]);
                    
                    writetable(m3dBuf, meas3DistFname);
                end
            end
        end
    
        % debug purposes
        W__viconBody = W__dataV.togrBody(idx+1, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                         'xyzColor', {'m', 'y', 'c'}});
        PV__viconBody = W__viconBody.changeRefFrame('MIDPEL');
    end
    
    if ~isempty(dataX)
        nSamples = min(dataX.nSamples, dataS.nSamples);
        qXsensV2W = rotm2quat([0 0 1; 1 0 0; 0 1 0]);
        
        dataX = dataX.toWorldFrame(qXsensV2W);
        W__dataX = dataX.getSubset(1:nSamples);
        W__dataX.changePosUnit('m', true);
        W__dataS = dataS.getSubset(1:nSamples); % .toViconFrame(calibW2V);
        % order is not important as calibW2V fixes only the ankle yaw offset
        W__dataS.Pelvis.acc = W__dataS.Pelvis.acc - bias.w__x; 
        % apply yaw offset to orientation
        for i=calibYawFixList
            if ~isempty(calibYawFix.(i))
                W__dataS.(i).ori = quatmultiply(calibYawFix.(i).ori, W__dataS.(i).ori);
            end
        end
       
        sIdx = startFrame;
        eIdx = min(length(W__dataX.Hips(:,1)) - 1, endFrame);
        idx = sIdx:eIdx; idx0 = 1:(eIdx-sIdx+1);
        allIdx.w__x = idx;
        xsensCalibSB = W__dataS.calcCalibSB(W__dataX.togrBody(sIdx+1:sIdx+1, {}), sIdx(1)); 
        
        % debug purposes
        W__xsensBody = W__dataX.togrBody(idx+1, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', fs, ...
                         'xyzColor', {'m', 'y', 'c'}});
        PV__xsensBody = W__xsensBody.changeRefFrame('MIDPEL');
    end
    
    
    %% Save processing
    if ~strcmp(savedir, '')
        if ~isempty(dataX)
            save(sprintf("%s/%s-debug.mat", savedir, name), ...
                 'W__viconBody', 'W__xsensBody', ...
                 'gfrAcc', 'qOri', 'bodyAcc', 'wbodyOri', ...
                 'x0', 'allIdx')
        else
            save(sprintf("%s/%s-debug.mat", savedir, name), ...
                 'W__viconBody', ...
                 'gfrAcc', 'qOri', 'bodyAcc', 'wbodyOri', ...
                 'x0', 'allIdx')
        end
    end
            
    %% Run Experiment            
    resultsIdx = 1; clear results;
    
    for sI=1:setupN
        t0 = cputime;
        
        cs = setups{sI};
        
        idx = allIdx.(cs.initSrc);
        sIdx = idx(1);
        eIdx = idx(end);
        idx0 = 1:(eIdx-sIdx+1);
        
        if cs.accData(end) == 'v'
            csGfrAcc = gfrAcc.(getVLabel(cs.accData, cs.accDataNoise));
        elseif ( strcmp(cs.accData, 'w__s') || strcmp(cs.accData, 'v__s') || ...
           strcmp(cs.accData, 'w__sf') || strcmp(cs.accData, 'v__sf') )
            csGfrAcc = gfrAcc.(strcat(cs.accData, cs.initSrc(end)));
        else
            csGfrAcc = gfrAcc.(cs.accData);
        end
        
        if strcmp(cs.oriData, 'w__s') || strcmp(cs.oriData, 'v__s')
            csQOri = qOri.(strcat(cs.oriData, cs.initSrc(end)));
        else
            csQOri = qOri.(cs.oriData);
        end
        
        if ( strcmp(cs.accData, 'w__s') || strcmp(cs.accData, 'v__s') || ...
           strcmp(cs.accData, 'w__sf') || strcmp(cs.accData, 'v__sf') )
            % init velocity adjustment so pos at t=1 is equal to
            % vicon/xsens body
            csx0 = x0.(cs.initSrc);
            dt = 1.0/fs;
            csx0(4:6, :) = (csx0(4:6, :)*dt - 0.5*csGfrAcc.MP(1,:)'*dt^2)/dt;
            csx0(14:16, :) = (csx0(14:16, :)*dt - 0.5*csGfrAcc.LA(1,:)'*dt^2)/dt;
            csx0(24:26, :) = (csx0(24:26, :)*dt - 0.5*csGfrAcc.RA(1,:)'*dt^2)/dt;
        else
            csx0 = x0.(cs.initSrc);
        end
        
        body = struct();
        if strcmp(cs.initSrc, 'w__v')
            csActBody = W__viconBody;
            csActBodyRel = PV__viconBody;
        elseif strcmp(cs.initSrc, 'v__v')
            csActBody = V__viconBody;
            csActBodyRel = PV__viconBody;
        else
            csActBody = W__xsensBody;
            csActBodyRel = PV__xsensBody;
        end
        body = struct('PV_d', csActBody.calcPelvisLength(), ...
                      'LT_d', csActBody.calcLFemurLength(), ...
                      'RT_d', csActBody.calcRFemurLength(), ...
                      'LS_d', csActBody.calcLShankLength(), ...
                      'RS_d', csActBody.calcRShankLength() );
                          
        % step detection
        bIsStat = struct();
        if strcmp(cs.stepDetection, 'av03')
            csNSamples = size(csGfrAcc.MP, 1);
            bIsStat.MP = false(csNSamples, 1);
            bIsStat.LA = revStepDetect.stepL(idx);
            bIsStat.RA = revStepDetect.stepR(idx);
            bIsStat.LF = revStepDetect.stepL(idx);
            bIsStat.RF = revStepDetect.stepR(idx);
        else
            csNSamples = size(csGfrAcc.MP, 1);
            for i = bodyKeyList.all
                bIsStat.(i.kp) = false(csNSamples, 1);
            end
        end
        
        % measurement 3 distance
        if (cs.sim3Dist > 0) && isMeas3DistDir
            meas3DistFname = sprintf("%s/%s-I%s-%.3f-%03d.csv", ...
                                     meas3DistDir, name, cs.initSrc,  ...
                                     cs.sim3DistSigma, cs.sim3Dist);
            meas3Dist = readtable(meas3DistFname);
            meas3Dist = meas3Dist(idx,:);
        else
            meas3Dist = struct('PV_LA', [], 'PV_RA', [], 'LA_RA', []);
        end
                    
        % if the init position is negative knee angle, allow it
        alphaLKmin = csActBodyRel.calcJointAnglesLKnee(1);
        alphaLKmin = min(alphaLKmin(2), 0);
        alphaRKmin = csActBodyRel.calcJointAnglesRKnee(1);
        alphaRKmin = min(alphaRKmin(2), 0);
        
        if strcmp(cs.est, 'ckfdist')
            v3Options = struct('fs', fs, 'applyMeas', cs.applyMeas, ...
                'applyCstr', cs.applyCstr, 'sigma2QAccMP', cs.sigma2QAcc, ...
                'sigma2QAccLA', cs.sigma2QAcc, 'sigma2QAccRA', cs.sigma2QAcc, ...
                'alphaLKmin', alphaLKmin, 'alphaRKmin', alphaRKmin);
%                 disp(sprintf('%.2f %.2f', rad2deg(alphaLKmin), rad2deg(alphaRKmin)));
            csx0 = csx0([1:6 11:16 21:26], :);
            [ x_pri_v2, x_pos_v2, t_dat_v2 ] = papers.ckfdist2020.ckfdist_3imus( ...
                csx0, cs.P, csGfrAcc.MP, bIsStat.MP, csQOri.PELV, ...
                csGfrAcc.LA, bIsStat.LA, csQOri.LTIB, ...
                csGfrAcc.RA, bIsStat.RA, csQOri.RTIB, ...
                body, meas3Dist, v3Options);

            idx1EndIdx = csActBody.nSamples;
            found = find(any(isnan(x_pos_v2), 2), 1);
            if ~isempty(found), idx1EndIdx = found; end
            found = find(any(isnan(csActBody.MIDPEL), 2), 1);
            if ~isempty(found), idx1EndIdx = min(idx1EndIdx,found); end
            found = find(any(isnan(csActBody.LTIO), 2), 1);
            if ~isempty(found), idx1EndIdx = min(idx1EndIdx,found); end
            found = find(any(isnan(csActBody.RTIO), 2), 1);
            if ~isempty(found), idx1EndIdx = min(idx1EndIdx,found); end

            if idx1EndIdx < csActBody.nSamples
                idx1 = idx0(1):idx0(idx1EndIdx-1); 
                csActBody = csActBody.getSubset(idx1);
                csActBodyRel = csActBodyRel.getSubset(idx1);
                disp('NaN found in estimator result. Only evaluating until the last non NaN value');
            else
                idx1 = idx0;
            end

            estBody = papers.ckfdist2020.exportGrBody(x_pos_v2, t_dat_v2);
            if estBody.nSamples > csActBody.nSamples
                estBody = estBody.getSubset(1:csActBody.nSamples);
            end
            if cs.c3d
                [acq, acqDebug] = papers.ckfdist2020.exportc3d( ...
                    x_pos_v2, t_dat_v2, ...
                    sprintf("%s/%s-%s", savedir, name, cs.label), ...
                    csActBody, bIsStat, dataS, idx(idx1));
            end

            estState = x_pos_v2;
            estState2 = t_dat_v2;
        end

        runtime = cputime-t0;
        estBody.nSamples = length(idx1);
        estBodyRel = estBody.changeRefFrame('MIDPEL');
        if ~strcmp(savedir, '')
            save(sprintf("%s/%s-%s.mat", savedir, name, cs.label), ...
                 'estBody', 'estState', 'estState2', 'runtime', 'cs')
        end

        estBody2 = estBodyRel.toWorldFrame(csActBody.MIDPEL, estBody.qRPV);
        csActBody2 = csActBodyRel.toWorldFrame(csActBody.MIDPEL, csActBody.qRPV);
%         results(resultsIdx) = estBody.diffRMSE(csActBody);
        results0a = estBody.diffRMSEandMean(csActBody);
        results0 = estBody2.diffRMSEandMean(csActBody2);

        intervals = struct('MIDPEL', logical(bIsStat.LA(idx1)), ...% false(n, 1)
                           'LTIO', logical(bIsStat.LA(idx1)), ...
                           'RTIO', logical(bIsStat.RA(idx1)) );
        results0 = estBody.calcTTD(csActBody, intervals, results0);
        
        results0.dPosW = results0a.dPos;   
        results0.name = name;
        results0.runtime = runtime;
        results0.label = cs.label;
        results0.sim3Dist = cs.sim3Dist;
        results0.sim3DistSigma = cs.sim3DistSigma;
        results(resultsIdx) = results0;
        fprintf("Index %3d/%3d: Running time: %.4f\n", resultsIdx, setupN, cputime-t0);
        resultsIdx = resultsIdx + 1;
    end
end

function label = getVLabel(data_source, sigma)
    if data_source(end) == 'v'
        if sigma == 0
            label = data_source;
        else
            label = strrep(sprintf('%s%.1f', data_source, sigma), '.', '');
        end
    else
        label = data_source;
    end
end
