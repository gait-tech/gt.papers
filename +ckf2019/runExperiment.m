function results = runExperiment(imu, gfrAcc, step, actBody, dataX, ...
                        setups, savedir, stepParamDir)
    % Run experiment on an instance of Vicon and Xsens dataset taken at NeuRA institute
    % 
    % Setup parameters:
    % - label: data instance name (e.g. s1-acting1 or s2-walking1)
    % - est: filter type to be used.
    %      - ckf: pelib.est.ckf_3imus
    % - accData: acceleration data to be used
    %      - w__v: vicon (world frame)
    %      - w__s: sparse (world frame)
    %      - w__x: xsens
    % - oriData: orientation data to be used
    %      - w__v: vicon (world frame)
    %      - w__s: sparse (world frame)
    %      - w__x: xsens
    % - stepDetection: step detection algorithm to be used
    %      - false: turn off
    %      - av03: use reviewed step detect file. Dimension similar to dataS
    % - initSrc: source of sensor to body orientation and position init
    %      - w__v: vicon (world frame) (default)
    %      - w__x: xsens
    % - applyMeas: measurement configuration number
    % - applyCstr: constraint configuration number
    % - sigmaQAcc: Q acceleration sigma (variance)
    % - P: initial P matrix
    %
    % :param imu: imu measurement, structs of mocapdb.XsensBody
    % :param gfrAcc: structs of inertial acceleration for diff. modes.
    % :param step: step detection table table if stepL and stepR was detected
    % :param actBody: reference body (pelib.grBody)
    % :param dataX: output of xsens reconstruction (BVHBody)
    % :param meas3DistDir: simulation distance measurement directory
    % :param setups: list of experiment parameters (struct) to be run. see details above
    % :param savedir: filepath to save .mat output/debug files (optional)
    % :param stepParamDir: directory to dump step parameters

    %% Inputs and Input Check
    validateattributes(actBody, {'pelib.grBody'}, {});
    validateattributes(dataX, {'mocapdb.BVHBody', 'numeric'}, {});
    
    if nargin <= 6, savedir = ''; end
    if nargin <= 7, stepParamDir = ''; end
    
    %% Initialization   
    % Initialize other variables
    name = actBody.name;
    setupDefault = struct('label', 'ckf', 'est', 'ckf', ...
        'accData', 'w__v', 'accDataNoise', 0.0, 'oriData', 'w__v', ...
        'initSrc', 'w__v', 'stepDetection', 'av03', ...
        'stepDetectWindow', 0.25, 'stepDetectThreshold', 1, ...
        'applyPred', 0, 'applyMeas', 0, 'applyCstr', 0, 'accPelvBias', false, ...
        'sigmaQAcc', 0.5, 'P', 100, 'c3d', false);
    
    setupN = length(setups);
    setupDefaultFN = fieldnames(setupDefault);
        
    for i=1:setupN
        for j=1:length(setupDefaultFN)
            if ~isfield(setups{i}, setupDefaultFN{j})
                setups{i}.(setupDefaultFN{j}) = setupDefault.(setupDefaultFN{j});
            end
        end
    end
            
    %% Run Experiment            
    resultsIdx = 1; clear results;
    
    body = struct('PV_d', actBody.calcPelvisLength(), ...
              'LT_d', actBody.calcLFemurLength(), ...
              'RT_d', actBody.calcRFemurLength(), ...
              'LS_d', actBody.calcLShankLength(), ...
              'RS_d', actBody.calcRShankLength() );
    
    % remove foot data from vicon for consistency
    actBody.LTOE = []; actBody.RTOE = [];
    actBody.qLFT = []; actBody.qRFT = [];

    for sI=1:setupN
        t0 = cputime;
        
        % load setup instances to cs with some preprocessing
        cs = setups{sI};
        if ~(cs.accData(end) == 'v')
            cs.accData = strcat(cs.accData, cs.initSrc(end));
        end
        if ~(cs.oriData(end) == 'v')
            cs.oriData = strcat(cs.oriData, cs.initSrc(end));
        end
        
        %% calculate initial state
        % get initial state from vicon grBody
        vel = actBody.calcJointVel({'MIDPEL', 'LTIO', 'RTIO'});
        x0 = [actBody.MIDPEL(1,:) vel.MIDPEL(1,:) ...
              actBody.LTIO(1,:) vel.LTIO(1,:) ...
              actBody.RTIO(1,:) vel.RTIO(1,:) ]'; 
          
        % Backtrack x0 such that x_pred at t=1 is x0 from viconBody
        dt = 1.0/actBody.fs;
        dt2 = dt*dt;
        x0(4:6,:) = x0(4:6,:) - gfrAcc.(cs.accData).Pelvis(1,:)'*dt;
        x0(10:12,:) = x0(10:12,:) - gfrAcc.(cs.accData).L_LowLeg(1,:)'*dt;
        x0(16:18,:) = x0(16:18,:) - gfrAcc.(cs.accData).R_LowLeg(1,:)'*dt;
        x0(1:3,:) = x0(1:3,:) - x0(4:6,:)*dt - 0.5*gfrAcc.(cs.accData).Pelvis(1,:)'*dt2;
        x0(7:9,:) = x0(7:9,:) - x0(10:12,:)*dt - 0.5*gfrAcc.(cs.accData).L_LowLeg(1,:)'*dt2;
        x0(13:15,:) = x0(13:15,:) - x0(16:18,:)*dt - 0.5*gfrAcc.(cs.accData).R_LowLeg(1,:)'*dt2;
%         x0(4:6, :) = (x0(4:6, :)*dt - 0.5*gfrAcc.(cs.accData).Pelvis(1,:)'*dt^2)/dt;
%         x0(10:12, :) = (x0(10:12, :)*dt - 0.5*gfrAcc.(cs.accData).L_LowLeg(1,:)'*dt^2)/dt;
%         x0(16:18, :) = (x0(16:18, :)*dt - 0.5*gfrAcc.(cs.accData).R_LowLeg(1,:)'*dt^2)/dt;

        % step detection
        isStep = struct;
        if strcmp(cs.stepDetection, 'av03')
            isStep.MIDPEL = false(actBody.nSamples, 1);
            isStep.LTIO = logical(step.stepL);
            isStep.RTIO = logical(step.stepR);
        else
            isStep.MIDPEL = false(actBody.nSamples, 1);
            isStep.LTIO = false(actBody.nSamples, 1);
            isStep.RTIO = false(actBody.nSamples, 1);
        end
        
        % if the init position is negative knee angle, allow it
        alphaLKmin = actBody.calcJointAnglesLKnee(1);
        alphaLKmin = min(alphaLKmin(2), 0);
        alphaRKmin = actBody.calcJointAnglesRKnee(1);
        alphaRKmin = min(alphaRKmin(2), 0);
        
        %% big switch statement for different algorithms
        if strcmp(cs.est, 'ckf')
            v3Options = struct('fs', actBody.fs, 'applyMeas', cs.applyMeas, ...
                'applyCstr', cs.applyCstr, 'sigma2QAccMP', cs.sigma2QAcc, ...
                'sigma2QAccLA', cs.sigma2QAcc, 'sigma2QAccRA', cs.sigma2QAcc, ...
                'alphaLKmin', alphaLKmin, 'alphaRKmin', alphaRKmin);

            [~, x_pos_v2, t_dat_v2 ] = papers.ckf2019.ckf_3imus(x0, cs.P, ...
                gfrAcc.(cs.accData).Pelvis, isStep.MIDPEL, imu.(cs.oriData).Pelvis.ori, ...
                gfrAcc.(cs.accData).L_LowLeg, isStep.LTIO, imu.(cs.oriData).L_LowLeg.ori, ...
                gfrAcc.(cs.accData).R_LowLeg, isStep.RTIO, imu.(cs.oriData).R_LowLeg.ori, ...
                body.PV_d, body.LT_d, body.RT_d, body.LS_d, body.RS_d, ...
                v3Options);
            
            estBody = pelib.grBody('name', 'est', 'posUnit', 'm', 'oriUnit', 'deg', ...
               'lnSymbol', '--', 'ptSymbol', 'o', 'frame', 'world', ...
               'xyzColor', {'r', 'g', 'b'}, 'fs', actBody.fs, ...
               'MIDPEL', x_pos_v2(:, 1:3), 'LFEP', t_dat_v2.LFEP, ...
               'LFEO', t_dat_v2.LFEO, 'LTIO', x_pos_v2(:, 7:9), ...
               'RFEP', t_dat_v2.RFEP, ...
               'RFEO', t_dat_v2.RFEO, 'RTIO', x_pos_v2(:, 13:15), ...
               'qRPV', imu.(cs.oriData).Pelvis.ori, ...
               'qLTH', t_dat_v2.qLTH, 'qRTH', t_dat_v2.qRTH, ...
               'qLSK', imu.(cs.oriData).L_LowLeg.ori, ...
               'qRSK', imu.(cs.oriData).R_LowLeg.ori);

            estState = x_pos_v2;
            estState2 = t_dat_v2;
        elseif strcmp(cs.est, 'vsxsens')
            estBody = dataX.togrBody(1:actBody.nSamples, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', actBody.fs, ...
                         'xyzColor', {'m', 'y', 'c'}});
            estBody.LTOE = []; estBody.RTOE = [];
            estBody.qLFT = []; estBody.qRFT = [];
            estState = [];
            estState2 = [];
        end

        runtime = cputime-t0;

        % save result
        if ~strcmp(savedir, '')
            save(sprintf("%s/%s-%s.mat", savedir, name, cs.label), ...
                 'estBody', 'estState', 'estState2', 'runtime', 'cs')
        end
        
        if cs.c3d
            estBody.exportc3d(sprintf("%s/%s-%s.c3d", savedir, name, cs.label), ...
                struct(), actBody);
        end
        
        %% evaluation
        targetSeg = {'qLTH', 'qRTH'};
        
        % convert grBody to midpel frame
        estBodyRel = estBody.changeRefFrame('MIDPEL');
        % revert back to global frame but with MIDPEL position equal to actBody
        estBody2 = estBodyRel.toWorldFrame(actBody.MIDPEL, estBody.qRPV);

        % calculate kinematic gait parameters against actBody (vicon)
        r_buf = estBody2.diffRMSEandMean(actBody, true, targetSeg);
        
        % calculate some spatiotemporal parameters
        r_bufW = estBody.diffRMSEandMean(actBody, true, targetSeg);
        r_buf = estBody.calcTTDandStepParams(actBody, isStep, r_buf);
        estBody.dumpStepParams(actBody, isStep, ...
            sprintf('%s/%s-%s-stepParam', stepParamDir, name, cs.label));
        for i=["MIDPEL", "LTIO", "RTIO"]
            for j=["RMSE", "Std", "Mean"]
                r_buf.(sprintf("%sW%s",i,j)) = r_bufW.(sprintf("%s%s",i,j));
            end
        end
        r_buf.dPosW = r_bufW.dPos;   
        
        % labeling
        r_buf.name = name;
        r_buf.label = cs.label;
        r_buf.runtime = runtime;
        
        results(resultsIdx) = r_buf;
        fprintf("Index %3d/%3d: Running time: %.4f\n", resultsIdx, setupN, cputime-t0);
        resultsIdx = resultsIdx + 1;
    end
end