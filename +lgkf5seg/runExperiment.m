function results = runExperiment(imu, gfrAcc, step, actBody, dataX, ...
                        setups, savedir, stepParamDir, meas3DistDir)
    % Run experiment on an instance of Vicon and Xsens dataset taken at NeuRA institute
    % 
    % Setup parameters:
    % - label: data instance name (e.g. s1-acting1 or s2-walking1)
    % - est: filter type to be used.
    %      - ckf: papers.ckf2019.ckf_3imus
    %      - ckfdist: papers.ckfdist2020.ckfdist_3imus
    %      - lieekf: papers.lgkf5seg.lieekf_3_kmus
    % - accData: acceleration data to be used
    %      - w__v: vicon (world frame)
    %      - w__s: sparse (world frame)
    % - oriData: orientation data to be used
    %      - w__v: vicon (world frame)
    %      - w__s: sparse (world frame)
    % - stepDetection: step detection algorithm to be used
    %      - false: turn off
    %      - av03: use reviewed step detect file. Dimension similar to dataS
    % - initSrc: source of sensor to body orientation and position init
    %      - w__v: vicon (world frame) (default)
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
    % :param setups: list of experiment parameters (struct) to be run. see details above
    % :param savedir: filepath to save .mat output/debug files (optional)
    % :param stepParamDir: directory to dump step parameters
    % :param meas3DistDir: simulation distance measurement directory


    %% Inputs and Input Check
    validateattributes(actBody, {'pelib.grBody'}, {});
    validateattributes(dataX, {'mocapdb.BVHBody', 'numeric'}, {});
    
    if nargin <= 6, savedir = ''; end
    if nargin <= 7, stepParamDir = ''; end
    if nargin <= 8, meas3DistDir = ''; end
    
    %% Initialization   
    % Initialize other variables
    name = actBody.name;
    setupDefault = struct('label', 'ekfv3', 'est', 'ekfv3', ...
        'accData', 'w__v', 'accDataNoise', 0.0, 'oriData', 'w__v', ...
        'initSrc', 'w__v', 'stepDetection', 'av01', ...
        'stepDetectWindow', 0.25, ...
        'stepDetectThreshold', struct('MP', 0, 'LA', 1, 'RA', 1, 'LF', 1, 'RF', 1), ...
        'applyPred', 0, 'applyMeas', 0, 'applyCstr', 0, 'accPelvBias', false, ...
        'sigma2DistLeg', 1e0, 'sigma2DistLARA', 1e0, 'c3d', false, ...
        'sigmaQAcc', 0.5, 'P', 100, ...
        'sigma2QGyr', 1, ...
        'sim3Dist', 0, 'sim3DistSigma', 0.0);
    
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
        % load setup instances to cs with some preprocessing
        cs = setups{sI};
        if ~(cs.accData(end) == 'v')
            cs.accData = strcat(cs.accData, cs.initSrc(end));
        end
        if ~(cs.oriData(end) == 'v')
            cs.oriData = strcat(cs.oriData, cs.initSrc(end));
        end
        
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
        
        %% measurement of 3 distance (pelvis to ankles, ankle to ankle)
        %  Simulate distance measurement by generating pairwise combinations
        if (cs.sim3Dist > 0) && ~(strcmp(meas3DistDir, ''))
            if ~exist(meas3DistDir, 'dir')
                mkdir(meas3DistDir);
            end
            
            meas3DistFname = sprintf("%s/%s-I%s-%.3f-%03d.csv", ...
                                     meas3DistDir, name, cs.initSrc,  ...
                                     cs.sim3DistSigma, mod(cs.sim3Dist, 10));
            if ~exist(meas3DistFname, 'file')
                meas3Dist = mocapdb.dist.simulateDistMeas(actBody, ...
                                        cs.sim3DistSigma, meas3DistFname);
            else
                meas3Dist = mocapdb.dist.loadCSV(meas3DistFname, ...
                                [actBody.ftStartIndex, actBody.ftEndIndex]);
            end
        else
            meas3Dist = false;
        end
        
        %% big switch statement for different algorithms
        t0 = cputime;
        
        try
            if strncmp(cs.est, 'ckf', 3)
                % calculate initial state
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

                v3Options = struct('fs', actBody.fs, 'applyMeas', cs.applyMeas, ...
                    'applyCstr', cs.applyCstr, 'sigma2QAccMP', cs.sigma2QAcc, ...
                    'sigma2QAccLA', cs.sigma2QAcc, 'sigma2QAccRA', cs.sigma2QAcc, ...
                    'sigma2DistLLeg', cs.sigma2DistLeg, 'sigma2DistRLeg', cs.sigma2DistLeg, ...
                    'sigma2DistLARA', cs.sigma2DistLARA, ...
                    'alphaLKmin', alphaLKmin, 'alphaRKmin', alphaRKmin);
                if isfield(cs, 'sigma2RXYPosPVLSRS')
                    v3Options.sigma2RPosMPLARA = double(cs.sigma2RXYPosPVLSRS);
                end
                if isfield(cs, 'sigma2RLim')
                    v3Options.sigma2RPosMPLimit = double(cs.sigma2RLim);
                    v3Options.sigma2RPosLALimit = double(cs.sigma2RLim);
                    v3Options.sigma2RPosRALimit = double(cs.sigma2RLim);
                end
                if strcmp(cs.est, 'ckf')
                    [~, x_pos_v2, t_dat_v2 ] = papers.ckf2019.ckf_3imus(x0, cs.P, ...
                        gfrAcc.(cs.accData).Pelvis, isStep.MIDPEL, imu.(cs.oriData).Pelvis.ori, ...
                        gfrAcc.(cs.accData).L_LowLeg, isStep.LTIO, imu.(cs.oriData).L_LowLeg.ori, ...
                        gfrAcc.(cs.accData).R_LowLeg, isStep.RTIO, imu.(cs.oriData).R_LowLeg.ori, ...
                        body.PV_d, body.LT_d, body.RT_d, body.LS_d, body.RS_d, ...
                        v3Options);
                elseif strcmp(cs.est, 'ckfdist')
                    [~, x_pos_v2, t_dat_v2 ] = papers.ckfdist2020.ckfdist_3imus(x0, cs.P, ...
                        gfrAcc.(cs.accData).Pelvis, isStep.MIDPEL, imu.(cs.oriData).Pelvis.ori, ...
                        gfrAcc.(cs.accData).L_LowLeg, isStep.LTIO, imu.(cs.oriData).L_LowLeg.ori, ...
                        gfrAcc.(cs.accData).R_LowLeg, isStep.RTIO, imu.(cs.oriData).R_LowLeg.ori, ...
                        body, meas3Dist, v3Options);
                else
                    error("est mode %s not supported", cs.est);
                end

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

                if cs.c3d
                    if strcmp(cs.est, 'ckfdist')
                        acq = papers.ckfdist2020.exportc3d( ...
                            x_pos_v2, t_dat_v2, ...
                            sprintf("%s/%s-%s", savedir, name, cs.label), ...
                            actBody, isStep, imu.raw);
                    else
                        estBody.exportc3d(sprintf("%s/%s-%s.c3d", savedir, name, cs.label), ...
                            struct(), actBody, isStep.LTIO, isStep.RTIO);
                    end
                end
            elseif strncmp(cs.est, 'lieekf', 6)
                tmp = struct();
                % a bit unintuitive but all I am doing here is change the key
                % from Pelvis, L_LowLeg, R_LowLeg to PV, LS, RS because that
                % is what the algo wants
                for i=[struct('xsn', 'Pelvis', 'key', 'PV'), ...
                       struct('xsn', 'L_LowLeg', 'key', 'LS'), ...
                       struct('xsn', 'R_LowLeg', 'key', 'RS')]
                    for j=["acc", "ori", "gyr"]
                       tmp.(j).(i.key) = imu.(cs.oriData).(i.xsn).(j);
                    end
                end
                tmp.step.PV = isStep.MIDPEL;
                tmp.step.LS = isStep.LTIO;
                tmp.step.RS = isStep.RTIO;

                % calculate initial state
                vel = actBody.calcJointVel({'MIDPEL', 'LTIO', 'RTIO'});                          
                x0 = struct( ...
                        'W_T_PV', [quat2rotm(actBody.qRPV(1,:)) actBody.MIDPEL(1,:)'; ...
                                   zeros(1,3) 1], ...
                        'W_T_LS', [quat2rotm(actBody.qLSK(1,:)) actBody.LTIO(1,:)'; ...
                                   zeros(1,3) 1], ...
                        'W_T_RS', [quat2rotm(actBody.qRSK(1,:)) actBody.RTIO(1,:)'; ...
                                   zeros(1,3) 1], ...
                        'vec', [vel.MIDPEL(1,:), vel.LTIO(1,:), vel.RTIO(1,:)]' );

                % Backtrack x0 such that x_pred at t=1 is x0 from viconBody
                dt = 1.0/actBody.fs; dt2 = dt*dt;
                x0.vec = x0.vec - [gfrAcc.(cs.accData).Pelvis(1,:) ...
                                   gfrAcc.(cs.accData).L_LowLeg(1,:) ...
                                   gfrAcc.(cs.accData).R_LowLeg(1,:) ]' *dt;
                x0.W_T_PV(1:3,4) = x0.W_T_PV(1:3,4) - x0.vec(1:3,:)*dt ...
                                - 0.5*gfrAcc.(cs.accData).Pelvis(1,:)'*dt2;
                x0.W_T_LS(1:3,4) = x0.W_T_LS(1:3,4) - x0.vec(4:6,:)*dt ...
                                - 0.5*gfrAcc.(cs.accData).L_LowLeg(1,:)'*dt2;
                x0.W_T_RS(1:3,4) = x0.W_T_RS(1:3,4) - x0.vec(7:9,:)*dt ...
                                - 0.5*gfrAcc.(cs.accData).R_LowLeg(1,:)'*dt2;

                v3Options = struct('fs', actBody.fs, 'applyPred', cs.applyPred, ...
                    'applyMeas', cs.applyMeas, 'applyCstr', cs.applyCstr, ...
                    'sigma2DistLLeg', cs.sigma2DistLeg, 'sigma2DistRLeg', cs.sigma2DistLeg, ...
                    'sigma2DistLARA', cs.sigma2DistLARA, ...
                    'sigma2QAccPV', cs.sigma2QAcc, 'sigma2QAccLS', cs.sigma2QAcc, ...
                    'sigma2QAccRS', cs.sigma2QAcc, 'sigma2QGyrPV', cs.sigma2QGyr, ...
                    'sigma2QGyrLS', cs.sigma2QGyr, 'sigma2QGyrRS', cs.sigma2QGyr );
                if isfield(cs, 'sigma2RXYPosPVLSRS')
                    v3Options.sigma2RXYPosPVLSRS = double(cs.sigma2RXYPosPVLSRS);
                end
                if isfield(cs, 'sigma2RLim')
                    v3Options.sigma2RLimPos = double(cs.sigma2RLim);
                    v3Options.sigma2RLimOri = double(cs.sigma2RLim);
                end
                if isfield(cs, 'sigma2ROri')
                    v3Options.sigma2ROriPV = double(cs.sigma2ROri);
                    v3Options.sigma2ROriLS = double(cs.sigma2ROri);
                    v3Options.sigma2ROriRS = double(cs.sigma2ROri);
                end

                if strncmp(cs.est, 'lieekfdist', 10)
                    [ x_pos_v2, t_dat_v2 ] = papers.lgkf5seg.lieekfdist_3_kmus( ...
                            x0, cs.P, tmp.acc, tmp.step, tmp.ori, tmp.gyr, ...
                            body, meas3Dist, v3Options);
                else
                    [ x_pos_v2, t_dat_v2 ] = papers.lgkf5seg.lieekf_3_kmus( ...
                            x0, cs.P, tmp.acc, tmp.step, tmp.ori, tmp.gyr, ...
                            body, v3Options);
                end

                estBody = papers.lgkf5seg.exportGrBody(x_pos_v2, t_dat_v2);
                if cs.c3d
                    acq = papers.lgkf5seg.exportc3d( ...
                        x_pos_v2, t_dat_v2, ...
                        sprintf("%s/%s-%s", savedir, name, cs.label), ...
                        actBody, isStep, imu.raw);
                end
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

                if cs.c3d
                    estBody.exportc3d(sprintf("%s/%s-%s.c3d", savedir, name, cs.label), ...
                        struct(), actBody);
                end
            end

            runtime = cputime-t0;

            % save result
            if ~strcmp(savedir, '')
                save(sprintf("%s/%s-%s.mat", savedir, name, cs.label), ...
                     'estBody', 'estState', 'estState2', 'runtime', 'cs')
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
            if ~strcmp(stepParamDir, '')
                if ~exist(stepParamDir, 'dir'), mkdir(stepParamDir); end
                estBody.dumpStepParams(actBody, isStep, ...
                    sprintf('%s/%s-%s-stepParam', stepParamDir, name, cs.label));
            end
        
        catch e
            runtime = cputime-t0;
            fprintf("Error %s!\n", e.message);
            
            % calculate kinematic gait parameters against actBody (vicon)
            r_buf = estBody2.diffRMSEandMean(nan);
            % calculate some spatiotemporal parameters
            r_bufW = estBody.diffRMSEandMean(nan);
            r_buf = estBody.calcTTDandStepParams(nan, [], r_buf);
        end
        
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
        r_buf.sim3Dist = cs.sim3Dist;
        r_buf.sim3DistSigma = cs.sim3DistSigma;
        
        results(resultsIdx) = r_buf;
        fprintf("Index %3d/%3d: Running time: %.4f\n", resultsIdx, setupN, cputime-t0);
        resultsIdx = resultsIdx + 1;
    end
end
