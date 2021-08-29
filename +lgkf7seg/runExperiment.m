function results = runExperiment(imu, gfrAcc, step, actBody0, dataX, ...
                        setups, savedir, stepParamDir, meas3DistDir)
    % Run experiment on an instance of Vicon and Xsens dataset taken at NeuRA institute
    % 
    % Setup parameters:
    % - label: data instance name (e.g. s1-acting1 or s2-walking1)
    % - est: filter type to be used.
    %      - ckf: papers.ckf2019.ckf_3imus
    %      - ckfdist: papers.ckfdist2020.ckfdist_3imus
    %      - lieekf: papers.lgkf5seg.lieekf_3_kmus
    %      - vsxsens: vs BVH output from xsens (5 seg)
    %      - vsxsens7s: vs BVH output from xsens (7 seg)
    %      - lgckf7stwoimuXY: lie group 7 segment using 2 IMU
    %           X (1st bit): is 0 pitch and roll PV ori
    %             (2nd bit): is PV acc mean of foot acc
    %             (3rd bit): is no ori input
    %           Y: balance point
    % - accData: acceleration data to be used
    %      - w__v: vicon (world frame)
    %      - w__s: sparse (world frame)
    % - oriData: orientation data to be used
    %      - w__v: vicon (world frame)
    %      - w__s: sparse (world frame)
    % - stepDetection: step detection algorithm to be used
    %      - false: turn off
    %      - av01: fixed acceleration variance on accData (var = 1)
    %      - av02: fixed acceleration variance on vicon accData (var = 1)
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
    validateattributes(actBody0, {'pelib.grBody'}, {});
    validateattributes(dataX, {'mocapdb.BVHBody', 'numeric'}, {});
    
    if nargin <= 6, savedir = ''; end
    if nargin <= 7, stepParamDir = ''; end
    if nargin <= 8, meas3DistDir = ''; end
    
    %% Initialization   
    % Initialize other variables
    name = actBody0.name;
    setupDefault = struct('label', 'ekfv3', 'est', 'ekfv3', ...
        'accData', 'w__v', 'accDataNoise', 0.0, 'oriData', 'w__v', ...
        'initSrc', 'w__v', 'debug', false, 'stepDetection', 'av01', ...
        'stepDetectWindow', 0.25, ...
        'stepDetectThreshold', struct('Pelvis', 0, 'L_LowLeg', 1, ...
                            'R_LowLeg', 1, 'L_Foot', 1, 'R_Foot', 1), ...
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
    printIdx = 1;
    body = struct('PV_d', actBody0.calcPelvisLength(), ...
                  'LT_d', actBody0.calcLFemurLength(), ...
                  'RT_d', actBody0.calcRFemurLength(), ...
                  'LS_d', actBody0.calcLShankLength(), ...
                  'RS_d', actBody0.calcRShankLength(), ...
                  'LF_d', actBody0.calcLFootLength(), ...
                  'RF_d', actBody0.calcRFootLength() );

    for sI=1:setupN
        actBody = actBody0.copy();
        
        % load setup instances to cs with some preprocessing
        cs = setups{sI};
        if ~(cs.accData(end) == 'v')
            cs.accData = strcat(cs.accData, cs.initSrc(end));
        end
        if ~(cs.oriData(end) == 'v')
            cs.oriData = strcat(cs.oriData, cs.initSrc(end));
        end
        
        % step detection
        tmp.mapList = [struct('xsn', 'Pelvis', 'key', 'MIDPEL'), ...
                       struct('xsn', 'L_LowLeg', 'key', 'LTIO'), ...
                       struct('xsn', 'R_LowLeg', 'key', 'RTIO'), ...
                       struct('xsn', 'L_Foot', 'key', 'LFT'), ...
                       struct('xsn', 'R_Foot', 'key', 'RFT')];
        if strcmp(cs.stepDetection, 'av01')
            VAR_WIN  = floor(actBody0.fs*cs.stepDetectWindow); % NUM_SAMPLES
            ACC_VAR_THRESH = cs.stepDetectThreshold;

            for i = tmp.mapList
                movVarAccBuf = movingvar(sqrt( sum(gfrAcc.(cs.accData).(i.xsn) .^2, 2)), VAR_WIN);
                isStep.(i.key) = movVarAccBuf < ACC_VAR_THRESH.(i.xsn);
            end
        elseif strcmp(cs.stepDetection, 'av02')            
            VAR_WIN  = floor(actBody0.fs*cs.stepDetectWindow); % NUM_SAMPLES
            ACC_VAR_THRESH = cs.stepDetectThreshold;

            for i = tmp.mapList
                movVarAccBuf = movingvar(sqrt( sum(gfrAcc.w__v.(i.xsn) .^2, 2)), VAR_WIN);
                isStep.(i.key) = movVarAccBuf < ACC_VAR_THRESH.(i.xsn);
            end
        elseif strcmp(cs.stepDetection, 'av03')
            isStep = struct('MIDPEL', false(actBody.nSamples, 1), ...
                'LTIO', logical(step.stepL), 'RTIO', logical(step.stepR), ...
                'LFT', logical(step.stepL), 'RFT', logical(step.stepR) );
        else
            isStep = struct('MIDPEL', false(actBody.nSamples, 1), ...
                'LTIO', false(actBody.nSamples, 1), 'RTIO', false(actBody.nSamples, 1), ...
                'LFT', false(actBody.nSamples, 1), 'RFT', false(actBody.nSamples, 1) );
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
        
        if strncmp(cs.est, 'ckf', 3)
            % target segments for calculating dOri
            targetSeg = {'qLTH', 'qRTH'};
            % remove foot data from vicon since we're only tracking 5 segs
            actBody.LTOE = []; actBody.RTOE = [];
            actBody.qLFT = []; actBody.qRFT = [];

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
            % target segments for calculating dOri
            targetSeg = {'qLTH', 'qRTH'};
            % remove foot data from vicon since we're only tracking 5 segs
            actBody.LTOE = []; actBody.RTOE = [];
            actBody.qLFT = []; actBody.qRFT = [];
            
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
%             else
%                 [ x_pos_v2, t_dat_v2 ] = papers.lgkf5seg.lieekf_3_kmus( ...
%                         x0, cs.P, tmp.acc, tmp.step, tmp.ori, tmp.gyr, ...
%                         body, v3Options);
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
        elseif strncmp(cs.est, 'lgckf7s', 7)
            tmp = struct('LF_p_LA', cs.LFA_p_LF', 'RF_p_RA', cs.RFA_p_RF');
            tmp.LF_p_LA(3) = -tmp.LF_p_LA(3);
            tmp.RF_p_RA(3) = -tmp.RF_p_RA(3);
            v3Options = struct('fs', actBody.fs, 'applyPred', cs.applyPred, ...
                'applyMeas', cs.applyMeas, 'applyCstr', cs.applyCstr, ...
                'sigma2QAccPV', cs.sigma2QAcc, 'sigma2QAccLF', cs.sigma2QAcc, ...
                'sigma2QAccRF', cs.sigma2QAcc, 'sigma2QGyrPV', cs.sigma2QGyr, ...
                'sigma2QGyrLF', cs.sigma2QGyr, 'sigma2QGyrRF', cs.sigma2QGyr, ...
                'LF_p_LA', tmp.LF_p_LA, 'RF_p_RA', tmp.RF_p_RA, ...
                'fulldebug', cs.debug);
            
            if isfield(cs, 'sigma2RXYPosPVLSRS')
                v3Options.sigma2RxyposPVLARA = double(cs.sigma2RXYPosPVLSRS);
            end
            if isfield(cs, 'sigma2RLim')
                v3Options.sigma2RLimPos = double(cs.sigma2RLim);
                v3Options.sigma2RLimOri = double(cs.sigma2RLim);
            end
            if isfield(cs, 'sigma2ROri')
                v3Options.sigma2ROriPV = double(cs.sigma2ROri);
                v3Options.sigma2ROriLF = double(cs.sigma2ROri);
                v3Options.sigma2ROriRF = double(cs.sigma2ROri);
            end
            if isfield(cs, 'sigma2RZPosPV')
                v3Options.sigma2RZPosPV = double(cs.sigma2RZPosPV);
            end
            if isfield(cs, 'sigma2RYawAve')
                v3Options.sigma2Ryawave = double(cs.sigma2RYawAve);
                v3Options.sigma2Rpelvpitchroll = double(cs.sigma2RYawAve);
            end
            
            if strncmp(cs.est, 'lgckf7stwoimu', 13)
                targetSeg = {'qRPV', 'qLTH', 'qRTH', 'qLSK', 'qRSK'};
                tmp.mapList = [struct('xsn', 'L_Foot', 'key', 'LF'), ...
                               struct('xsn', 'R_Foot', 'key', 'RF')];
                tmp.ori.PV = imu.(cs.oriData).calcSegMeanRot({'L_Foot', 'R_Foot'}, -pi/2);

                tmp.mode = struct('val', str2double(cs.est(14)));
                tmp.mode.isFlatOri = bitand(tmp.mode.val, 1);
                tmp.mode.isAccMean = bitand(tmp.mode.val, 2);
                tmp.mode.isNoOri = bitand(tmp.mode.val, 4);
                
                if tmp.mode.isNoOri
                    tmp.ori.PV = nan(actBody.nSamples, 4);
                elseif tmp.mode.isFlatOri
                    % lgckf7sv2twoimu
                    tmp.pvOriEuler = quat2eul(tmp.ori.PV, 'ZYX');
                    tmp.pvOriYaw = zeros(actBody.nSamples,4);
                    tmp.pvOriYaw(:,3) = 1; tmp.pvOriYaw(:,4) = tmp.pvOriEuler(:,1);
                    tmp.ori.PV = axang2quat(tmp.pvOriYaw);
                else
                end
                
                if tmp.mode.isAccMean
                    tmp.buf = (gfrAcc.(cs.accData).L_Foot + gfrAcc.(cs.accData).R_Foot)/2 + [0 0 9.81];
                    tmp.W_PVacc0 = tmp.buf(1,:) - [0 0 9.81];
                    tmp.W_acc.PV = tmp.buf;
                else % lgckf7sv2twoimu
                    % tmp.acc.PV = quatrotate(tmp.ori.PV, repmat([0 0 9.81], actBody.nSamples, 1));
                    tmp.W_acc.PV = repmat([0 0 9.81], actBody.nSamples, 1);
                    tmp.W_PVacc0 = [0 0 0];
                end
                tmp.acc.PV = quatrotate(tmp.ori.PV, tmp.W_acc.PV);
                tmp.gyr.PV = zeros(actBody.nSamples,3);
%                 v3Options.sigma2RxyposPVLARA = 1e0;
            else
                % target segments for calculating dOri
                targetSeg = {'qLTH', 'qRTH', 'qLSK', 'qRSK'};
            
                tmp.mapList = [struct('xsn', 'Pelvis', 'key', 'PV'), ...
                               struct('xsn', 'L_Foot', 'key', 'LF'), ...
                               struct('xsn', 'R_Foot', 'key', 'RF')];
            end
            % a bit unintuitive but all I am doing here is change the key
            % from Pelvis, L_Foot, R_Foot to PV, LF, RF because that
            % is what the algo wants
            for i=tmp.mapList
                for j=["acc", "ori", "gyr"]
                   tmp.(j).(i.key) = imu.(cs.oriData).(i.xsn).(j);
                end
                tmp.W_acc.(i.key) = quatrotate(quatconj(tmp.ori.(i.key)), tmp.acc.(i.key));
            end
            
            tmp.step.PV = isStep.MIDPEL;
            tmp.step.LF = isStep.LFT;
            tmp.step.RF = isStep.RFT;
            
            
            % setting up left and right balance point
            if strcmp(cs.est(end), 'b') % LBP = foot center
                v3Options.LF_p_LBP = [0; 0; 0; 1];
                v3Options.RF_p_RBP = [0; 0; 0; 1];
            elseif strcmp(cs.est(end), 'c')
                tmp.alpha = -0.105760597323741; % obtained experimentally from neura data
                v3Options.LF_p_LBP = [0; 0; tmp.alpha*body.LF_d; 1];
                v3Options.RF_p_RBP = [0; 0; tmp.alpha*body.RF_d; 1];
            else % LBP = ankles
                v3Options.LF_p_LBP = v3Options.LF_p_LA; % left balance point
                v3Options.RF_p_RBP = v3Options.RF_p_LA; % right balance point
            end
                
            % calculate initial state
            % LFA and RFA = LF and RF frame centered at ankle
            tmp.LFA_T_LF = eye(4,4); tmp.LFA_T_LF(1:4,4) = cs.LFA_p_LF;
            tmp.RFA_T_RF = eye(4,4); tmp.RFA_T_RF(1:4,4) = cs.RFA_p_RF;
            
            vel = actBody.calcJointVel(struct('MIDPEL', [0 0 0 1], ...
                                              'LFT', cs.LFA_p_LF, ...
                                              'RFT', cs.RFA_p_RF));                          
            x0 = struct( ...
                    'W_T_PV', actBody.getTMatrix('MIDPEL', 1), ...
                    'W_T_LF', actBody.getTMatrix('LFT', 1) * tmp.LFA_T_LF, ...
                    'W_T_RF', actBody.getTMatrix('RFT', 1) * tmp.RFA_T_RF, ...
                    'vec', [vel.MIDPEL(1,:)'; vel.LFT(1,:)'; vel.RFT(1,:)' ] );

            % Backtrack x0 such that x_pred at t=1 is x0 from viconBody
            dt = 1.0/actBody.fs; dt2 = dt*dt;
            if strncmp(cs.est, 'lgckf7stwoimu', 13)
                x0.vec = x0.vec - [tmp.W_PVacc0 ...
                                   gfrAcc.(cs.accData).L_Foot(1,:) ...
                                   gfrAcc.(cs.accData).R_Foot(1,:) ]' *dt;
                x0.W_T_PV(1:3,4) = x0.W_T_PV(1:3,4) - x0.vec(1:3,:)*dt ...
                                - 0.5*tmp.W_PVacc0'*dt2;
            else
                x0.vec = x0.vec - [gfrAcc.(cs.accData).Pelvis(1,:) ...
                                   gfrAcc.(cs.accData).L_Foot(1,:) ...
                                   gfrAcc.(cs.accData).R_Foot(1,:) ]' *dt;
                x0.W_T_PV(1:3,4) = x0.W_T_PV(1:3,4) - x0.vec(1:3,:)*dt ...
                                - 0.5*gfrAcc.(cs.accData).Pelvis(1,:)'*dt2;

            end
            x0.W_T_LF(1:3,4) = x0.W_T_LF(1:3,4) - x0.vec(4:6,:)*dt ...
                            - 0.5*gfrAcc.(cs.accData).L_Foot(1,:)'*dt2;
            x0.W_T_RF(1:3,4) = x0.W_T_RF(1:3,4) - x0.vec(7:9,:)*dt ...
                            - 0.5*gfrAcc.(cs.accData).R_Foot(1,:)'*dt2;
                        
            [ x_pos_v2, t_dat_v2 ] = papers.lgkf7seg.lgcekf7seg_3imus_v3( ...
                    x0, cs.P, tmp.W_acc, tmp.step, tmp.ori, tmp.gyr, ...
                    body, meas3Dist, v3Options);
            
            estBody = papers.lgkf7seg.exportGrBody(x_pos_v2, t_dat_v2);
            if cs.c3d
                if strncmp(cs.est, 'lgckf7stwoimu', 13)
                    tmp.simPV = struct('acc', tmp.acc.PV, 'gyr', tmp.gyr.PV, ...
                                       'ori', tmp.ori.PV);
                    if cs.debug
                        [acq, acqDebug] = papers.lgkf7seg.exportc3d( ...
                            x_pos_v2, t_dat_v2, ...
                            sprintf("%s/%s-%s", savedir, name, cs.label), ...
                            actBody, isStep, imu.raw, -1, tmp.simPV);
                    else
                        acq = papers.lgkf7seg.exportc3d( ...
                            x_pos_v2, t_dat_v2, ...
                            sprintf("%s/%s-%s", savedir, name, cs.label), ...
                            actBody, isStep, imu.raw, -1, tmp.simPV);
                    end
                else
                    if cs.debug
                        [acq, acqDebug] = papers.lgkf7seg.exportc3d( ...
                            x_pos_v2, t_dat_v2, ...
                            sprintf("%s/%s-%s", savedir, name, cs.label), ...
                            actBody, isStep, imu.raw);
                    else
                        acq = papers.lgkf7seg.exportc3d( ...
                            x_pos_v2, t_dat_v2, ...
                            sprintf("%s/%s-%s", savedir, name, cs.label), ...
                            actBody, isStep, imu.raw);
                    end
                end
            end
            estState = x_pos_v2;
            estState2 = t_dat_v2;
        elseif strncmp(cs.est, 'vsxsens', 7)
            if isempty(dataX)
                fprintf("Index %3d/%3d: Input dataX (xsens body) is empty. Skipping %s\n", printIdx, setupN, cs.est);
                printIdx = printIdx + 1;
                continue
            end
            estBody = dataX.togrBody(1:actBody.nSamples, {'name', 'act', 'oriUnit', 'deg', ...
                         'lnSymbol', '-', 'ptSymbol', '*', 'fs', actBody.fs, ...
                         'xyzColor', {'m', 'y', 'c'}});
                     
            if strcmp(cs.est, 'vsxsens7s')
                % target segments for calculating dOri
                targetSeg = {'qLTH', 'qRTH', 'qLSK', 'qRSK'};
            else
                % target segments for calculating dOri
                targetSeg = {'qLTH', 'qRTH'};
                % remove foot data from vicon since we're only tracking 5 segs
                actBody.LTOE = []; actBody.RTOE = [];
                actBody.qLFT = []; actBody.qRFT = [];

                estBody.LTOE = []; estBody.RTOE = [];
                estBody.qLFT = []; estBody.qRFT = [];
            end
            
            estState = [];
            estState2 = [];
            
            if cs.c3d
                estBody.exportc3d(sprintf("%s/%s-%s.c3d", savedir, name, cs.label), ...
                    struct(), actBody);
            end
        end

        runtime = cputime-t0;
        fprintf("Index %3d/%3d: Running time: %.4f\n", printIdx, setupN, cputime-t0);
        
        t1 = cputime;
        % save result
        if ~strcmp(savedir, '')
            save(sprintf("%s/%s-%s.mat", savedir, name, cs.label), ...
                 'estBody', 'estState', 'estState2', 'runtime', 'cs')
        end
        
        %% evaluation
        % convert grBody to midpel frame
        estBodyRel = estBody.changeRefFrame('MIDPEL');
        % revert back to global frame but with MIDPEL position equal to actBody
        estBody2 = estBodyRel.toWorldFrame(actBody.MIDPEL, estBody.qRPV);

        % calculate kinematic gait parameters against actBody (vicon)
        r_buf = estBody2.diffRMSEandMean(actBody, true, targetSeg);
        
        % calculate some spatiotemporal parameters
        r_bufW = estBody.diffRMSEandMean(actBody, true, targetSeg);
        r_buf = estBody.calcTTDandStepParams(actBody, isStep, r_buf);
        for i=["MIDPEL", "LTIO", "RTIO"]
            for j=["RMSE", "Std", "Mean"]
                r_buf.(sprintf("%sW%s",i,j)) = r_bufW.(sprintf("%s%s",i,j));
            end
        end
        r_buf.dPosW = r_bufW.dPos;   
        if ~strcmp(stepParamDir, '')
            if ~exist(stepParamDir, 'dir'), mkdir(stepParamDir); end
            estBody.dumpStepParams(actBody, isStep, ...
                sprintf('%s/%s-%s-stepParam', stepParamDir, name, cs.label));
        end
        % labeling
        r_buf.name = name;
        r_buf.label = cs.label;
        r_buf.runtime = runtime;
        r_buf.sim3Dist = cs.sim3Dist;
        r_buf.sim3DistSigma = cs.sim3DistSigma;
        
        results(resultsIdx) = r_buf;
        fprintf("Index %3d/%3d: Evaluation time: %.4f\n", printIdx, setupN, cputime-t1);
        resultsIdx = resultsIdx + 1;
        printIdx = printIdx + 1;
    end
end
