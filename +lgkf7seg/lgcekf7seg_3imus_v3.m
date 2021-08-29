function [ xtilde, debug_dat ] = lgcekf7seg_3imus_v3(x0, P0, ...
    W_a_, step, W_R_, B_w_, body, uwb_mea, options)
% Constrained EKF using Lie group/algebra representation
% Representation: SE(3)^3 + R^9
% 3 IMUs presumably worn on the body in the following configuration: 
% mid pelvis, left foot, right foot
% 
% Coding notation based on <a href="http://paulfurgale.info/news/2014/6/9/representing-robot-pose-the-good-the-bad-and-the-ugly">link</a>
%
% More detail about options
%      fs: sampling frequency of the magnetic and inertial measurement units
%      applyPred: 3 digit ZYX
%          X:   1st bit use velocity to predict position e
%               2nd bit use angular velocity to predict orientation
%          Y:   1st bit clip pelvis X and Y orientation Q growth
%          Z:   1st bit track angular velocity (not yet implemented)
%               2nd bit track acc bias in sensor frame
%               3rd bit track ang bias in sensor frame (not yet implemented)
%      applyMeas: 4 digit AZYX
%          X:   1st bit Orientation update (PV, LF, RF)
%               2nd bit Zupt and Ankle zpos
%               3rd bit Track bias
%          Y:   1st bit Pelvis xy=ankle average
%               2nd bit Pelvis z=initial height
%          Z:   1st bit covariance limiter
%               2nd and 3rd bit: ankle z pos settings
%                   0: z pos = floor, constant R cov matrix
%                   1: z pos = external pos (e.g., guided by smoothed foot
%                   dead reckoning)
%                   2: z pos = floor but with varying R cov matrix 
%                   (constant exp decay)
%                   3: z pos = floor but with varying R cov matrix 
%                   (look ahead and exp decay by next foot step)
%                   4: z pos = floor but with varying R cov matrix 
%                   (look ahead and linear decay by next foot step)
%          A:   1st bit Orientation update (LF and RF only)
%               2nd-3rd bit 
%                   0: off
%                   1: pelvis yaw = ankle yaw average
%                   2: pelvis yaw = ankle yaw average + separate zero pitch
%                   3: pelvis yaw = ankle yaw average with pitch and roll zero
%                   through exp(theta_yaw [0 0 1]'^) formulation
%                      theta yaw ave = (atan(lf) + atan(rf))/2
%                   4: same as 3 but theta yaw ave = atan(lf+rf)
%      applyCstr: 3 digit ZYX
%          X:   1st bit enforce max leg length
%               2nd bit enforce hinge knee and ankle joint
%               3rd bit enforce pelvis ori no pitch and roll
%                       (not 100% working...)
%          Z:   1st bit do P update
%
% :param x0: initial state in the GFR
% :param P0: initial covariance
% :param W_a_: acceleration of PV, LS, RS in the body frame
% :param step: boolean vector of PV, LS, RS indicating step detection
% :param W_R_: PV, LS, RS orientation in the GFR (rotm)
% :param B_w_: PV, LS, RS angular velocity in the body frame
% :param body: Length of PV_d (pelvis), RT_d and LT_d (r/l femur), RS_d and LS_d (r/l tibia)
% :param uwb_mea: a structure containing the range measurements (m) between
% :param options: struct containing the estimator settings:
%
% :returns: xtilde and debug_dat
% .. Author: - Luke Wicent Sy (GSBME, 2019 Oct 31)
    N = {};
    [N.samples, ~] = size(W_a_.PV);
    
    %% input parsing
    fOpt = struct('fs', 100, ...
          'applyPred', 1, 'applyMeas', 1, 'applyCstr', 1, ...
          'sigma2QAccPV', 1, 'sigma2QAccLF', 1, 'sigma2QAccRF', 1, ...
          'sigma2QGyrPV', 1, 'sigma2QGyrLF', 1, 'sigma2QGyrRF', 1, ...
          'sigma2QOriPV', 1e3, 'sigma2QOriLF', 1e3, 'sigma2QOriRF', 1e3, ...
          'sigma2QClipPVXYOri', 0, ...
          'sigma2QVelMP', 1, 'sigma2QVelLA', 1, 'sigma2QVelRA', 1, ...
          'sigma2QAccBMP', 1, 'sigma2QAccBLF', 1, 'sigma2QAccBRF', 1, ...
          'sigma2ROriPV', 1e1, 'sigma2ROriLF', 1e1, 'sigma2ROriRF', 1e1, ...
          'sigma2RZPosPV', 1e0, 'sigma2RZPosLF', 1e1, 'sigma2RZPosRF', 1e1, ...
          'sigma2RZPosLFStep', 1e-4, 'sigma2RZPosRFStep', 1e-4, ...
          'decaytZPosLF', 1, 'decaytZPosRF', 1, ... 
          'sigma2RZuptLF', 1e-2, 'sigma2RZuptRF', 1e-2, ...
          'sigma2RAccBiasLF', 1e-2, 'sigma2RAccBiasRF', 1e-2, ...
          'sigma2RxyposPVLARA', 1e0, 'sigma2Ryawave', 1e-1, ...
          'sigma2Rpelvpitchroll', 1e-1, ...
          'sigma2RLimPos', 1e1, 'sigma2RLimOri', 1e1, ...
          'yZPosLF', nan, 'yZPosRF', nan, ...
          'LF_p_LBP', [0; 0; -0.5*body.LF_d; 1], ... % left balance point
          'RF_p_RBP', [0; 0; -0.5*body.RF_d; 1], ... % right balance point
          'LF_p_LA', [0; 0; -0.5*body.LF_d; 1], ...
          'RF_p_RA', [0; 0; -0.5*body.RF_d; 1], ...
          'fulldebug', false);
    optionFieldNames = fieldnames(options);
    for i=1:length(optionFieldNames)
        if ~isfield(fOpt, optionFieldNames{i})
            error("Field %s is not a valid option", optionFieldNames{i});
        end
        fOpt.(optionFieldNames{i}) = options.(optionFieldNames{i});
    end
    
    %% Feature on off
    knob = struct();
    for i=["applyPred", "applyMeas", "applyCstr"]
        buf = 1;
        for j=["One", "Ten", "Hun", "Tho", "TTh"]
            sname = sprintf('a%c%sDig', i{1}(6), j);
            knob.(sname) = mod(idivide(int32(fOpt.(i)), buf, 'floor'), 10);
            buf = buf * 10;
        end
    end
    
    % Prediction knobs
    knob.pred = struct('any', fOpt.applyPred > 0, ...
        'posFromVel', bitand(knob.aPOneDig, 1), ...
        'oriFromAVel', bitand(knob.aPOneDig, 2), ...
        'trackAVel', bitand(knob.aPHunDig, 1), ...
        'trackAccBias', bitand(knob.aPHunDig, 2) > 0, ...
        'trackAVelBias', bitand(knob.aPHunDig, 4) > 0, ...
        'clipPVXYOri', bitand(knob.aPTenDig, 1) );
    
    % Measurement knobs
    knob.meas = struct('ori', ...
        struct('LF', bitand(knob.aMOneDig, 1) | bitand(knob.aMThoDig, 1), ...
               'RF', bitand(knob.aMOneDig, 1) | bitand(knob.aMThoDig, 1), ...
               'PV', bitand(knob.aMOneDig, 1) > 0), ...
        'zpos', struct('LF', false, 'RF', false, 'PV', false), ...
        'bias', bitand(knob.aMOneDig, 4)>0, ...
        'zupt', false, 'xyposPVLARA', bitand(knob.aMTenDig, 1), ...
        'covlim', bitand(knob.aMHunDig, 1), ...
        'zposmode', idivide(knob.aMHunDig,2), ...
        'yawavemode', idivide(knob.aMThoDig,2) );

    if bitand(knob.aMOneDig, 2)
        knob.meas.zupt = true;
        knob.meas.zpos.LF = true;
        knob.meas.zpos.RF = true;
    else
        step.LF = false(N.samples, 1);
        step.RF = false(N.samples, 1);
    end
    
    if bitand(knob.aMTenDig, 2)
        knob.meas.zpos.PV = true;
        step.PV = true(N.samples, 1);
    else
        knob.meas.zpos.PV = false;
        step.PV = false(N.samples, 1);
    end
    
    % Constraint knobs
    knob.cstr = struct('any', fOpt.applyCstr > 0, ...
        'maxleglength', bitand(knob.aCOneDig, 1), ...
        'hingekneeankle', bitand(knob.aCOneDig, 2), ...
        'pelvnopitchroll', bitand(knob.aCOneDig, 4), ...
        'Pupdate', bitand(knob.aCHunDig, 1));
    
    %% State and error covariance initialization
    list.se3State = ["W_T_PV", "W_T_LF", "W_T_RF"];
    list.body = ["PV", "LF", "RF"];
    N.list.state = length(list.se3State);
    N.se3State = N.list.state*6;    
    N.r3State = N.list.state*3 + ...
                N.list.state*3*knob.pred.trackAccBias;
    N.state = N.se3State;
    
    % calculate indices
    idx = struct('se3State', 1:N.se3State, ...
                 'W_T_PV', 1:6, 'W_T_LF', 7:12, 'W_T_RF', 13:18, ...
                 'vecState', (1:N.r3State)+N.se3State, ...
                 'vec', struct() );
    for i=1:N.list.state
        idx.(list.body(i)) = (1:6)+(i-1)*6;
    end
    % 'posvel', [1:3 7:9 13:15 (N.se3State+1):N.state], ...
    idx.posvel = [idx.W_T_PV(1:3), idx.W_T_LF(1:3), idx.W_T_RF(1:3), ...
                  idx.vecState];
    
    if knob.pred.trackAccBias, iAccB="accB"; else, iAccB=[]; end
    for i=["vel", iAccB]
        for j=list.body
            k = sprintf('%s%s', i, j);
            idx.(k) = (1:3) + N.state;
            idx.vec.(k) = (1:3) + N.state - N.se3State;
            N.state = N.state + 3;
        end
    end
    
    g = [0 0 9.81]';        % gravity
    dt = 1/(fOpt.fs);       % assume constant sampling interval
    dt2 = 0.5*dt^2;
    I_N = eye(N.state);
    
    list.kfStep = ["hatPri", "hatPos", "tilde"];
    x = struct();
    P = struct();
    for i=list.kfStep
        x.(i) = struct();
        P.(i) = struct();
        for j=list.se3State
            x.(i).(j) = nan(4,4,N.samples+1);
        end
        x.(i).vec = nan(N.r3State,N.samples+1);
        P.(i) = nan(N.state,N.state,N.samples+1);
    end
    
    %% Set initial states (k=1)
    for i=list.se3State
        x.tilde.(i)(:,:,1) = x0.(i);
    end
    x.tilde.vec(:,1) = x0.vec(1:N.r3State);
    if isscalar(P0)
        P0 = I_N*P0;
    end
    if knob.pred.clipPVXYOri
        P0(idx.W_T_PV(4:5), idx.W_T_PV(4:5)) = 100*P0(idx.W_T_PV(4:5), idx.W_T_PV(4:5));
    end
    P.tilde(:,:,1) = P0;
    body.LL_d = body.LT_d + body.LS_d;
    body.RL_d = body.RT_d + body.RS_d;
    
    %% Prepare input
    for i=list.body
        if ndims(W_R_.(i)) == 2 && size(W_R_.(i), 2) == 4 % quaternion rep
            W_R_.(i) = quat2rotm(W_R_.(i));
        end
    end
    
    %% Prediction step initialization
    F.vec0 = eye(N.r3State,N.r3State);
    F.curlyC0 = zeros(N.state,N.state);
    if knob.pred.trackAccBias
        for i=list.body
            sname = sprintf('W_T_%s',i);
            bname = sprintf('accB%s',i);
            F.curlyC0(idx.(sname)(1:3),idx.(bname)) = -dt2*eye(3,3);
        end
    end
    
    if knob.pred.trackAccBias, iAccB=zeros(9,18); else, iAccB=[]; end
    G.vec = [dt*eye(9,9) zeros(9,9); iAccB];
    G.naive = zeros(N.state, 9);
    G.naive(idx.W_T_PV(1:3), 1:3) = dt2*eye(3,3);
    G.naive(idx.W_T_LF(1:3), 4:6) = dt2*eye(3,3);
    G.naive(idx.W_T_RF(1:3), 7:9) = dt2*eye(3,3);
    G.naive(idx.vecState, 1:9) = dt*eye(9,9);
    
%     Q.val0 = diag(repelem([fOpt.sigma2QAccPV, fOpt.sigma2QAccLF, ...
%                       fOpt.sigma2QAccRF, ...
%                       fOpt.sigma2QGyrPV, fOpt.sigma2QGyrLF, ...
%                       fOpt.sigma2QGyrRF], 3));
%     Q.vec = G.vec*Q.val0*G.vec';
    if knob.pred.trackAccBias
        iAccB=[fOpt.sigma2QAccBMP, fOpt.sigma2QAccBLF, fOpt.sigma2QAccBRF]; 
    else
        iAccB=[];
    end
    if knob.pred.oriFromAVel
        Q.comb = diag(repelem([fOpt.sigma2QAccPV*(dt2^2), fOpt.sigma2QGyrPV*(dt^2), ...
                          fOpt.sigma2QAccLF*(dt2^2), fOpt.sigma2QGyrLF*(dt^2), ...
                          fOpt.sigma2QAccRF*(dt2^2), fOpt.sigma2QGyrRF*(dt^2), ...
                          fOpt.sigma2QAccPV*(dt^2), fOpt.sigma2QAccLF*(dt^2), ...
                          fOpt.sigma2QAccRF*(dt^2), iAccB], 3));        
    else
        Q.comb = diag(repelem([fOpt.sigma2QAccPV*(dt2^2), fOpt.sigma2QOriPV, ...
                          fOpt.sigma2QAccLF*(dt2^2), fOpt.sigma2QOriLF, ...
                          fOpt.sigma2QAccRF*(dt2^2), fOpt.sigma2QOriRF, ...
                          fOpt.sigma2QAccPV*(dt^2), fOpt.sigma2QAccLF*(dt^2), ...
                          fOpt.sigma2QAccRF*(dt^2), iAccB], 3));
    end
    Q.comb(idx.posvel, idx.posvel) = G.naive(idx.posvel, :) * ...
                                diag(repelem([fOpt.sigma2QAccPV ...
                                              fOpt.sigma2QAccLF ...
                                              fOpt.sigma2QAccRF], 3)) ...
                                * G.naive(idx.posvel, :)';
    if knob.pred.clipPVXYOri
        Q.comb(idx.W_T_PV(4:5),idx.W_T_PV(4:5)) = diag(repelem(fOpt.sigma2QClipPVXYOri, 2));
    end
    
    idx.u = struct('accPV', 1:3, 'accLF', 4:6, 'accRF', 7:9, ...
                   'avelPV', 10:12, 'avelLF', 13:15, 'avelRF', 16:18);
    u = zeros(18, N.samples);
    
    %% Measurement step initialization
    H = struct('i_x', [1 0 0 0]', 'i_z', [0 0 1 0]', 'i_0', [0 0 0 1]', ...
               'icirc_x', point2fs([1 0 0 0]'), ...
               'icirc_0', point2fs([0 0 0 1]'));
    R = {}; y = {};

    % update orientation measurement
    H.ori = struct();
    for i = list.body
        H.ori.(i) = zeros(3,N.state);
        j = sprintf('W_T_%s', i);
        H.ori.(i)(:, idx.(j)) = [zeros(3,3) eye(3,3)];
        R.ori.(i) = repelem([fOpt.(sprintf('sigma2ROri%s', i))], 3);
        % y.ori is variable and will be generated on the spot
    end
    % Zero velocity and angular velocity update
%     idx.lzupt = [idx.velLS idx.avelLS];
%     idx.rzupt = [idx.velRS idx.avelRS];
%     idx.vec.lzupt = [idx.vec.velLS idx.vec.avelLS];
%     idx.vec.rzupt = [idx.vec.velRS idx.vec.avelRS];
    idx.lzupt = idx.velLF;
    idx.rzupt = idx.velRF;    
    N.zupt = length(idx.lzupt);
    
    H.lzupt = zeros(N.zupt, N.state);
    H.lzupt(:, idx.lzupt) = eye(N.zupt, N.zupt);
    R.lzupt = repelem(fOpt.sigma2RZuptLF, N.zupt);
    y.lzupt = zeros(N.zupt, 1);
        
    H.rzupt = zeros(N.zupt, N.state);
    H.rzupt(:, idx.rzupt) = eye(N.zupt, N.zupt);
    R.rzupt = repelem(fOpt.sigma2RZuptRF, N.zupt);
    y.rzupt = zeros(N.zupt, 1);
    
    % bias
    if knob.meas.bias
        if knob.pred.trackAccBias
            idx.lbias = idx.accBLF;
            idx.rbias = idx.accBRF;
        else
            idx.lbias = [];
            idx.rbias = [];
        end
        N.bias = length(idx.lbias);
               
        H.lbias = zeros(N.bias, N.state);
        H.lbias(:,idx.lbias) = eye(N.bias, N.bias);
        R.lbias = repelem(fOpt.sigma2RAccBiasLF, N.bias);
        
        H.rbias = zeros(N.bias, N.state);
        H.rbias(:,idx.rbias) = eye(N.bias, N.bias);
        R.rbias = repelem(fOpt.sigma2RAccBiasRF, N.bias);
        
        % y.bias is variable and will be generated on the spot
    else
        N.bias = 0;
    end
    
    % zpos = some value assumption (e.g., flat floor)
    H.zpos = {}; y.zpos = {}; 
    R.zpos = struct('PV', repelem(fOpt.sigma2RZPosPV,N.samples,1), ...
                    'LF', repelem(fOpt.sigma2RZPosLFStep,N.samples,1), ...
                    'RF', repelem(fOpt.sigma2RZPosRFStep,N.samples,1) );
    y.zpos.PV = repelem(x.tilde.W_T_PV(3,4,1), N.samples);
    if knob.meas.zposmode == 1
        if isnan(fOpt.yZPosLF), error('fOpt yZPosLF not provided'); 
        else, y.zpos.LF = fOpt.yZPosLF; end
        if isnan(fOpt.yZPosRF), error('fOpt yZPosRF not provided'); 
        else, y.zpos.RF = fOpt.yZPosRF; end
    else
        y.zpos.LF = repelem(min([x.tilde.W_T_LF(3,4,1), ...
                                 x.tilde.W_T_RF(3,4,1)]), N.samples);
        y.zpos.RF = y.zpos.LF;
    end

    if knob.meas.zposmode >= 2
        for i = ["LF", "RF"]
            if ~knob.meas.zpos.(i), continue; end

            events = (step.(i) - step.(i)([1, 1:end-1]));
            decayt = fOpt.(sprintf("decaytZPos%s",i))*fOpt.fs;
            R0 = fOpt.(sprintf("sigma2RZPos%s",i));

            [idxFS, ~] = find(events == 1);
            [idxTF, ~] = find(events == -1);

            if isempty(idxFS) && isempty(idxTF), continue; end

            % if foot strike first, we add an imaginary TF at index 1
            if idxFS(1) < idxTF(1), idxTF = [1; idxTF]; end
            % if toe off last, we add an imaginary FS at index end
            if idxFS(end) < idxTF(end), idxFS = [idxFS; N.samples]; end

            Rbuf = step.(i).*R.zpos.(i);
            for j = 1:length(idxFS)
                sIdx = idxTF(j); eIdx = idxFS(j);
                if knob.meas.zposmode == 2
                    tau = decayt/4;
                    Rbuf(sIdx:eIdx) = exp(-((sIdx:eIdx)-sIdx)/tau)*R0;
                elseif knob.meas.zposmode == 3
                    tau = (eIdx-sIdx+1)/4;
                    Rbuf(sIdx:eIdx) = exp(-((sIdx:eIdx)-sIdx)/tau)*R0;
                elseif knob.meas.zposmode == 4
                    % m = (Rbuf(eIdx)-1) / (decayt);
                    m = (Rbuf(eIdx)-1) / (eIdx-sIdx+1);
                    Rbuf(sIdx:eIdx) = (m*((sIdx:eIdx)-sIdx)+1)*R0;
                end
            end
            R.zpos.(i) = Rbuf;
        end
    end
    % pelvis = ankle x y pos
    H.xyposE = [eye(2,2) zeros(2,2)];
    H.LF_pcirc_LAF = point2fs(fOpt.LF_p_LBP);
    H.RF_pcirc_RAF = point2fs(fOpt.RF_p_RBP);
    y.xypos = zeros(2, 1);
    R.xypos = repelem(fOpt.sigma2RxyposPVLARA, 2);
    
    % yaw average
    if knob.meas.yawavemode == 1
        y.yawave = zeros(1, 1);
        R.yawave = fOpt.sigma2Ryawave;
    elseif knob.meas.yawavemode == 2
        y.yawave = zeros(2, 1);
        R.yawave = [fOpt.sigma2Ryawave, fOpt.sigma2Rpelvpitchroll];
    elseif (knob.meas.yawavemode == 3) || (knob.meas.yawavemode == 4)
        y.yawave = eye(3, 3);
        R.yawave = [repelem(fOpt.sigma2Rpelvpitchroll, 2), ...
                            fOpt.sigma2Ryawave];
    else
        y.yawave = []; R.yawave = [];
    end
    
    H.E = [eye(3,3) zeros(3,1)];
    % covariance limiter
    H.covlim = [eye(3) zeros(3,3+12) zeros(3, N.state-18);
                zeros(3,6) eye(3) zeros(3,3+6) zeros(3, N.state-18);
                zeros(3,12) eye(3) zeros(3,3) zeros(3, N.state-18) ];
    R.covlim = repelem([fOpt.sigma2RLimPos, fOpt.sigma2RLimPos, ...
                        fOpt.sigma2RLimPos ], 3); 
%     H.covlim = [eye(18) zeros(18, N.state-18)];
%     R.covlim = repelem([fOpt.sigma2RLimPos, fOpt.sigma2RLimOri, ...
%                          fOpt.sigma2RLimPos, fOpt.sigma2RLimOri, ...
%                          fOpt.sigma2RLimPos, fOpt.sigma2RLimOri], 3);
    
    %% Constraint step initialization
    D = struct('i_x', [1 0 0 0]', 'i_y', [0 1 0 0]', 'i_z', [0 0 1 0]', ...
               'icirc_x', point2fs([1 0 0 0]'), 'icirc_y', point2fs([0 1 0 0]'), ...
               'icirc_z', point2fs([0 0 1 0]'), ...
               'i_0', [0 0 0 1]', 'icirc_0', point2fs([0 0 0 1]'), ...
               'I30', [eye(3) zeros(3,1)]);
    D.I30TI30 = D.I30'*D.I30;
    
%     D.ihat_x = hat(D.i_x); D.ihat_y = hat(D.i_y); D.ihat_z = hat(D.i_z); 
    d = {};
    
    % leg length constraint
    D.PV_p_LH = [0 body.PV_d/2 0 1]';   D.PV_pcirc_LH = point2fs(D.PV_p_LH);
    D.PV_p_RH = [0 -body.PV_d/2 0 1]';  D.PV_pcirc_RH = point2fs(D.PV_p_RH);
    D.LF_p_LA = fOpt.LF_p_LA;           D.LF_pcirc_LA = point2fs(D.LF_p_LA);  
    D.RF_p_RA = fOpt.RF_p_RA;           D.RF_pcirc_RA = point2fs(D.RF_p_RA);
    d.lll = 0; d.rll = 0;
    
    % hinge knee and ankle joint
    d.lkah = 0;
    d.rkah = 0;

    % pelvis no pitch and roll
    d.pnpr = 0; % zeros(2,1);
    
    %% Iteration
    for k=2:(N.samples+1)
        kPast = k-1;
        
        %% Prediction update
        if knob.pred.any
            u(:,kPast) = [W_a_.PV(kPast,:)' - g; ...
                     W_a_.LF(kPast,:)' - g; ...
                     W_a_.RF(kPast,:)' - g; ...
                     B_w_.PV(kPast,:)'; ...
                     B_w_.LF(kPast,:)'; ...
                     B_w_.RF(kPast,:)'];
            
            F.AdG = eye(N.state, N.state);
            bigphi = eye(N.state, N.state);
            F.curlyC = F.curlyC0;
            F.vec = F.vec0;
            for i = list.se3State
                sname = char(i);
                bname = sname(end-1:end);
                accName = sprintf('acc%s', bname);
                vecName = sprintf('vel%s', bname);
                accBName = sprintf('accB%s', bname);
                
                % accelerometer bias
                if knob.pred.trackAccBias
                    B_b_a = x.tilde.vec(idx.vec.(accBName),kPast);
                    F.vec(idx.vec.(vecName),idx.vec.(accBName)) = ...
                        -W_R_.(bname)(:,:,kPast)*dt;
                    F.curlyC(idx.(vecName),idx.(accBName)) = -W_R_.(bname)(:,:,kPast)*dt;
                else
                    B_b_a = zeros(3,1);
                end
                
                % convert W_vel to B_vel
                B_R_W = x.tilde.(sname)(1:3,1:3,kPast)';
                % B_R_W = W_R_.(bname)(:,:,kPast)';

                xi = zeros(6,1);
                if knob.pred.posFromVel
                    W_dp = (x.tilde.vec(idx.vec.(vecName),kPast) ...
                            + 0.5*dt*u(idx.u.(accName),kPast));
                    xi(1:3) = B_R_W*W_dp - 0.5*dt*B_b_a;
                
%                     F.curlyC(idx.(sname)(1:3),idx.(sname)(4:6)) = -B_R_W*hat(W_dp);
                end
                if knob.pred.oriFromAVel
                    xi(4:6) = B_w_.(bname)(kPast,:)';
                end
                F.curlyC(idx.(sname)(1:3),idx.(vecName)) = dt*B_R_W;

                bigxi = vec2tran(xi*dt);
                bigphi(idx.(sname),idx.(sname)) = vec2jac(-xi*dt);
                F.AdG(idx.(sname),idx.(sname)) = tranAd(bigxi);

                x.hatPri.(sname)(:,:,k) = x.tilde.(sname)(:,:,kPast)*bigxi;
            end
            F.comb = F.AdG + bigphi*F.curlyC;
            x.hatPri.vec(:,k) = F.vec*x.tilde.vec(:,kPast) + G.vec*u(:,kPast);
            P.hatPri(:,:,k) = F.comb*P.tilde(:,:,kPast)*F.comb' + bigphi*Q.comb*bigphi';
        else
            for i=list.se3State
                x.hatPri.(i)(:,:,k) = x.tilde.(i)(:,:,kPast);
            end
            x.hatPri.vec(:,k) = x.tilde.vec(:,kPast);
            P.hatPri(:,:,k) = P.tilde(:,:,kPast);
        end
        
        %% Measurement update
        N.meas_k = (knob.meas.ori.PV)*3 + ...  % orientation
            (knob.meas.ori.LF)*3 + (knob.meas.ori.RF)*3 + ...
            (knob.meas.zpos.PV && step.PV(kPast)) + ... %zpos PV
            ((knob.meas.zposmode>0) || (knob.meas.zpos.LF && step.LF(kPast))) + ... %zpos LS
            ((knob.meas.zposmode>0) || (knob.meas.zpos.RF && step.RF(kPast))) + ... %zpos RS
            (knob.meas.xyposPVLARA > 0) * 2 + ... % PV = ankle xy pos
            (knob.meas.zupt && step.LF(kPast)) * N.zupt + ... %lzupt
            (knob.meas.zupt && step.RF(kPast)) * N.zupt + ... %rzupt
            (knob.meas.bias && step.LF(kPast)) * N.bias + ... %lbias
            (knob.meas.bias && step.RF(kPast)) * N.bias + ... %rbias
            (knob.meas.yawavemode > 0) * size(y.yawave, 1);
        
        H.comb = zeros(N.meas_k, N.state);
        dy = zeros(N.meas_k, 1);
        R.comb = zeros(N.meas_k, 1);
        N.meas_k = 0;

        % Update orientation from measurement
        for i=list.body
            if knob.meas.ori.(i)
                H.comb(N.meas_k+1:N.meas_k+3, :) = H.ori.(i);
                R.comb(N.meas_k+1:N.meas_k+3, :) = R.ori.(i);
                j = sprintf('W_T_%s', i);
                dy(N.meas_k+1:N.meas_k+3, :) = ...
                    real(rot2vec(x.hatPri.(j)(1:3,1:3,k)' * W_R_.(i)(:,:,kPast)));
                N.meas_k = N.meas_k + 3;
            end
        end
        
        % z pos assumption
        for i=list.se3State
            j = i{1}(end-1:end);
            if (knob.meas.zposmode>0) || (knob.meas.zpos.(j) && step.(j)(kPast)) % step detected
                H.comb(N.meas_k+1, :) = zeros(1, N.state);
                H.comb(N.meas_k+1, idx.(i)) = ...
                    H.i_z' * x.hatPri.(i)(:,:,k) * H.icirc_0;
                R.comb(N.meas_k+1, :) = R.zpos.(j)(kPast);
                dy(N.meas_k+1, :) = y.zpos.(j)(kPast) - H.i_z'*x.hatPri.(i)(:,:,k)*H.i_0;
                N.meas_k = N.meas_k + 1;
            end
        end
        
        % ZUPT
        if knob.meas.zupt && step.LF(kPast)
            H.comb(N.meas_k+1:N.meas_k+N.zupt, :) = H.lzupt;
            R.comb(N.meas_k+1:N.meas_k+N.zupt, :) = R.lzupt;
            dy(N.meas_k+1:N.meas_k+N.zupt, :) = y.lzupt - x.hatPri.vec(idx.vec.velLF,k);
            N.meas_k = N.meas_k + N.zupt;
        end
        if knob.meas.zupt && step.RF(kPast)
            H.comb(N.meas_k+1:N.meas_k+N.zupt, :) = H.rzupt;
            R.comb(N.meas_k+1:N.meas_k+N.zupt, :) = R.rzupt;
            dy(N.meas_k+1:N.meas_k+N.zupt, :) = y.rzupt - x.hatPri.vec(idx.vec.velRF,k);
            N.meas_k = N.meas_k + N.zupt;
        end
        
        if knob.meas.bias && step.LF(kPast)
            H.comb(N.meas_k+1:N.meas_k+N.bias, :) = H.lbias;
            R.comb(N.meas_k+1:N.meas_k+N.bias, :) = R.lbias;
            dy(N.meas_k+1:N.meas_k+N.bias, :) = ...
                W_R_.LF(:,:,kPast)'*u(idx.u.accLF,kPast) ...
                - x.hatPri.vec(idx.vec.accBLF,k);
            N.meas_k = N.meas_k + N.bias;
        end
        
        if knob.meas.bias && step.RF(kPast)
            H.comb(N.meas_k+1:N.meas_k+N.bias, :) = H.rbias;
            R.comb(N.meas_k+1:N.meas_k+N.bias, :) = R.rbias;
            dy(N.meas_k+1:N.meas_k+N.bias, :) = ...
                W_R_.RF(:,:,kPast)'*u(idx.u.accRF,kPast) ...
                - x.hatPri.vec(idx.vec.accBRF,k);
            N.meas_k = N.meas_k + N.bias;
        end
        
        % PV = ankle XY position assumption
        if knob.meas.xyposPVLARA
            H.comb(N.meas_k+1:N.meas_k+2,:) = zeros(2,N.state);
            H.comb(N.meas_k+1:N.meas_k+2,idx.W_T_PV) = - H.xyposE * x.hatPri.W_T_PV(:,:,k) ...
                                     * H.icirc_0;
            H.comb(N.meas_k+1:N.meas_k+2,idx.W_T_LF) = H.xyposE * x.hatPri.W_T_LF(:,:,k) ...
                                     * H.LF_pcirc_LAF / 2;
            H.comb(N.meas_k+1:N.meas_k+2,idx.W_T_RF) = H.xyposE * x.hatPri.W_T_RF(:,:,k) ...
                                     * H.RF_pcirc_RAF / 2;
            dy(N.meas_k+1:N.meas_k+2,:) = y.xypos - H.xyposE * ...
                ( x.hatPri.W_T_LF(:,:,k) * fOpt.RF_p_RBP / 2 ...
                  + x.hatPri.W_T_RF(:,:,k) * fOpt.LF_p_LBP / 2 ...
                  - x.hatPri.W_T_PV(:,:,k) * H.i_0); 
            R.comb(N.meas_k+1:N.meas_k+2,:) = R.xypos;
            N.meas_k = N.meas_k + 2;
        end
        
        % PV yaw = Foot yaw average
        if (knob.meas.yawavemode == 1) || (knob.meas.yawavemode == 2)
            H.comb(N.meas_k+1,:) = zeros(1,N.state);
            H.comb(N.meas_k+1,idx.W_T_PV) = calcYawJac(x.hatPri.W_T_PV(:,:,k));
            H.comb(N.meas_k+1,idx.W_T_LF) = -calcYawJacFT(x.hatPri.W_T_LF(:,:,k))/2.0;
            H.comb(N.meas_k+1,idx.W_T_RF) = -calcYawJacFT(x.hatPri.W_T_RF(:,:,k))/2.0;
            dy(N.meas_k+1,:) = y.yawave(1) - ...
                ( calcYaw(x.hatPri.W_T_PV(:,:,k)) ...
                  - calcYawFT(x.hatPri.W_T_LF(:,:,k))/2.0 ...
                  - calcYawFT(x.hatPri.W_T_RF(:,:,k))/2.0 ); 
            
            if knob.meas.yawavemode == 1
                R.comb(N.meas_k+1:N.meas_k+1,:) = R.yawave;
                N.meas_k = N.meas_k + 1;
            else
                H.comb(N.meas_k+2,:) = zeros(1,N.state);
                H.comb(N.meas_k+2,idx.W_T_PV) = calcPitchJac(x.hatPri.W_T_PV(:,:,k));
            
                dy(N.meas_k+2,:) = y.yawave(2) - calcPitch(x.hatPri.W_T_PV(:,:,k));
                
                R.comb(N.meas_k+1:N.meas_k+2,:) = R.yawave;
                N.meas_k = N.meas_k + 2;
            end
        elseif (knob.meas.yawavemode == 3) || (knob.meas.yawavemode == 4)
            H.comb(N.meas_k+1:N.meas_k+3,:) = H.ori.PV; % identity of PV ori
            W_r_PVT_z = H.E * x.hatPri.W_T_PV(:,:,k)' * D.i_z;
            
            if knob.meas.yawavemode == 3
                H.comb(N.meas_k+1:N.meas_k+3,idx.W_T_LF) = - W_r_PVT_z * calcYawJacFT(x.hatPri.W_T_LF(:,:,k))/2.0;
                H.comb(N.meas_k+1:N.meas_k+3,idx.W_T_RF) = - W_r_PVT_z * calcYawJacFT(x.hatPri.W_T_RF(:,:,k))/2.0;
                
                % approach 1
                LFYaw = calcYawFT(x.hatPri.W_T_LF(:,:,k));
%                 RFYaw = calcYawFT2(x.hatPri.W_T_RF(:,:,k), LFYaw);
                RFYaw = calcYawFT(x.hatPri.W_T_RF(:,:,k));
                W_Ryaw_PV = vec2rot( ( LFYaw/2.0 + RFYaw/2.0 ) *D.i_z(1:3) );
                % approach 2
                % W_Ryaw_PV = vec2rot( calcYawFT( vec2rot( 0.5*rot2vec(x.hatPri.W_T_RF(1:3,1:3,k) * x.hatPri.W_T_LF(1:3,1:3,k)') ) * x.hatPri.W_T_LF(1:3,1:3,k) ) * D.i_z(1:3) );
            elseif knob.meas.yawavemode == 4
                H.comb(N.meas_k+1:N.meas_k+3,idx.W_T_LF) = - W_r_PVT_z * ...
                    calcYawJacFTAve(x.hatPri.W_T_LF(:,:,k), x.hatPri.W_T_RF(:,:,k));
                H.comb(N.meas_k+1:N.meas_k+3,idx.W_T_RF) = - W_r_PVT_z * ...
                    calcYawJacFTAve(x.hatPri.W_T_RF(:,:,k), x.hatPri.W_T_LF(:,:,k));
                W_Ryaw_PV = vec2rot( calcYawFT( x.hatPri.W_T_LF(:,:,k) + x.hatPri.W_T_RF(:,:,k) ) *D.i_z(1:3) );
            end

            % approach 1
            % if knob.meas.yawavezeropelvpitchroll2
            %    LFYaw = calcYawFT(x.hatPri.W_T_LF(:,:,k));
            %    RFYaw = calcYawFT2(x.hatPri.W_T_RF(:,:,k), LFYaw);
            %    W_Ryaw_PV = vec2rot( ( LFYaw/2.0 + RFYaw/2.0 ) *D.i_z(1:3) );
            %else
            
            %end
            
            R.comb(N.meas_k+1:N.meas_k+3,:) = R.yawave;
            dy(N.meas_k+1:N.meas_k+3, :) = ...
                real(rot2vec( x.hatPri.W_T_PV(1:3,1:3,k)' * W_Ryaw_PV * y.yawave));
            
            N.meas_k = N.meas_k + 3;
        end
        
        if N.meas_k > 0 || knob.meas.covlim
            K = P.hatPri(:,:,k)*H.comb'/(H.comb*P.hatPri(:,:,k)*H.comb' + diag(R.comb));
            bigphi = eye(N.state, N.state);
            measUpt = K*(dy);
            debug_dat.measUptPos(:,kPast) = measUpt;
            
            for i=list.se3State
                x.hatPos.(i)(:,:,k) = x.hatPri.(i)(:,:,k)*vec2tran(measUpt(idx.(i)));
                bigphi(idx.(i), idx.(i)) = vec2jac(-measUpt(idx.(i)));
            end
            x.hatPos.vec(:,k) = x.hatPri.vec(:,k) + measUpt(idx.vecState);
            
            if knob.meas.covlim
                H.comb2 = [H.comb; H.covlim];
                R.comb2 = [R.comb' R.covlim];
                K2 = P.hatPri(:,:,k)*H.comb2'/(H.comb2*P.hatPri(:,:,k)*H.comb2' + diag(R.comb2));
                P.hatPos(:,:,k) = bigphi*(I_N-K2*H.comb2)*P.hatPri(:,:,k)*bigphi';
            else
                P.hatPos(:,:,k) = bigphi*(I_N-K*H.comb)*P.hatPri(:,:,k)*bigphi';
            end
        else
            for j=list.se3State
                x.hatPos.(j)(:,:,k) = x.hatPri.(j)(:,:,k);
            end
            x.hatPos.vec(:,k) = x.hatPri.vec(:,k);
            P.hatPos(:,:,k) = P.hatPri(:,:,k);
        end
    
        %% Constraint update
        if knob.cstr.any
            n_LL = D.I30*(x.hatPos.W_T_PV(:,:,k)*D.PV_p_LH - ...
                                x.hatPos.W_T_LF(:,:,k)*D.LF_p_LA);
            n_RL = D.I30*(x.hatPos.W_T_PV(:,:,k)*D.PV_p_RH - ...
                                x.hatPos.W_T_RF(:,:,k)*D.RF_p_RA);

            N.cstr_k = (knob.cstr.maxleglength  && norm(n_LL) > body.LL_d ) + ...
                    (knob.cstr.maxleglength  && norm(n_RL) > body.RL_d ) + ...
                    (knob.cstr.hingekneeankle > 0) * 2 + ...
                    (knob.cstr.pelvnopitchroll > 0) * 1;
            
            D.comb = zeros(N.cstr_k, N.state);
            dy = zeros(N.cstr_k, 1);
            N.cstr_k = 0;
            
            % leg length constraint (left)
            if (knob.cstr.maxleglength  && norm(n_LL) > body.LL_d )
                D.comb(N.cstr_k+1,idx.W_T_PV) = 2*n_LL'*D.I30*...
                                x.hatPos.W_T_PV(:,:,k)*D.PV_pcirc_LH;
                D.comb(N.cstr_k+1,idx.W_T_LF) = -2*n_LL'*D.I30*...
                                x.hatPos.W_T_LF(:,:,k)*D.LF_pcirc_LA;
                dy(N.cstr_k+1,:) = d.lll - n_LL'*n_LL + body.LL_d.^2;
                N.cstr_k = N.cstr_k+1;
            end
            
            % leg length constraint (right)
            if (knob.cstr.maxleglength  && norm(n_RL) > body.RL_d )
                D.comb(N.cstr_k+1,idx.W_T_PV) = 2*n_RL'*D.I30*...
                                x.hatPos.W_T_PV(:,:,k)*D.PV_pcirc_RH;
                D.comb(N.cstr_k+1,idx.W_T_RF) = -2*n_RL'*D.I30*...
                                x.hatPos.W_T_RF(:,:,k)*D.RF_pcirc_RA;
                dy(N.cstr_k+1,:) = d.rll - n_RL'*n_RL + body.RL_d.^2;
                N.cstr_k = N.cstr_k+1;
            end
            
            % hinge knee and ankle joint
            if knob.cstr.hingekneeankle
                W_r_LF_y = x.hatPos.W_T_LF(:,:,k) * D.i_y;
                W_r_RF_y = x.hatPos.W_T_RF(:,:,k) * D.i_y;
                D.comb(N.cstr_k+1,idx.W_T_PV) = W_r_LF_y' * D.I30TI30 ...
                                  * x.hatPos.W_T_PV(:,:,k) * D.PV_pcirc_LH;
                D.comb(N.cstr_k+1,idx.W_T_LF) = -W_r_LF_y' * D.I30TI30 ...
                                  * x.hatPos.W_T_LF(:,:,k)*D.LF_pcirc_LA ...
                                  + n_LL'*D.I30*x.hatPos.W_T_LF(:,:,k)*D.icirc_y;
                D.comb(N.cstr_k+2,idx.W_T_PV) = W_r_RF_y'*D.I30TI30 ...
                                  * x.hatPos.W_T_PV(:,:,k) * D.PV_pcirc_RH;
                D.comb(N.cstr_k+2,idx.W_T_RF) = -W_r_RF_y'*D.I30TI30 ...
                                  * x.hatPos.W_T_RF(:,:,k)*D.RF_pcirc_RA ...
                                  + n_RL'*D.I30*x.hatPos.W_T_RF(:,:,k)*D.icirc_y;
                dy(N.cstr_k+(1:2),:) = [d.lkah; d.rkah] - ...
                    [(W_r_LF_y)'*D.I30'*n_LL; (W_r_RF_y)'*D.I30'*n_RL];
                N.cstr_k = N.cstr_k + 2;
            end
            
            % pelvis no pitch and roll
            if knob.cstr.pelvnopitchroll
                W_r_PV_x = x.hatPos.W_T_PV(:,:,k) * D.i_x;
%                 W_r_PV_y = x.hatPos.W_T_PV(:,:,k) * D.i_y;
%                 W_r_PV_z = x.hatPos.W_T_PV(:,:,k) * D.i_z;
                D.comb(N.cstr_k+1,idx.W_T_PV) = D.i_z' * x.hatPos.W_T_PV(:,:,k) * D.icirc_x;
%                 D.comb(N.cstr_k+2,idx.W_T_PV) = D.i_z' * x.hatPos.W_T_PV(:,:,k) * D.icirc_y;               
%                 D.comb(N.cstr_k+2,idx.W_T_PV) = D.i_z' * x.hatPos.W_T_PV(:,:,k) * D.icirc_z;
                dy(N.cstr_k+1,:) = d.pnpr - D.i_z'*W_r_PV_x;
%                 dy(N.cstr_k+(1:2),:) = d.pnpr - [D.i_z'*W_r_PV_x; D.i_z'*W_r_PV_y];
                N.cstr_k = N.cstr_k + 1;
            end
        
            K = P.hatPos(:,:,k)*D.comb'/(D.comb*P.hatPos(:,:,k)*D.comb');
            if knob.cstr.Pupdate
                P.tilde(:,:,k) = (I_N-K*D.comb)*P.hatPos(:,:,k);
            else
                P.tilde(:,:,k) = P.hatPos(:,:,k);
            end
            
            measUpt = K*dy;
            debug_dat.measUptTilde(:,kPast) = measUpt;
                        
            for i=list.se3State
                x.tilde.(i)(:,:,k) = x.hatPos.(i)(:,:,k)*...
                                        vec2tran(measUpt(idx.(i)));
            end
            x.tilde.vec(:,k) = x.hatPos.vec(:,k) + measUpt(idx.vecState);
        else
            for i=list.se3State
                x.tilde.(i)(:,:,k) = x.hatPos.(i)(:,:,k);
            end
            x.tilde.vec(:,k) = x.hatPos.vec(:,k);
            P.tilde(:,:,k) = P.hatPos(:,:,k);
        end
    end
    
    %% remove offset state (k=0)
    for i=list.kfStep
        for j=list.se3State
            x.(i).(j) = x.(i).(j)(:,:,2:end);
        end
        x.(i).vec = x.(i).vec(:,2:end);
        P.(i) = P.(i)(:,:,2:end);
    end
    
    %% Postprocessing + Compute debug information
    post = struct();
    post.LFTy = squeeze(x.tilde.W_T_LF(1:3,2,:));
    post.LFTz = squeeze(x.tilde.W_T_LF(1:3,3,:))';
    post.RFTy = squeeze(x.tilde.W_T_RF(1:3,2,:));
    post.RFTz = squeeze(x.tilde.W_T_RF(1:3,3,:))';
    post.PELVy = squeeze(x.tilde.W_T_PV(1:3,2,:))';
    post.MIDPEL = squeeze(x.tilde.W_T_PV(1:3,4,:))';
    
    debug_dat.LTIO = zeros(4,N.samples);
    debug_dat.RTIO = zeros(4,N.samples);
    for k=1:N.samples
        debug_dat.LTIO(:,k) = x.tilde.W_T_LF(:,:,k)*fOpt.LF_p_LA;
        debug_dat.RTIO(:,k) = x.tilde.W_T_RF(:,:,k)*fOpt.RF_p_RA;
    end
    debug_dat.LTIO = debug_dat.LTIO(1:3,:)';
    debug_dat.RTIO = debug_dat.RTIO(1:3,:)';
    
    debug_dat.LTOE = debug_dat.LTIO + body.LF_d * post.LFTz;
    debug_dat.RTOE = debug_dat.RTIO + body.RF_d * post.RFTz;
    debug_dat.LFEP = post.MIDPEL + body.PV_d/2 * post.PELVy;
    debug_dat.RFEP = post.MIDPEL - body.PV_d/2 * post.PELVy;
    
    % calculate LFEO and RFEO
    post.v_LFEPLTIO = debug_dat.LFEP - debug_dat.LTIO;
    post.v_RFEPRTIO = debug_dat.RFEP - debug_dat.RTIO;
    post.d0_LFEPLTIO = vecnorm(post.v_LFEPLTIO, 2, 2);
    post.d0_RFEPRTIO = vecnorm(post.v_RFEPRTIO, 2, 2);
    post.d_LFEPLTIO = min(post.d0_LFEPLTIO, body.LT_d+body.LS_d);
    post.d_RFEPRTIO = min(post.d0_RFEPRTIO, body.RT_d+body.RS_d);
    % without real(), theta_LKNE will have a very small complex component
    post.theta_LKNE = real(acos((body.LT_d.^2 + body.LS_d.^2 - post.d_LFEPLTIO.^2)/(2*body.LT_d*body.LS_d)));
    post.theta_RKNE = real(acos((body.RT_d.^2 + body.RS_d.^2 - post.d_RFEPRTIO.^2)/(2*body.RT_d*body.RS_d)));
    post.theta_LSHK = asin(body.LT_d.*sin(post.theta_LKNE)./post.d_LFEPLTIO);
    post.theta_RSHK = asin(body.RT_d.*sin(post.theta_RKNE)./post.d_RFEPRTIO);
    
    % Matlab quat library version
    post.qLFEPLFEO = axang2quat([post.LFTy' post.theta_LSHK]);
    post.qRFEPRFEO = axang2quat([post.RFTy' post.theta_RSHK]);
    debug_dat.LFEO = debug_dat.LTIO + quatrotate(quatconj(post.qLFEPLFEO), ...
                            post.v_LFEPLTIO./post.d0_LFEPLTIO.*body.LS_d);
    debug_dat.RFEO = debug_dat.RTIO + quatrotate(quatconj(post.qRFEPRFEO), ...
                            post.v_RFEPRTIO./post.d0_RFEPRTIO.*body.RS_d);
    % Lie group version (the above and this one is equivalent. chose the
    % above version as built in MATLAB library should be faster
%     post.RLFEPLFEO = zeros(3,3,N.samples);
%     post.RRFEPRFEO = zeros(3,3,N.samples);
%     debug_dat.LFEO = zeros(N.samples,3);
%     debug_dat.RFEO = zeros(N.samples,3);
%     for i=1:N.samples
%         post.RLFEPLFEO(:,:,i) = vec2rot(post.LFTy(:,i)*post.theta_LSHK(i));
%         post.RRFEPRFEO(:,:,i) = vec2rot(post.RFTy(:,i)*post.theta_RSHK(i));
%         debug_dat.LFEO(i,:) = debug_dat.LTIO(i,:) + ...
%             (post.RLFEPLFEO(:,:,i)*(post.v_LFEPLTIO(i,:)'./post.d0_LFEPLTIO(i).*body.LS_d))';
%         debug_dat.RFEO(i,:) = debug_dat.RTIO(i,:) + ...
%             (post.RRFEPRFEO(:,:,i)*(post.v_RFEPRTIO(i,:)'./post.d0_RFEPRTIO(i).*body.RS_d))';  
%     end                       
    
    post.W_R_LT = zeros(3,3,N.samples);
    post.W_R_RT = zeros(3,3,N.samples);
    post.LFEM_z = (debug_dat.LFEP-debug_dat.LFEO)';
    post.LFEM_z = post.LFEM_z ./ vecnorm(post.LFEM_z, 2, 1);
    post.W_R_LT(:,3,:) = post.LFEM_z;
    post.W_R_LT(:,2,:) = post.LFTy;
    post.W_R_LT(:,1,:) = cross(post.LFTy, post.LFEM_z); 
    post.RFEM_z = (debug_dat.RFEP-debug_dat.RFEO)';
    post.RFEM_z = post.RFEM_z ./ vecnorm(post.RFEM_z, 2, 1);
    post.W_R_RT(:,3,:) = post.RFEM_z;
    post.W_R_RT(:,2,:) = post.RFTy;
    post.W_R_RT(:,1,:) = cross(post.RFTy, post.RFEM_z); 
    debug_dat.qLTH = rotm2quat(post.W_R_LT);
    debug_dat.qRTH = rotm2quat(post.W_R_RT);
    
    post.W_R_LS = zeros(3,3,N.samples);
    post.W_R_RS = zeros(3,3,N.samples);  
    
    post.LSHK_z = (debug_dat.LFEO-debug_dat.LTIO)';
    post.LSHK_z = post.LSHK_z ./ vecnorm(post.LSHK_z, 2, 1);
    post.W_R_LS(:,3,:) = post.LSHK_z;
    post.W_R_LS(:,2,:) = post.LFTy;
    post.W_R_LS(:,1,:) = cross(post.LFTy, post.LSHK_z); 
    post.RSHK_z = (debug_dat.RFEO-debug_dat.RTIO)';
    post.RSHK_z = post.RSHK_z ./ vecnorm(post.RSHK_z, 2, 1);
    post.W_R_RS(:,3,:) = post.RSHK_z;
    post.W_R_RS(:,2,:) = post.RFTy;
    post.W_R_RS(:,1,:) = cross(post.RFTy, post.RSHK_z);
    debug_dat.qLSK = rotm2quat(post.W_R_LS);
    debug_dat.qRSK = rotm2quat(post.W_R_RS);
    
    if fOpt.fulldebug
%     debug_dat.u = u;
        debug_dat.x.hatPri = x.hatPri;
        debug_dat.x.hatPos = x.hatPos;
        debug_dat.x.tilde = x.tilde;
        debug_dat.P.hatPri = P.hatPri;
        debug_dat.P.hatPos = P.hatPos;
    end
    debug_dat.P.tilde = P.tilde;
    debug_dat.fs = fOpt.fs;
    debug_dat.body = body;
    
    xtilde = x.tilde;
end

function out = calcYaw(W_T_B)
    out = atan2(W_T_B(2,1), W_T_B(1,1));
end

function out = calcPitch(W_T_B)
    out = atan2(W_T_B(3,1), W_T_B(1,1));
end

function out = calcYawFT(W_T_B)
    out = atan2(W_T_B(2,3), W_T_B(1,3));
end

function minOut = calcYawFT2(W_T_B, refAngle)
    if nargin <= 1, refAngle = 0; end
    
    out = atan2(W_T_B(2,3), W_T_B(1,3));
    
    minOut = out; minDiff = abs(out-refAngle);
    if(abs(out+2*pi-refAngle)<minDiff)
        minOut = out + 2*pi;
        minDiff = abs(out+2*pi-refAngle);
    end
    if(abs(out-2*pi-refAngle)<minDiff)
        minOut = out - 2*pi;
        minDiff = abs(out-2*pi-refAngle);
    end
end

function out = calcYawJac(W_T_B)
    D = struct('i_x', [1 0 0 0]', 'i_y', [0 1 0 0]', ...
               'icirc_x', point2fs([1 0 0 0]') );
    out = (W_T_B(1,1).*(D.i_y'*W_T_B*D.icirc_x) - W_T_B(2,1).*(D.i_x'*W_T_B*D.icirc_x) ) ...
            ./ (W_T_B(2,1)^2 + W_T_B(1,1)^2);
end

function out = calcPitchJac(W_T_B)
    D = struct('i_x', [1 0 0 0]', 'i_z', [0 0 1 0]', ...
               'icirc_x', point2fs([1 0 0 0]') );
    out = (W_T_B(1,1).*(D.i_z'*W_T_B*D.icirc_x) - W_T_B(3,1).*(D.i_x'*W_T_B*D.icirc_x) ) ...
            ./ (W_T_B(3,1)^2 + W_T_B(1,1)^2);
end

function out = calcYawJacFT(W_T_B)
    D = struct('i_x', [1 0 0 0]', 'i_y', [0 1 0 0]', 'i_z', [0 0 1 0]', ...
               'icirc_z', point2fs([0 0 1 0]') );
    out = (W_T_B(1,3).*(D.i_y'*W_T_B*D.icirc_z) - W_T_B(2,3).*(D.i_x'*W_T_B*D.icirc_z) ) ...
            ./ (W_T_B(2,3)^2 + W_T_B(1,3)^2);
end

function out = calcYawJacFTAve(W_T_A, W_T_B)
    D = struct('i_x', [1 0 0 0]', 'i_y', [0 1 0 0]', 'i_z', [0 0 1 0]', ...
               'icirc_z', point2fs([0 0 1 0]') );
    tx = W_T_A(1,3)+W_T_B(1,3);
    ty = W_T_A(2,3)+W_T_B(2,3);
    out = (tx.*(D.i_y'*W_T_A*D.icirc_z) - ty.*(D.i_x'*W_T_A*D.icirc_z) ) ...
            ./ (tx^2 + ty^2);
end