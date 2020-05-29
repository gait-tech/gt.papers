function [ xtilde, debug_dat ] = lieekf_3_kmus(x0, P0, ...
    B_a_, step, W_R_, B_w_, body, options)
    % @brief Constrained Lie group based extended Kalman filter implementation
    % @author Luke Wicent Sy
    % @date 14 May 2019
    %
    % EKF using Lie group/algebra representation
    % 3 KMUs presumably worn on the body in the following configuration: 
    % mid pelvis, left ankle, right ankle
    % 
    % Coding notation based on <a href="http://paulfurgale.info/news/2014/6/9/representing-robot-pose-the-good-the-bad-and-the-ugly">link</a>
    %
    % More detail about options
    %      fs: sampling frequency of the magnetic and inertial measurement units
    %      applyPred: 3 digit ZYX
    %          X:   1st bit use velocity to predict position
    %      applyMeas: 3 digit ZYX
    %          X:   1st bit Zupt and Ankle zpos
    %               2nd bit Update angular velocity
    %               3rd bit Orientation update
    %          Y:   1st bit Pelvis xy=ankle average
    %               2nd bit Pelvis z=initial height
    %          Z:   1st bit covariance limiter
    %      applyCstr: 3 digit ZYX
    %          X:   1st bit enforce thigh length
    %               2nd bit enforce hinge knee joint
    %               3rd bit enforce knee range of motion
    %          Z:   1st bit do P update
    %
    % :param x0 initial state in the GFR
    % :param P0 initial covariance
    % :param bodyAcc acceleration of PV, LS, RS in the body frame
    % :param step boolean vector of PV, LS, RS indicating step detection
    % :param qOri PV, LS, RS orientation in the GFR (rotm)
    % :param wbody PV, LS, RS angular velocity in the body frame
    % :param body Length of PV_d (pelvis), RT_d and LT_d (r/l femur), RS_d and LS_d (r/l tibia)
    % :param uwb_mea a structure containing the range measurements (m) between
    % :param options struct containing the estimator settings:

    N = {};
    [N.samples, ~] = size(B_a_.PV);
    
    %% input parsing
    fOpt = struct('fs', 100, ...
          'applyPred', 1, 'applyMeas', 1, 'applyCstr', 1, ...
          'sigma2QAccPV', 1, 'sigma2QAccLS', 1, 'sigma2QAccRS', 1, ...
          'sigma2QGyrPV', 1, 'sigma2QGyrLS', 1, 'sigma2QGyrRS', 1, ...
          'sigma2QPosMP', 1, 'sigma2QPosLA', 1, 'sigma2QPosRA', 1, ...
          'sigma2QOriPV', 1e3, 'sigma2QOriLS', 1e3, 'sigma2QOriRS', 1e3, ...
          'sigma2QVelMP', 1, 'sigma2QVelLA', 1, 'sigma2QVelRA', 1, ...
          'sigma2ROriPV', 1e-2, 'sigma2ROriLS', 1e-2, 'sigma2ROriRS', 1e-2, ...
          'sigma2RZPosPV', 1e-1, 'sigma2RZPosLS', 1e-4, 'sigma2RZPosRS', 1e-4, ...
          'sigma2RZuptLA', 1e-2, 'sigma2RZuptRA', 1e-2, ...
          'sigma2RXYPosPVLSRS', 1e2, ...
          'sigma2RLimPos', 1e1, 'sigma2RLimOri', 1e1, ...
          'alphaLKmin', 0, 'alphaLKmax', pi*8/9, ...
          'alphaRKmin', 0, 'alphaRKmax', pi*8/9 );
    optionFieldNames = fieldnames(options);
    for i=1:length(optionFieldNames)
        if ~isfield(fOpt, optionFieldNames{i})
            error("Field %s is not a valid option", optionFieldNames{i});
        end
        fOpt.(optionFieldNames{i}) = options.(optionFieldNames{i});
    end
    
    %% Feature on off
    knob = struct();
    knob.apModTen = mod(fOpt.applyPred, 10);
    knob.apTenDig = mod(idivide(int32(fOpt.applyPred), 10, 'floor'), 10);
    knob.apHunDig = mod(idivide(int32(fOpt.applyPred), 100, 'floor'), 10);
    knob.amModTen = mod(fOpt.applyMeas, 10);
    knob.amTenDig = mod(idivide(int32(fOpt.applyMeas), 10, 'floor'), 10);
    knob.amHunDig = mod(idivide(int32(fOpt.applyMeas), 100, 'floor'), 10);
    knob.acModTen = mod(fOpt.applyCstr, 10);
    knob.acTenDig = mod(idivide(int32(fOpt.applyCstr), 10, 'floor'), 10);
    knob.acHunDig = mod(idivide(int32(fOpt.applyCstr), 100, 'floor'), 10);
    
    % Prediction knobs
    % zero velocity otherwise
    knob.pred.posfromvel = bitand(knob.apModTen, 1); 
    
    % Measurement knobs
    knob.meas = struct('ori', bitand(knob.amModTen, 4), ...
        'zpos', struct('LS', false, 'RS', false, 'PV', false), ...
        'zupt', false, 'xyposPVLSRS', false );   
%          X: Ori always on
%               1st bit Zupt and Ankle zpos
%               2nd bit Update angular velocity
    if bitand(knob.amModTen, 1)
        knob.meas.zupt = true;
        knob.meas.zpos.LS = true;
        knob.meas.zpos.RS = true;
    else
        step.LS = false(N.samples, 1);
        step.RS = false(N.samples, 1);
    end
    knob.meas.covlim = knob.amHunDig;
%          Y:   1st bit Pelvis assumption (xy=ankle average, z=initial height)
%               2nd bit Vel cstr Y and Z

    knob.meas.xyposPVLSRS = bitand(knob.amTenDig, 1);
    if bitand(knob.amTenDig, 2)
        knob.meas.zpos.PV = true;
        step.PV = true(N.samples, 1);
    else
        knob.meas.zpos.PV = false;
        step.PV = false(N.samples, 1);
    end
    
    % Constraint
    knob.cstr.thighlength = bitand(knob.acModTen, 1);
    knob.cstr.hingeknee = bitand(knob.acModTen, 2);
    knob.cstr.kneerom = bitand(knob.acModTen, 4);
    knob.cstr.on = knob.cstr.thighlength | knob.cstr.hingeknee | knob.cstr.kneerom;
    knob.cstr.Pupdate = bitand(knob.acHunDig, 1);
        
    % wOri = struct('RPV', wMP, 'LSK', wLA, 'RSK', wRA);
    g = [0 0 9.81]';        % gravity
    dt = 1/(fOpt.fs);       % assume constant sampling interval
    dt2 = 0.5*dt^2;
    
    %% State and error covariance initialization
    bodyList = ["PV", "LS", "RS"];
    se3StateList = {'W_T_PV', 'W_T_LS', 'W_T_RS'};
    N.se3StateList = length(se3StateList);
    N.se3State = N.se3StateList*6;
    N.r3State = 9;
    
    N.state = N.se3State + N.r3State;
    
    idx = struct('se3state', 1:18, ...
                 'W_T_PV', 1:6, 'W_T_LS', 7:12, 'W_T_RS', 13:18, ...
                 'vec', (N.se3State+1):N.state, ...
                 'vec0', 1:N.r3State, ...
                 'avelState', (N.se3State+10):N.state);
    r3StateList = {'vecVelPV', 'vecVelLS', 'vecVelRS'};
    r3StateList2 = {'velPV', 'velLS', 'velRS'};
    for i=1:length(r3StateList)
        idx.(r3StateList{i}) = (1:3) + ((i-1)*3);
        idx.(r3StateList2{i}) = (1:3) + ((i-1)*3) + N.se3State;
    end
    I_N = eye(N.state);
    
    x.hatPri = struct();     % prediction update state
    x.hatPos = struct();     % measurement update state
    x.tilde = struct();      % constraint update state
    for i=1:N.se3StateList
        sname = se3StateList{i};
        x.hatPri.(sname) = nan(4,4,N.samples+1);
        x.hatPos.(sname) = nan(4,4,N.samples+1);
        x.tilde.(sname) = nan(4,4,N.samples+1);
    end
    x.hatPri.vec = nan(N.r3State,N.samples+1);
    x.hatPos.vec = nan(N.r3State,N.samples+1);
    x.tilde.vec = nan(N.r3State,N.samples+1);
    
    P.hatPri = nan(N.state,N.state,N.samples+1); % prediction update covariance
    P.hatPos = nan(N.state,N.state,N.samples+1); % measurement update covariance
    P.tilde = nan(N.state,N.state,N.samples+1); % constraint update covariance
    
    if isscalar(P0)
        P0 = I_N*P0;
    end
    P.tilde(:,:,1) = P0;
    
    F = struct(); G = struct(); Q = struct();
    F.vec = eye(9,9);
    G.vec = [dt*eye(9,9) zeros(9,9)];
    F.curlyC = zeros(N.state, N.state);
    for i=1:N.se3StateList
        sname = se3StateList{i};
        vname = sprintf('vel%s', sname(end-1:end));
        
        F.curlyC(idx.(sname)(1:3),idx.(vname)) = dt*eye(3,3);
    end
    G.naive = zeros(N.state, 18);
    G.naive(idx.W_T_PV(1:3), 1:3) = dt2*eye(3,3);
    G.naive(idx.W_T_LS(1:3), 4:6) = dt2*eye(3,3);
    G.naive(idx.W_T_RS(1:3), 7:9) = dt2*eye(3,3);
    G.naive(idx.vec, 1:18) = G.vec;
    Q.comb = diag(repelem([fOpt.sigma2QAccPV*dt2, fOpt.sigma2QOriPV, ...
                      fOpt.sigma2QAccLS*dt2, fOpt.sigma2QOriLS, ...
                      fOpt.sigma2QAccRS*dt2, fOpt.sigma2QOriRS, ...
                      fOpt.sigma2QVelMP*dt, fOpt.sigma2QVelLA*dt, ...
                      fOpt.sigma2QVelRA*dt], 3));
    
    %% Prepare input
    for i=bodyList
        if ndims(W_R_.(i)) == 2 && size(W_R_.(i), 2) == 4 % quaternion rep
            W_R_.(i) = quat2rotm(W_R_.(i));
        end
    end
    
    %% Set k=1 states (initial state taken from input)
    % SE(3) state initilization (i.e., position and orientation)
    x.tilde.W_T_PV(:,:,1) = x0.W_T_PV;
    x.tilde.W_T_LS(:,:,1) = x0.W_T_LS;
    x.tilde.W_T_RS(:,:,1) = x0.W_T_RS;
    x.tilde.vec(:,1) = x0.vec;
    
    %% H matrix initialization
    H0 = {}; y0 = {}; R0 = {};
    
    % Orientation check
    H0.ori = {};    R0.ori = {};
    for i=1:N.se3StateList
        sname = se3StateList{i};
        H0.ori.(sname) = zeros(3, N.state);
        H0.ori.(sname)(:, idx.(sname)) = [zeros(3,3) eye(3,3)];
        fOptparam = sprintf('sigma2ROri%s', sname(end-1:end));
        R0.ori.(sname) = repelem(fOpt.(fOptparam), 3);
    end
    
    % ZUPT
    % zero velocity and angular velocity
%     idx.lzupt = [idx.velLS idx.aVelLS];
%     idx.rzupt = [idx.velRS idx.aVelRS];
%     idx.vecLZupt = [idx.vecVelLS idx.vecAVelLS];
%     idx.vecRZupt = [idx.vecVelRS idx.vecAVelRS];
    idx.lzupt = idx.velLS;
    idx.rzupt = idx.velRS;
    idx.vecLZupt = idx.vecVelLS;
    idx.vecRZupt = idx.vecVelRS;

    N.zupt = length(idx.lzupt);
    H0.lzupt = zeros(N.zupt, N.state);
    H0.lzupt(:,idx.lzupt) = eye(N.zupt, N.zupt);
    y0.lzupt = zeros(N.zupt, 1);
    R0.lzupt = repelem(fOpt.sigma2RZuptLA, N.zupt);
    H0.rzupt = zeros(N.zupt, N.state);
    H0.rzupt(:,idx.rzupt) = eye(N.zupt, N.zupt);
    y0.rzupt = zeros(N.zupt, 1);
    R0.rzupt = repelem(fOpt.sigma2RZuptRA, N.zupt);
        
    % zpos = some value assumption (e.g., flat floor)
    % H.zposMP, H.zposLS, H.zposRS varies with t=k
    H0.zposAnkPCircdot = point2fs([0 0 0 1]');
    y0.zpos = {}; R0.zpos = {};
    for i=1:N.se3StateList
        sname = se3StateList{i};
        fOptparam = sprintf('sigma2RZPos%s', sname(end-1:end));
        R0.zpos.(sname) = fOpt.(fOptparam);
    end
    y0.zpos.W_T_PV = x.tilde.W_T_PV(3,4,1);
    y0.zpos.W_T_LS =  min([x.tilde.W_T_LS(3,4,1), x.tilde.W_T_RS(3,4,1)]);
    y0.zpos.W_T_RS =  y0.zpos.W_T_LS;
    
    % pelvis = ankle x y pos
    H0.p0 = [0 0 0 1]';
    H0.p0Circdot = point2fs([0 0 0 1]');
    H0.xyposDT = [eye(2,2) zeros(2,2)];
    y0.xypos = zeros(2, 1);
    R0.xypos = repelem(fOpt.sigma2RXYPosPVLSRS, 2);
    
    % covariance limiter
    H0.covlim = [eye(18) zeros(18, N.state-18)];
    R0.covlim = repelem([fOpt.sigma2RLimPos, fOpt.sigma2RLimOri, ...
                         fOpt.sigma2RLimPos, fOpt.sigma2RLimOri, ...
                         fOpt.sigma2RLimPos, fOpt.sigma2RLimOri], 3);
    
    %% D matrix initialization
    D0 = {}; d0 = {};
    % thigh length constraint
    D0.H2C = [eye(3,3); zeros(1,3)]; % homogenous to cartesian
    D0.H2CT = D0.H2C';
    D0.H2H = D0.H2C*D0.H2CT;
    D0.PV_p_LH = [0 body.PV_d/2 0 1]'; 
    D0.PV_p_RH = [0 -body.PV_d/2 0 1]';
    D0.LS_p_LK = [0 0 body.LS_d 1]';    
    D0.RS_p_RK = [0 0 body.RS_d 1]';
    D0.PV_p_LH_Circdot = point2fs(D0.PV_p_LH);
    D0.PV_p_RH_Circdot = point2fs(D0.PV_p_RH);
    D0.LS_p_LK_Circdot = point2fs(D0.LS_p_LK);
    D0.RS_p_RK_Circdot = point2fs(D0.RS_p_RK);
    d0.ltl = (body.LT_d).^2;
    d0.rtl = (body.RT_d).^2;
    
    % hinge knee joint
    D0.p_y = [0 1 0 0]'; 
    D0.p_y_Circdot = point2fs(D0.p_y);
    d0.lkh = 0;
    d0.rkh = 0;
    
    % knee range of motion
    d0.lkrom = 0;
    d0.rkrom = 0;
    
    %% Debug
    debug_dat.measUptPos = zeros(N.state, N.samples);
    debug_dat.measUptTilde = zeros(N.state, N.samples);
    
    %% Iteration
    se32vecIdxs = {1:3, 4:6, 7:9};
    u = zeros(18, N.samples);
    for k=2:(N.samples+1)
        kPast = k-1;
        %% Prediction update
        u(:,kPast) = [W_R_.PV(:,:,kPast)*B_a_.PV(kPast,:)' - g; ...
             W_R_.LS(:,:,kPast)*B_a_.LS(kPast,:)' - g; ...
             W_R_.RS(:,:,kPast)*B_a_.RS(kPast,:)' - g; ...
             W_R_.PV(:,:,kPast)*B_w_.PV(kPast,:)'; ...
             W_R_.LS(:,:,kPast)*B_w_.LS(kPast,:)'; ...
             W_R_.RS(:,:,kPast)*B_w_.RS(kPast,:)'];
        
        F.AdG = eye(N.state, N.state);
        bigphi = eye(N.state, N.state);
        F.curlyC = zeros(N.state, N.state);

        for i=1:N.se3StateList
            sname = se3StateList{i};
            bname = sname(end-1:end);

            xi = [x.tilde.vec(se32vecIdxs{i},kPast) ; 0; 0; 0];
            B_R_W = W_R_.(bname)(:,:,kPast)';
                           
            if knob.pred.posfromvel
                xi(1:3) = B_R_W*(xi(1:3) + 0.5*dt*u(se32vecIdxs{i}(1:3),kPast));
%                 xi(1:3) = B_R_W*xi(1:3);
                vname = sprintf('vel%s', sname(end-1:end));
                F.curlyC(idx.(sname)(1:3),idx.(vname)) = dt*B_R_W;
            else
                xi(1:3) = 0;
            end

            bigxi = vec2tran(xi*dt);
            bigphi(idx.(sname),idx.(sname)) = vec2jac(-xi*dt);
            F.AdG(idx.(sname),idx.(sname)) = tranAd(bigxi);
            
            x.hatPri.(sname)(:,:,k) = x.tilde.(sname)(:,:,kPast)*bigxi;
        end
        F.comb = F.AdG + bigphi*F.curlyC;
        x.hatPri.vec(:,k) = F.vec*x.tilde.vec(:,kPast) + G.vec*u(:,kPast);
        P.hatPri(:,:,k) = F.comb*P.tilde(:,:,kPast)*F.comb' + bigphi*Q.comb*bigphi';
        
        %% Measurement update
        H = {}; deltay = {}; R = {};
        N.meas_k = 0;
        
        H.ori_k = {};   R.ori_k = {};
        H.zpos_k = {};  R.zpos_k = {};

        for i=1:N.se3StateList
            sname = se3StateList{i};
            bname = sname(end-1:end);

            % Orientation
            if knob.meas.ori
                H.ori_k.(sname) = H0.ori.(sname);
                deltay.ori_k.(sname) = real(rot2vec(x.hatPri.(sname)(1:3,1:3,k)' ...
                                * W_R_.(bname)(:,:,kPast)));
                R.ori_k.(sname) = R0.ori.(sname);
            else
                H.ori_k.(sname) = []; 
                deltay.ori_k.(sname) = []; 
                R.ori_k.(sname) = [];
            end
            
            N.meas_k = N.meas_k + size(H.ori_k.(sname), 1);
            
            % zPos assumption
            if knob.meas.zpos.(bname) && step.(bname)(kPast) % step detected
                H.zpos_k.(sname) = zeros(1, N.state);
                H.zpos_k.(sname)(1, idx.(sname)) = ...
                    [0 0 1 0] * x.hatPri.(sname)(:,:,k) * H0.zposAnkPCircdot;
                deltay.zpos_k.(sname) = y0.zpos.(sname) ...
                    - [0 0 1 0]*x.hatPri.(sname)(:,:,k)*[0; 0; 0; 1];
                R.zpos_k.(sname) = R0.zpos.(sname);
            else
                H.zpos_k.(sname) = []; 
                deltay.zpos_k.(sname) = []; 
                R.zpos_k.(sname) = [];
            end
                       
            N.meas_k = N.meas_k + size(H.zpos_k.(sname), 1);
        end
        
        % XY pos
        if knob.meas.xyposPVLSRS
            H.xypos = zeros(2, N.state);
            H.xypos(:, idx.W_T_PV) = - H0.xyposDT * x.hatPri.W_T_PV(:,:,k) ...
                                     * H0.p0Circdot;
            H.xypos(:, idx.W_T_LS) = H0.xyposDT * x.hatPri.W_T_LS(:,:,k) ...
                                     * H0.p0Circdot / 2;
            H.xypos(:, idx.W_T_RS) = H0.xyposDT * x.hatPri.W_T_RS(:,:,k) ...
                                     * H0.p0Circdot / 2;
            deltay.xypos = y0.xypos - H0.xyposDT * ...
                ( x.hatPri.W_T_LS(:,:,k) * H0.p0 / 2 ...
                  + x.hatPri.W_T_RS(:,:,k) * H0.p0 / 2 ...
                  - x.hatPri.W_T_PV(:,:,k) * H0.p0); 
            R.xypos = R0.xypos;
        else
            H.xypos = []; deltay.xypos = []; R.xypos = [];
        end
        N.meas_k = N.meas_k + size(H.xypos, 1);
        
        % ZUPT
        if knob.meas.zupt && step.LS(kPast)
            H.lzupt_k = H0.lzupt; R.lzupt_k = R0.lzupt;
            deltay.lzupt_k = y0.lzupt - x.hatPri.vec(idx.vecLZupt,k); 
        else
            H.lzupt_k = []; deltay.lzupt_k = []; R.lzupt_k = [];
        end
        if knob.meas.zupt && step.RS(kPast)
            H.rzupt_k = H0.rzupt; R.rzupt_k = R0.rzupt;
            deltay.rzupt_k = y0.rzupt - x.hatPri.vec(idx.vecRZupt,k);
        else
            H.rzupt_k = []; deltay.rzupt_k = []; R.rzupt_k = [];
        end
        N.meas_k = N.meas_k + size(H.lzupt_k, 1);
        N.meas_k = N.meas_k + size(H.rzupt_k, 1);
        
        if N.meas_k > 0
            H.comb = [H.ori_k.W_T_PV; H.zpos_k.W_T_PV; ...
                      H.ori_k.W_T_LS; H.zpos_k.W_T_LS; ...
                      H.ori_k.W_T_RS; H.zpos_k.W_T_RS; ...
                      H.xypos; H.lzupt_k; H.rzupt_k];
            deltay.comb = [deltay.ori_k.W_T_PV; deltay.zpos_k.W_T_PV; ...
                           deltay.ori_k.W_T_LS; deltay.zpos_k.W_T_LS; ...
                           deltay.ori_k.W_T_RS; deltay.zpos_k.W_T_RS; ...
                           deltay.xypos; deltay.lzupt_k; deltay.rzupt_k];
            R.comb = diag([R.ori_k.W_T_PV R.zpos_k.W_T_PV ...
                           R.ori_k.W_T_LS R.zpos_k.W_T_LS ...
                           R.ori_k.W_T_RS R.zpos_k.W_T_RS ...
                           R.xypos R.lzupt_k R.rzupt_k ]);
                
            K = P.hatPri(:,:,k)*H.comb'/(H.comb*P.hatPri(:,:,k)*H.comb' + R.comb);
            bigphi = eye(N.state, N.state);
            measUpt = K*(deltay.comb);
            debug_dat.measUptPos(:,kPast) = measUpt;
            
            for i=1:N.se3StateList
                sname = se3StateList{i};
                x.hatPos.(sname)(:,:,k) = x.hatPri.(sname)(:,:,k)*...
                                            vec2tran(measUpt(idx.(sname)));
                bigphi(idx.(sname), idx.(sname)) = vec2jac(-measUpt(idx.(sname)));
            end
            x.hatPos.vec(:,k) = x.hatPri.vec(:,k) + measUpt(idx.vec);
            
            if knob.meas.covlim
                H.comb2 = [H.comb; H0.covlim];
                R.comb2 = diag([diag(R.comb)' R0.covlim]);
                K2 = P.hatPri(:,:,k)*H.comb2'/(H.comb2*P.hatPri(:,:,k)*H.comb2' + R.comb2);
                P.hatPos(:,:,k) = bigphi*(I_N-K2*H.comb2)*P.hatPri(:,:,k)*bigphi';
            else
                P.hatPos(:,:,k) = bigphi*(I_N-K*H.comb)*P.hatPri(:,:,k)*bigphi';
            end
        else
            for i=1:N.se3StateList
                sname = se3StateList{i};
                x.hatPos.(sname)(:,:,k) = x.hatPri.(sname)(:,:,k);
            end
            x.hatPos.vec(:,k) = x.hatPri.vec(:,k);
            P.hatPos(:,:,k) = P.hatPri(:,:,k);
        end

        %% Constraint update
        if knob.cstr.on
            n_LT = D0.H2CT*(x.hatPos.W_T_PV(:,:,k)*D0.PV_p_LH - ...
                                x.hatPos.W_T_LS(:,:,k)*D0.LS_p_LK);
            n_RT = D0.H2CT*(x.hatPos.W_T_PV(:,:,k)*D0.PV_p_RH - ...
                                x.hatPos.W_T_RS(:,:,k)*D0.RS_p_RK);
            N.cstr_k = 0;
            if knob.cstr.thighlength
                N.cstr_k = N.cstr_k + 2;
            end
            if knob.cstr.hingeknee
                N.cstr_k = N.cstr_k + 2;
            end
            if knob.cstr.kneerom
                alphaLK = atan2(-dot(n_LT, x.hatPos.W_T_LS(1:3,3,k)), ...
                                -dot(n_LT, x.hatPos.W_T_LS(1:3,1,k))) + 0.5*pi;
                alphaRK = atan2(-dot(n_RT, x.hatPos.W_T_RS(1:3,3,k)), ...
                                -dot(n_RT, x.hatPos.W_T_RS(1:3,1,k))) + 0.5*pi;
                if (alphaLK < fOpt.alphaLKmin) || (alphaLK > fOpt.alphaLKmax)
                    N.cstr_k = N.cstr_k + 1;
                end
                if (alphaRK < fOpt.alphaRKmin) || (alphaRK > fOpt.alphaRKmax)
                    N.cstr_k = N.cstr_k + 1;
                end
            end
            
            D = zeros(N.cstr_k, N.state);
            d = zeros(N.cstr_k, 1);
            dhat = zeros(N.cstr_k, 1);
            N.cstr_k = 0;
            
            % thigh length constraint
            if knob.cstr.thighlength
                D(N.cstr_k+1,idx.W_T_PV) = 2*n_LT'*D0.H2CT*x.hatPos.W_T_PV(:,:,k)*...
                                    D0.PV_p_LH_Circdot;
                D(N.cstr_k+1,idx.W_T_LS) = -2*n_LT'*D0.H2CT*x.hatPos.W_T_LS(:,:,k)*...
                                    D0.LS_p_LK_Circdot;
                D(N.cstr_k+2,idx.W_T_PV) = 2*n_RT'*D0.H2CT*x.hatPos.W_T_PV(:,:,k)*...
                                    D0.PV_p_RH_Circdot;
                D(N.cstr_k+2,idx.W_T_RS) = -2*n_RT'*D0.H2CT*x.hatPos.W_T_RS(:,:,k)*...
                                    D0.RS_p_RK_Circdot;
                d(N.cstr_k+(1:2),:) = [d0.ltl; d0.rtl];
                dhat(N.cstr_k+(1:2),:) = [n_LT'*n_LT; n_RT'*n_RT];
                N.cstr_k = N.cstr_k + 2;
            end
            
            % hinge knee joint constraint
            if knob.cstr.hingeknee
                W_r_LS_y = x.hatPos.W_T_LS(:,:,k) * D0.p_y;
                W_r_RS_y = x.hatPos.W_T_RS(:,:,k) * D0.p_y;
                D(N.cstr_k+1,idx.W_T_PV) = W_r_LS_y' * D0.H2H ...
                                  * x.hatPos.W_T_PV(:,:,k) * D0.PV_p_LH_Circdot;
                D(N.cstr_k+1,idx.W_T_LS) = - W_r_LS_y' * D0.H2H ...
                                  * x.hatPos.W_T_LS(:,:,k)*D0.LS_p_LK_Circdot ...
                                  + n_LT' * D0.H2CT * x.hatPos.W_T_LS(:,:,k) ...
                                  * D0.p_y_Circdot;
                D(N.cstr_k+2,idx.W_T_PV) = W_r_RS_y'*D0.H2H ...
                                  * x.hatPos.W_T_PV(:,:,k) * D0.PV_p_RH_Circdot;
                D(N.cstr_k+2,idx.W_T_RS) = - W_r_RS_y'*D0.H2H ...
                                  * x.hatPos.W_T_RS(:,:,k)*D0.RS_p_RK_Circdot ...
                                  + n_RT'*D0.H2CT*x.hatPos.W_T_RS(:,:,k) ...
                                  * D0.p_y_Circdot;
                d(N.cstr_k+(1:2),:) = [d0.lkh; d0.rkh];
                dhat(N.cstr_k+(1:2),:) = [(W_r_LS_y)'*D0.H2C*n_LT;
                               (W_r_RS_y)'*D0.H2C*n_RT];
                N.cstr_k = N.cstr_k + 2;
            end
            
            % knee range of motion
            if knob.cstr.kneerom
                if (alphaLK < fOpt.alphaLKmin) || (alphaLK > fOpt.alphaLKmax)
                    alphaLK2 = min(max(alphaLK, fOpt.alphaLKmin), fOpt.alphaLKmax);
                    a = [-sin(alphaLK2-pi/2); 0; cos(alphaLK2-pi/2); 0];
                    Ta = x.hatPos.W_T_LS(:,:,k) * a;

                    D(N.cstr_k+1,idx.W_T_PV) = Ta' * D0.H2H ...
                            * x.hatPos.W_T_PV(:,:,k) * D0.PV_p_LH_Circdot;
                    D(N.cstr_k+1,idx.W_T_LS) = n_LT' * D0.H2CT ...
                            * x.hatPos.W_T_LS(:,:,k) * point2fs(a) ...
                        - Ta' * D0.H2H * x.hatPos.W_T_LS(:,:,k) ...
                            * D0.LS_p_LK_Circdot;
                    d(N.cstr_k+1,:) = d0.lkrom;
                    dhat(N.cstr_k+1,:) = n_LT' * D0.H2CT * Ta;

                    N.cstr_k = N.cstr_k + 1;
                end
                if (alphaRK < fOpt.alphaRKmin) || (alphaRK > fOpt.alphaRKmax)
                    alphaRK2 = min(max(alphaRK, fOpt.alphaRKmin), fOpt.alphaRKmax);
                    a = [-sin(alphaRK2-pi/2); 0; cos(alphaRK2-pi/2); 0];
                    Ta = x.hatPos.W_T_RS(:,:,k) * a;

                    D(N.cstr_k+1,idx.W_T_PV) = Ta' * D0.H2H ...
                            * x.hatPos.W_T_PV(:,:,k) * D0.PV_p_RH_Circdot;
                    D(N.cstr_k+1,idx.W_T_RS) = n_RT' * D0.H2CT ...
                            * x.hatPos.W_T_RS(:,:,k) * point2fs(a) ...
                        - Ta' * D0.H2H * x.hatPos.W_T_RS(:,:,k) ...
                            * D0.RS_p_RK_Circdot;
                    d(N.cstr_k+1,:) = d0.rkrom;
                    dhat(N.cstr_k+1,:) = n_RT' * D0.H2CT * Ta;

                    N.cstr_k = N.cstr_k + 1;
                end
            end
        
            K = P.hatPos(:,:,k)*D'/(D*P.hatPos(:,:,k)*D');
            if knob.cstr.Pupdate
                P.tilde(:,:,k) = (I_N-K*D)*P.hatPos(:,:,k);
            else
                P.tilde(:,:,k) = P.hatPos(:,:,k);
            end
            
            measUpt = K*(d-dhat);
            debug_dat.measUptTilde(:,kPast) = measUpt;
                        
            for i=1:N.se3StateList
                sname = se3StateList{i};
                x.tilde.(sname)(:,:,k) = x.hatPos.(sname)(:,:,k)*...
                                            vec2tran(measUpt(idx.(sname)));
            end
            x.tilde.vec(:,k) = x.hatPos.vec(:,k) + measUpt(idx.vec);
        else
            for i=1:N.se3StateList
                sname = se3StateList{i};
                x.tilde.(sname)(:,:,k) = x.hatPos.(sname)(:,:,k);
            end
            x.tilde.vec(:,k) = x.hatPos.vec(:,k);
            P.tilde(:,:,k) = P.hatPos(:,:,k);
        end
    end
    
    for i=1:N.se3StateList
        sname = se3StateList{i};
        x.hatPri.(sname) = x.hatPri.(sname)(:,:,2:end);
        x.hatPos.(sname) = x.hatPos.(sname)(:,:,2:end);
        x.tilde.(sname) = x.tilde.(sname)(:,:,2:end);
    end
    x.hatPri.vec = x.hatPri.vec(:,2:end);
    x.hatPos.vec = x.hatPos.vec(:,2:end);
    x.tilde.vec = x.tilde.vec(:,2:end);
    P.hatPri = P.hatPri(:,:,2:end);
    P.hatPos = P.hatPos(:,:,2:end);
    P.tilde = P.tilde(:,:,2:end);
    
    xtilde = x.tilde;
    
    LTIBz = squeeze(x.tilde.W_T_LS(1:3,3,:))';
    RTIBz = squeeze(x.tilde.W_T_RS(1:3,3,:))';
    PELVy = squeeze(x.tilde.W_T_PV(1:3,2,:))';
    debug_dat.LFEO = squeeze(x.tilde.W_T_LS(1:3,4,:))' + body.LS_d * LTIBz;
    debug_dat.RFEO = squeeze(x.tilde.W_T_RS(1:3,4,:))' + body.RS_d * RTIBz;
    debug_dat.LFEP = squeeze(x.tilde.W_T_PV(1:3,4,:))' + body.PV_d/2 * PELVy;
    debug_dat.RFEP = squeeze(x.tilde.W_T_PV(1:3,4,:))' - body.PV_d/2 * PELVy;
    
    R_LFEM = zeros(3,3,N.samples);
    R_RFEM = zeros(3,3,N.samples);
    
    LFEM_z = (debug_dat.LFEP-debug_dat.LFEO)'; 
    LFEM_y = squeeze(x.tilde.W_T_LS(1:3,2,:));
    LFEM_x = cross(LFEM_y, LFEM_z);
    R_LFEM(:,3,:) = LFEM_z ./ vecnorm(LFEM_z, 2, 1);
    R_LFEM(:,2,:) = LFEM_y ./ vecnorm(LFEM_y, 2, 1);
    R_LFEM(:,1,:) = LFEM_x ./ vecnorm(LFEM_x, 2, 1);
    RFEM_z = (debug_dat.RFEP-debug_dat.RFEO)';
    RFEM_y = squeeze(x.tilde.W_T_RS(1:3,2,:));
    RFEM_x = cross(RFEM_y, RFEM_z);
    R_RFEM(:,3,:) = RFEM_z ./ vecnorm(RFEM_z, 2, 1);
    R_RFEM(:,2,:) = RFEM_y ./ vecnorm(RFEM_y, 2, 1);
    R_RFEM(:,1,:) = RFEM_x ./ vecnorm(RFEM_x, 2, 1);
    debug_dat.qLTH = rotm2quat(R_LFEM);
    debug_dat.qRTH = rotm2quat(R_RFEM);
    
    debug_dat.u = u;
    debug_dat.x.hatPri = x.hatPri;
    debug_dat.x.hatPos = x.hatPos;
    debug_dat.x.tilde = x.tilde;
    debug_dat.P.hatPri = P.hatPri;
    debug_dat.P.hatPos = P.hatPos;
    debug_dat.P.tilde = P.tilde;
    debug_dat.fs = fOpt.fs;
    debug_dat.body = body;
end