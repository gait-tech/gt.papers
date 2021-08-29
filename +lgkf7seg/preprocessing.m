function [imuData, gfrAcc, step] = preprocessing(dataS, dataV, B_p_TP, ...
                    calibSW2W, calibS2BD, revStepDetect, savefname)
    % Run preprocessing for this movement trial
    % 
    % - imuData contain multiple XsensBody of body wrt world frame
    % - gfrAcc contain multiple body inertial acceleration wrt world frame (w/out gravity)
    % - step is a table containing stepL and stepR
    % 
    % By default, imuData and gfrAcc are struct of modes (w__sv, w__v)
    %
    % :param dataS: loaded mocapdb.XsensBody
    % :param dataV: loaded mocapdb.ViconBody
    % :param B_p_TP: struct of target point wrt body segment. keys must be in :class:`+pelib.@grBody.grBody.TMap`
    % :param calibSW2W: mocapdb.XsensBody aligning sensor world to global world. E.g., fix ankle sensor yaw orientation 
    % :param calibS2B: mocapdb.XsensBody aligning sensor to body frame. If empty, infer from dataV
    % :param revStepDetect: manually reviewed table if stepL and stepR was detected
    % :param savefname: filepath to save .mat output/debug files (optional)
    % 
    % :return: imuData, gfrAcc, and step

    %% Inputs and Input Check
    validateattributes(dataS, {'mocapdb.XsensBody'}, {});
    validateattributes(dataV, {'mocapdb.ViconBody'}, {});
    validateattributes(calibSW2W, {'mocapdb.XsensBody'}, {});
    if nargin <= 6, savefname = ''; end
    
    %% Initialization   
    % Initialize other variables   
    imuData = {};
    gfrAcc = {};
    step = revStepDetect;
    
    segList = {'Pelvis', 'L_LowLeg', 'R_LowLeg', 'L_Foot', 'R_Foot'};
    dataS = dataS.getSegSubset(segList);
    calibSW2W = calibSW2W.getSegSubset(segList);
    % obtain the sensor to body frame offset from the initial vicon pose
    calibS2BV = dataS.adjustFrame(calibSW2W, [1 0 0 0], true).calcCalibSB(dataV.togrBody(2, {}), 1);
    % calibSW2W.qB * dataS.qB * calibS2B.qB for each body segment
    % dataB containts IMU measurement of body
    % dataB.Pelvis.ori = body in world frame
    % dataB.Pelvis.acc, gyr, mag = in body frame
    dataBV = dataS.adjustFrame(calibSW2W, calibS2BV.conj());
    
    %% orientation, body acceleration, angular velocity
    W__viconBody = dataV.togrBody(1:dataV.nSamples, {'name', 'act', 'oriUnit', 'deg', ...
                     'lnSymbol', '-', 'ptSymbol', '*', 'fs', dataV.fs, ...
                     'xyzColor', {'m', 'y', 'c'}});
    imuData.w__sv = dataBV;
    imuData.w__v = calcSparseSimIMU(W__viconBody, B_p_TP);
    imuData.raw = dataS.adjustFrame(calibSW2W, [1 0 0 0]);

    %% body acceleration (w/out gravity) in world frame
    gfrAcc.w__sv = dataBV.calcGfrAcc();
    gfrAcc.w__v = imuData.w__v.calcGfrAcc();
    
    if ~isempty(calibS2BD)
        calibS2BD = calibS2BD.getSegSubset(segList);
        dataBD = dataS.adjustFrame(calibSW2W, calibS2BD.conj());
        imuData.w__sd = dataBD;
        gfrAcc.w__sd = dataBD.calcGfrAcc();
    end
    
    %% Save processing
    if ~strcmp(savefname, '')
        save(savefname, 'W__viconBody', 'imuData', 'gfrAcc');
    end
end

function out = calcSparseSimIMU(grB, B_p_TP)
	% Calculate the sparse simulation imu of grBody
	%
	% :param grB: input grBody
    % :type grB: :class:`+pelib.@grBody`
    % :param B_p_TP: struct of target point wrt body segment. keys must be in :class:`+pelib.@grBody.grBody.TMap`
    % :type B_p_TP: struct
    %
	% :return: simulated IMU data of body segment
	% :rtype: :class:`+mocapdb.@XsensBody`
    %
	% .. Author: - Luke Sy (UNSW GSBME) - 2020 May 25
    if ~strcmp(grB.frame, 'world')
        error('grBody must be in world frame');
    end
    
    bodyList.all = [
        struct('xsnbody', 'Pelvis', 'grjoint', 'MIDPEL', 'grseg', 'qRPV'), ...
        struct('xsnbody', 'L_LowLeg', 'grjoint', 'LTIO', 'grseg', 'qLSK'), ...
        struct('xsnbody', 'R_LowLeg', 'grjoint', 'RTIO', 'grseg', 'qRSK'), ...
        struct('xsnbody', 'L_Foot', 'grjoint', 'LFT', 'grseg', 'qLFT'), ...
        struct('xsnbody', 'R_Foot', 'grjoint', 'RFT', 'grseg', 'qRFT'), ...
        ];
    
    bodyList.grjoint = struct();
    for i = bodyList.all
        if isfield(B_p_TP, i.grjoint)
            bodyList.grjoint.(i.grjoint) = B_p_TP.(i.grjoint);
        else
            bodyList.grjoint.(i.grjoint) = [0 0 0 1];
        end
    end
    bodyList.grseg = arrayfun(@(x) x.grseg, bodyList.all, 'UniformOutput', false);
    
    out = mocapdb.XsensBody('srcFileName', grB.name, ...
                            'nSamples', grB.nSamples, ...
                            'fs', grB.fs, 'frame', grB.frame);
    
    
    gfracc = grB.calcJointAcc(bodyList.grjoint);
    angvel = grB.calcSegAngVel(bodyList.grseg, 'B');
    
    for i = bodyList.all
        out.(i.xsnbody).ori = grB.(i.grseg);
        out.(i.xsnbody).acc = quatrotate(grB.(i.grseg), ...
                                         gfracc.(i.grjoint) + [0 0 9.81]);
        out.(i.xsnbody).gyr = angvel.(i.grseg);
    end
end