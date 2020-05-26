function [imuData, gfrAcc, step] = preprocessing(dataS, dataV, ...
                    calibYawFix, revStepDetect, savefname)
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
    % :param calibYawFix: mocapdb.XsensBody fix ankle sensor yaw orientation offset
    % :param revStepDetect: manually reviewed table if stepL and stepR was detected
    % :param savefname: filepath to save .mat output/debug files (optional)
    % 
    % :return: imuData, gfrAcc, and step

    %% Inputs and Input Check
    validateattributes(dataS, {'mocapdb.XsensBody'}, {});
    validateattributes(dataV, {'mocapdb.ViconBody'}, {});
    validateattributes(calibYawFix, {'mocapdb.XsensBody'}, {});
    if nargin <= 4, savefname = ''; end
    
    %% Initialization   
    % Initialize other variables
    fs = dataS.fs;
    
    imuData = {};
    gfrAcc = {};
    step = revStepDetect;
    
    segList = {'Pelvis', 'L_LowLeg', 'R_LowLeg'};
    grJointList = {'MIDPEL', 'LTIO', 'RTIO'};
    dataS = dataS.getSegSubset(segList);
    % In the CKF 2020 paper, only the ankle yaw offset angles were used.
    % Hence, pelvis yaw offset was set to identity
    calibYawFix = calibYawFix.getSegSubset({'L_LowLeg', 'R_LowLeg'});
    % apply yaw fix and obtain the sensor to body frame offset from the
    % initial vicon pose
    calibS2B = dataS.adjustFrame(calibYawFix, [1 0 0 0], true).calcCalibSB(dataV.togrBody(2, {}), 1); 
    % calibYawFix.qB * dataS.qB * calibS2B.qB for each body segment
    % dataB containts IMU measurement of body
    % dataB.Pelvis.ori = body in world frame
    % dataB.Pelvis.acc, gyr, mag = in body frame
    dataB = dataS.adjustFrame(calibYawFix, calibS2B.conj());
    
    %% orientation, body acceleration, angular velocity
    W__viconBody = dataV.togrBody(1:dataV.nSamples, {'name', 'act', 'oriUnit', 'deg', ...
                     'lnSymbol', '-', 'ptSymbol', '*', 'fs', dataV.fs, ...
                     'xyzColor', {'m', 'y', 'c'}});
    imuData.w__sv = dataB;
    imuData.w__v = calcSparseSimIMU(W__viconBody);
    imuData.raw = dataS.adjustFrame(calibYawFix, [1 0 0 0]);

    %% body acceleration (w/out gravity) in world frame
    gfrAcc.w__sv = dataB.calcGfrAcc();
    acc = W__viconBody.calcJointAcc(grJointList);
    for i = 1:3
        gfrAcc.w__v.(segList{i}) = acc.(grJointList{i});
    end
    
    %% Save processing
    if ~strcmp(savefname, '')
        save(savefname, 'W__viconBody', 'imuData', 'gfrAcc');
    end
end

function out = calcSparseSimIMU(grB)
	% Calculate the sparse simulation imu of grBody
	%
	% :param obj: this grBody
    %
	% :return: out mocapdb.XsensBody
	%
	% .. Author: - Luke Sy (UNSW GSBME) - 2020 May 25
    if ~strcmp(grB.frame, 'world')
        error('grBody must be in world frame');
    end
    
    bodyList.all = [
        struct('xsnbody', 'Pelvis', 'grjoint', 'MIDPEL', 'grseg', 'qRPV'), ...
        struct('xsnbody', 'L_LowLeg', 'grjoint', 'LTIO', 'grseg', 'qLSK'), ...
        struct('xsnbody', 'R_LowLeg', 'grjoint', 'RTIO', 'grseg', 'qRSK'), ...
        ];
    
    bodyList.grjoint = arrayfun(@(x) x.grjoint, bodyList.all, 'UniformOutput', false);
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