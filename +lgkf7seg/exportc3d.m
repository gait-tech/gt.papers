function [acq, acqDebug] = exportc3d(xtilde, debugData, fname, ...
    actBody, bIsStat, dataS, dataSIdx, simPVIMU)
% Export output of LG CEKF to c3d file
%
% :param xtilde: initial state in the GFR
% :param debug_dat: initial covariance
% :param fname: output file will be named <fname>.c3d and <fname>-Debug.c3d
% :param actBody: [Optional] reference grBody. Default: null
% :param bIsStat: [Optional] struct containing if step detect. Must be same size as xtilde.
% :param dataS: loaded mocapdb.XsensBody. Size may be bigger than xtilde.
% :param dataSIdx: Index of dataS to be used. Must be same size as xtilde.
% :param simPVIMU: simulated pelvis imu measurements
%
% :returns: acq - handle pointer to new btk c3d
%
% .. Author: - Luke Wicent Sy (GSBME, 2019 Nov 26)

if nargin <= 6 || dataSIdx == -1
    dataSIdx = 1:dataS.nSamples;    
end
if nargin <= 7
    simPVIMU = [];
end

exportDebug = nargout > 1;
if exportDebug
    [estBody, estBodyDebug] = papers.lgkf7seg.exportGrBody(xtilde, debugData, false);
else
    estBody = papers.lgkf7seg.exportGrBody(xtilde, debugData, false);
end
if actBody.nSamples > estBody.nSamples
    actBody = actBody.getSubset(1:estBody.nSamples);
    idx0 = 1:estBody.nSamples;
elseif actBody.nSamples < estBody.nSamples
    estBody = estBody.getSubset(1:actBody.nSamples);
    estBodyDebug = estBodyDebug.getSubset(1:actBody.nSamples*3);
    idx0 = 1:actBody.nSamples;
else
    idx0 = 1:estBody.nSamples;
end
idx0Dbg = 1:idx0(end)*3;

eMarkers = struct();
sensors = struct();

bodyKeyList = struct();
bodyKeyList.all = [ ...
   struct('ks', 'PELV', 'kp', 'MP', 'x', 'Pelvis', 'xs', 'qHips', 'vs', 'qRPV', 'vp', 'MIDPEL', 'alias', 'PV'), ...
   struct('ks', 'LFT', 'kp', 'LF', 'x', 'L_Foot', 'xs', 'qLeftFoot', 'vs', 'qLFT', 'vp', 'LTIO', 'alias', 'LF'), ...
   struct('ks', 'RFT', 'kp', 'RF', 'x', 'R_Foot', 'xs', 'qRightFoot', 'vs', 'qRFT', 'vp', 'RTIO', 'alias', 'RF'), ...
];
%    struct('ks', 'LTIB', 'kp', 'LA', 'x', 'L_LowLeg', 'xs', 'qLeftLeg', 'vs', 'qLSK', 'vp', 'LTIO', 'alias', 'LS'), ...
%    struct('ks', 'RTIB', 'kp', 'RA', 'x', 'R_LowLeg', 'xs', 'qRightLeg', 'vs', 'qRSK', 'vp', 'RTIO', 'alias', 'RS'), ...

for i = ["ks", "kp", "x", "vs", "vp", "alias"]
    bodyKeyList.(i) = arrayfun(@(x) x.(i), bodyKeyList.all, 'UniformOutput', false);
end

estBodyRel = estBody.changeRefFrame('MIDPEL');
if ~isempty(actBody)
    estBody2 = estBodyRel.toWorldFrame(actBody.MIDPEL, estBody.qRPV);
    actBodyRel = actBody.changeRefFrame('MIDPEL');
    actBody2 = actBodyRel.toWorldFrame(actBody.MIDPEL, actBody.qRPV);
end
    
%% Acceleration, velocity, angular velocity
% estimator
if ~isempty(dataS)
    dataS = dataS.getSubset(dataSIdx);
    
    % acc, angvel
%     sensors = dataS.exportRawMeasurementAsStruct(bodyKeyList.x, bodyKeyList.alias);
    sensors = dataS.exportRawMeasurementAsStruct({'Pelvis', ...
        'L_LowLeg', 'R_LowLeg', 'L_Foot', 'R_Foot'}, ...
        {'PV', 'LS', 'RS', 'LF', 'RF'});
end
vLFIMU2LA = quatrotate(quatconj(estBody.qLFT), [0 0 -0.5]);
vRFIMU2RA = quatrotate(quatconj(estBody.qRFT), [0 0 -0.5]);
sensors.W_LFAcc2 = sensors.W_LFAcc + cross(sensors.W_LFGyr, cross(sensors.W_LFGyr, vLFIMU2LA));
sensors.W_RFAcc2 = sensors.W_RFAcc + cross(sensors.W_RFGyr, cross(sensors.W_RFGyr, vRFIMU2RA));
sensors.W_LFAcc2Ref = sensors.W_LFAcc + cross(sensors.W_LFGyr, cross(sensors.W_LFGyr, vLFIMU2LA));
sensors.W_RFAcc2Ref = sensors.W_RFAcc + cross(sensors.W_RFGyr, cross(sensors.W_RFGyr, vRFIMU2RA));

% vel
sensors.W_PVVel = xtilde.vec(1:3,idx0)';
sensors.W_LFVel = xtilde.vec(4:6,idx0)';
sensors.W_RFVel = xtilde.vec(7:9,idx0)';

% simpv
if ~isempty(simPVIMU)
    sensors.W_PVAcc2IMUSim = quatrotate(quatconj(simPVIMU.ori), simPVIMU.acc) - [0 0 9.81];
    sensors.W_PVAVel2IMUSim = quatrotate(quatconj(simPVIMU.ori), simPVIMU.gyr);
    sensors.B_PVAcc2IMUSim = quatrotate(simPVIMU.ori, sensors.W_PVAcc2IMUSim);
    sensors.B_PVAVel2IMUSim = simPVIMU.gyr;
end

% global yaw angle
sensors.W_PVYaw = atan2(squeeze(xtilde.W_T_PV(2,1,:)), squeeze(xtilde.W_T_PV(1,1,:)));
sensors.W_LFYaw = atan2(squeeze(xtilde.W_T_LF(2,3,:)), squeeze(xtilde.W_T_LF(1,3,:)));
sensors.W_RFYaw = atan2(squeeze(xtilde.W_T_RF(2,3,:)), squeeze(xtilde.W_T_RF(1,3,:)));
sensors.W_LRFAveYaw = (sensors.W_LFYaw + sensors.W_RFYaw)/2;

% reference
if ~isempty(actBody)
    acc = actBody.calcJointAcc(bodyKeyList.vp);
    vel = actBody.calcJointVel(bodyKeyList.vp);
    angvel = actBody.calcSegAngVel(bodyKeyList.vs, 'W');
    for i=bodyKeyList.all
        sensors.(sprintf('W_%sAccRef', i.alias)) = acc.(i.vp);
        sensors.(sprintf('W_%sVelRef', i.alias)) = vel.(i.vp);
        sensors.(sprintf('W_%sAVelRef', i.alias)) = angvel.(i.vs);
    end
    
    sensors.ePos = estBody2.calcDPos(actBody2);
    sensors.eOri = estBody2.calcDOri(actBody2);
    
    W_T_PV = quat2rotm(actBody2.qRPV);
    W_T_LF = quat2rotm(actBody2.qLFT);
    W_T_RF = quat2rotm(actBody2.qRFT);
    sensors.W_PVYawRef = atan2(squeeze(W_T_PV(2,1,:)), squeeze(W_T_PV(1,1,:)));
    sensors.W_LFYawRef = atan2(squeeze(W_T_LF(2,3,:)), squeeze(W_T_LF(1,3,:)));
    sensors.W_RFYawRef = atan2(squeeze(W_T_RF(2,3,:)), squeeze(W_T_RF(1,3,:)));
    sensors.W_LRFAveYawRef = (sensors.W_LFYaw + sensors.W_RFYaw)/2;
end

%% covariance
dataSIdx = {{'PV',  1: 3,  4: 6, 19:21, 28:30}, ...
       {'LF',  7: 9, 10:12, 22:24, 31:33}, ...
       {'RF', 13:15, 16:18, 25:27, 34:36}};
suffix = '';
if size(xtilde.vec, 1) > 9
    idx2 = {'Pos', 'Ori', 'Vel', 'AVel'};
else
    idx2 = {'Pos', 'Ori', 'Vel'};
end
[nState, ~, ~] = size(debugData.P.tilde);
nSamples = length(idx0);
for i=1:3
    buf = zeros(nSamples, nState);
    for j=idx0
        buf(j,:) = diag(debugData.P.tilde(:,:,j));
    end
    for j=1:length(idx2)
        sname = sprintf('P%s%s%s', idx2{j}, dataSIdx{i}{1}, suffix);
        sensors.(sname) = buf(:,dataSIdx{i}{j+1});
%         sname2 = sprintf('vPos%s%s%s', idx2{j}, idx{i}{1}, suffix);
%         sensors.(sname2) = state2.measUptPos(idx{i}{j+1}, :)';
%         sname3 = sprintf('vTilde%s%s%s', idx2{j}, idx{i}{1}, suffix);
%         sensors.(sname3) = state2.measUptTilde(idx{i}{j+1}, :)';
    end
end

if isfield(bIsStat, 'LTIO')
    bIsStat = struct('LA', bIsStat.LTIO, 'RA', bIsStat.RTIO);
end
acq = estBody.exportc3d(sprintf('%s.c3d', fname), sensors, actBody, ...
                        bIsStat.LA(idx0), bIsStat.RA(idx0), eMarkers);
              
%% Export debug
if exportDebug
    actBodyDebug = actBody.repelem(3);
    
    %% Acceleration, velocity, angular velocity (debug)
    sensorsDebug = struct();
    fn = fieldnames(sensors);
    for i=1:length(fn)
        j = fn{i};
        sensorsDebug.(j) = repelem(sensors.(j), 3, 1);
    end
    
    % vel
    n2 = estBody.nSamples*3;           
    sensorsDebug.W_PVVel(1:3:n2, :) = debugData.x.hatPri.vec(1:3,idx0)';
    sensorsDebug.W_PVVel(2:3:n2, :) = debugData.x.hatPos.vec(1:3,idx0)';
    sensorsDebug.W_LFVel(1:3:n2, :) = debugData.x.hatPri.vec(4:6,idx0)';
    sensorsDebug.W_LFVel(2:3:n2, :) = debugData.x.hatPos.vec(4:6,idx0)';
    sensorsDebug.W_RFVel(1:3:n2, :) = debugData.x.hatPri.vec(7:9,idx0)';
    sensorsDebug.W_RFVel(2:3:n2, :) = debugData.x.hatPos.vec(7:9,idx0)';
    sensorsDebug = rmfield(sensorsDebug, {'ePos', 'eOri'});
    
    % yaw average
    sensorsDebug.W_PVYaw(1:3:n2, :) = atan2(squeeze(debugData.x.hatPri.W_T_PV(2,1,:)), ...
                                            squeeze(debugData.x.hatPri.W_T_PV(1,1,:)));
    sensorsDebug.W_PVYaw(2:3:n2, :) = atan2(squeeze(debugData.x.hatPos.W_T_PV(2,1,:)), ...
                                            squeeze(debugData.x.hatPos.W_T_PV(1,1,:)));
    sensorsDebug.W_LFYaw(1:3:n2, :) = atan2(squeeze(debugData.x.hatPri.W_T_LF(2,3,:)), ...
                                            squeeze(debugData.x.hatPri.W_T_LF(1,3,:)));
    sensorsDebug.W_LFYaw(2:3:n2, :) = atan2(squeeze(debugData.x.hatPos.W_T_LF(2,3,:)), ...
                                            squeeze(debugData.x.hatPos.W_T_LF(1,3,:)));
    sensorsDebug.W_RFYaw(1:3:n2, :) = atan2(squeeze(debugData.x.hatPri.W_T_RF(2,3,:)), ...
                                            squeeze(debugData.x.hatPri.W_T_RF(1,3,:)));
    sensorsDebug.W_RFYaw(2:3:n2, :) = atan2(squeeze(debugData.x.hatPos.W_T_RF(2,3,:)), ...
                                            squeeze(debugData.x.hatPos.W_T_RF(1,3,:)));
    sensorsDebug.W_LRFAveYaw = (sensorsDebug.W_LFYaw + sensorsDebug.W_RFYaw)/2;
    
    %% step detection
    bIsStatDebug = struct();
    fn = fieldnames(bIsStat);
    for i=1:length(fn)
        j = fn{i};
        bIsStatDebug.(j) = repelem(bIsStat.(j), 3, 1);
    end
    acqDebug = estBodyDebug.exportc3d(sprintf('%s-Debug.c3d', fname), ...
                    sensorsDebug, actBodyDebug, ...
                    bIsStatDebug.LA(idx0Dbg), bIsStatDebug.RA(idx0Dbg));
end

end