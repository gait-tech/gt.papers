function [acq, acqDebug] = exportc3d(xtilde, debugData, fname, ...
    actBody, bIsStat, dataS, dataSIdx)
% Export output of LG CEKF to c3d file
%
% :param xtilde: initial state in the GFR
% :param debug_dat: initial covariance
% :param fname: output file will be named <fname>.c3d and <fname>-Debug.c3d
% :param actBody: [Optional] reference grBody. Default: null
% :param bIsStat: [Optional] struct containing if step detect. Must be same size as xtilde.
% :param dataS: loaded mocapdb.XsensBody. Size may be bigger than xtilde.
% :param dataSIdx: Index of dataS to be used. Must be same size as xtilde.
%
% :returns: acq - handle pointer to new btk c3d
%
% .. Author: - Luke Wicent Sy (GSBME, 2019 Nov 26)

if nargin <= 6
    dataSIdx = 1:dataS.nSamples;
end

exportDebug = nargout > 1;
if exportDebug
    [estBody, estBodyDebug] = papers.ckfdist2020.exportGrBody(xtilde, debugData, false);
else
    estBody = papers.ckfdist2020.exportGrBody(xtilde, debugData, false);
end

eMarkers = struct();
sensors = struct();

bodyKeyList = struct();
bodyKeyList.all = [ ...
   struct('ks', 'PELV', 'kp', 'MP', 'x', 'Pelvis', 'xs', 'qHips', 'vs', 'qRPV', 'vp', 'MIDPEL', 'alias', 'PV'), ...
   struct('ks', 'LTIB', 'kp', 'LA', 'x', 'L_LowLeg', 'xs', 'qLeftLeg', 'vs', 'qLSK', 'vp', 'LTIO', 'alias', 'LS'), ...
   struct('ks', 'RTIB', 'kp', 'RA', 'x', 'R_LowLeg', 'xs', 'qRightLeg', 'vs', 'qRSK', 'vp', 'RTIO', 'alias', 'RS'), ...
];

for i = ["ks", "kp", "x", "vs", "vp", "alias"]
    bodyKeyList.(i) = arrayfun(@(x) x.(i), bodyKeyList.all, 'UniformOutput', false);
end

estBodyRel = estBody.changeRefFrame('MIDPEL');
if ~isempty(actBody)
    actBody = actBody.copy();
    actBody.qLFT = [];
    actBody.qRFT = [];
    
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
        'L_LowLeg', 'R_LowLeg'}, ...
        {'PV', 'LS', 'RS'});
end
% vel
sensors.W_PVVel = xtilde(:,4:6);
sensors.W_LAVel = xtilde(:,10:12);
sensors.W_RAVel = xtilde(:,16:18);

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
end

%% covariance
dataSIdx = {{'PV',  1: 3,  4: 6}, ...
       {'LA',  7: 9, 10:12}, ...
       {'RA', 13:15, 16:18}};
suffix = '';
idx2 = {'Pos', 'Vel'};

[nState, ~, nSamples] = size(debugData.cstrP);
for i=1:3
    buf = zeros(nSamples, nState);
    for j=1:nSamples
        buf(j,:) = diag(debugData.cstrP(:,:,j));
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
                        bIsStat.LA, bIsStat.RA, eMarkers);
              
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
    sensorsDebug.W_PVVel(1:3:n2, :) = debugData.predState(:,4:6);
    sensorsDebug.W_PVVel(2:3:n2, :) = debugData.zuptState(:,4:6);
    sensorsDebug.W_LAVel(1:3:n2, :) = debugData.predState(:,10:12);
    sensorsDebug.W_LAVel(2:3:n2, :) = debugData.zuptState(:,10:12);
    sensorsDebug.W_RAVel(1:3:n2, :) = debugData.predState(:,16:18);
    sensorsDebug.W_RAVel(2:3:n2, :) = debugData.zuptState(:,16:18);
    sensorsDebug = rmfield(sensorsDebug, {'ePos', 'eOri'});
    
    %% step detection
    bIsStatDebug = struct();
    fn = fieldnames(bIsStat);
    for i=1:length(fn)
        j = fn{i};
        bIsStatDebug.(j) = repelem(bIsStat.(j), 3, 1);
    end
    acqDebug = estBodyDebug.exportc3d(sprintf('%s-Debug.c3d', fname), ...
                    sensorsDebug, actBodyDebug, bIsStatDebug.LA, bIsStatDebug.RA);
end

end