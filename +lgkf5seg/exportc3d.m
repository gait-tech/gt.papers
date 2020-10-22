function acq = exportc3d(xtilde, debugData, fname, ...
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

    estBody = papers.lgkf5seg.exportGrBody(xtilde, debugData, false);
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
        actBody.qLFT = []; actBody.qRFT = [];
        actBody.LTOE = []; actBody.RTOE = [];
        estBody2 = estBodyRel.toWorldFrame(actBody.MIDPEL, estBody.qRPV);
        actBodyRel = actBody.changeRefFrame('MIDPEL');
        actBody2 = actBodyRel.toWorldFrame(actBody.MIDPEL, actBody.qRPV);
    end

    %% Acceleration, velocity, angular velocity
    % estimator
    if ~isempty(dataS)
        dataS = dataS.getSubset(dataSIdx);

        % acc, angvel
        sensors = dataS.exportRawMeasurementAsStruct(bodyKeyList.x, bodyKeyList.alias);
    end
    % vel
    sensors.W_PVVel = xtilde.vec(1:3,:)';
    sensors.W_LSVel = xtilde.vec(4:6,:)';
    sensors.W_RSVel = xtilde.vec(7:9,:)';

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
    dataSIdx = {{'PV',  1: 3,  4: 6, 19:21, 28:30}, ...
           {'LS',  7: 9, 10:12, 22:24, 31:33}, ...
           {'RS', 13:15, 16:18, 25:27, 34:36}};
    suffix = '';
    if size(xtilde.vec, 1) > 9
        idx2 = {'Pos', 'Ori', 'Vel', 'AVel'};
    else
        idx2 = {'Pos', 'Ori', 'Vel'};
    end
    [nState, ~, nSamples] = size(debugData.P.tilde);
    for i=1:3
        buf = zeros(nSamples, nState);
        for j=1:nSamples
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
                            bIsStat.LA, bIsStat.RA, eMarkers);

end