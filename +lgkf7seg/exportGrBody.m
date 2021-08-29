function [out, outDebug] = exportGrBody(xtilde, debugData, endAtFirstNaN)
% Export output of LG CEKF to grBody class
%
% :param xtilde: initial state in the GFR
% :param debugData: initial covariance
% :param endAtFirstNaN: [Optional] Boolean. If true, end at first nan state. Default: true
%
% :returns: [estBody, estoutDebug] - output and output debug grBody class
%
% .. Author: - Luke Wicent Sy (GSBME, 2019 Nov 26)

if nargin <= 2
    endAtFirstNaN = true;
end
if endAtFirstNaN
    idxEndIdx = find(any(isnan(xtilde.vec), 1), 1);
    if isempty(idxEndIdx)
        idx = 1:length(xtilde.vec(1,:));
    elseif idxEndIdx == 0
        idx = []; 
    else
        idx = 1:(idxEndIdx-1); 
    end
else
    idx = 1:length(xtilde.vec(1,:));
end

out = pelib.grBody('name', 'est', 'posUnit', 'm', 'oriUnit', 'deg', ...
   'lnSymbol', '--', 'ptSymbol', 'o', 'frame', 'world', ...
   'xyzColor', {'r', 'g', 'b'}, 'fs', debugData.fs, ...
   'nSamples', idx(end), ...
   'MIDPEL', squeeze(xtilde.W_T_PV(1:3,4,idx))', ...
   'LFEP', debugData.LFEP(idx, :), ...
   'LFEO', debugData.LFEO(idx, :), ...
   'LTIO', debugData.LTIO(idx, :), ...
   'LTOE', debugData.LTOE(idx, :), ...
   'RFEP', debugData.RFEP(idx, :), ...
   'RFEO', debugData.RFEO(idx, :), ...
   'RTIO', debugData.RTIO(idx, :), ...
   'RTOE', debugData.RTOE(idx, :), ...
   'qRPV', rotm2quat(xtilde.W_T_PV(1:3,1:3,idx)), ...
   'qLTH', debugData.qLTH(idx, :), 'qRTH', debugData.qRTH(idx, :), ...
   'qLSK', debugData.qLSK(idx, :), 'qRSK', debugData.qRSK(idx, :), ...
   'qLFT', rotm2quat(xtilde.W_T_LF(1:3,1:3,idx)), ...
   'qRFT', rotm2quat(xtilde.W_T_RF(1:3,1:3,idx)));

% calculate outDebug only when asked
if nargout > 1
    outDebug = out.repelem(3);
    n2 = outDebug.nSamples;
    body = debugData.body;
    
    outDebug.MIDPEL(1:3:n2, :) = squeeze(debugData.x.hatPri.W_T_PV(1:3,4,idx))';
    outDebug.MIDPEL(2:3:n2, :) = squeeze(debugData.x.hatPos.W_T_PV(1:3,4,idx))';
    outDebug.LTIO(1:3:n2, :) = squeeze(debugData.x.hatPri.W_T_LF(1:3,4,idx))';
    outDebug.LTIO(2:3:n2, :) = squeeze(debugData.x.hatPos.W_T_LF(1:3,4,idx))';
    outDebug.RTIO(1:3:n2, :) = squeeze(debugData.x.hatPri.W_T_RF(1:3,4,idx))';
    outDebug.RTIO(2:3:n2, :) = squeeze(debugData.x.hatPos.W_T_RF(1:3,4,idx))';
    outDebug.qRPV(1:3:n2, :) = rotm2quat(debugData.x.hatPri.W_T_PV(1:3,1:3,idx));
    outDebug.qRPV(2:3:n2, :) = rotm2quat(debugData.x.hatPos.W_T_PV(1:3,1:3,idx));
    outDebug.qLFT(1:3:n2, :) = rotm2quat(debugData.x.hatPri.W_T_LF(1:3,1:3,idx));
    outDebug.qLFT(2:3:n2, :) = rotm2quat(debugData.x.hatPos.W_T_LF(1:3,1:3,idx));
    outDebug.qRFT(1:3:n2, :) = rotm2quat(debugData.x.hatPri.W_T_RF(1:3,1:3,idx));
    outDebug.qRFT(2:3:n2, :) = rotm2quat(debugData.x.hatPos.W_T_RF(1:3,1:3,idx));
    
    % copy pasted from lgcekf7seg post processing
    post = struct();
    post.W_R_PV = quat2rotm(outDebug.qRPV);
    post.W_R_LF = quat2rotm(outDebug.qLFT);
    post.W_R_RF = quat2rotm(outDebug.qRFT);
    
    post.LFTy = squeeze(post.W_R_LF(1:3,2,:));
    post.LFTz = squeeze(post.W_R_LF(1:3,3,:))';
    post.RFTy = squeeze(post.W_R_RF(1:3,2,:));
    post.RFTz = squeeze(post.W_R_RF(1:3,3,:))';
    post.PELVy = squeeze(post.W_R_PV(1:3,2,:))';
    
    outDebug.LTOE = outDebug.LTIO + body.LF_d * post.LFTz;
    outDebug.RTOE = outDebug.RTIO + body.RF_d * post.RFTz;
    outDebug.LFEP = outDebug.MIDPEL + body.PV_d/2 * post.PELVy;
    outDebug.RFEP = outDebug.MIDPEL - body.PV_d/2 * post.PELVy;
    
    % calculate LFEO and RFEO
    post.v_LFEPLTIO = outDebug.LFEP - outDebug.LTIO;
    post.v_RFEPRTIO = outDebug.RFEP - outDebug.RTIO;
    post.d0_LFEPLTIO = vecnorm(post.v_LFEPLTIO, 2, 2);
    post.d0_RFEPRTIO = vecnorm(post.v_RFEPRTIO, 2, 2);
    post.d_LFEPLTIO = min(post.d0_LFEPLTIO, body.LT_d+body.LS_d);
    post.d_RFEPRTIO = min(post.d0_RFEPRTIO, body.RT_d+body.RS_d);
    % without real(), theta_LKNE will have a very small complex component
    post.theta_LKNE = real(acos((body.LT_d.^2 + body.LS_d.^2 - post.d_LFEPLTIO.^2)/(2*body.LT_d*body.LS_d)));
    post.theta_RKNE = real(acos((body.RT_d.^2 + body.RS_d.^2 - post.d_RFEPRTIO.^2)/(2*body.RT_d*body.RS_d)));
    post.theta_LSHK = asin(body.LT_d.*sin(post.theta_LKNE)./post.d_LFEPLTIO);
    post.theta_RSHK = asin(body.RT_d.*sin(post.theta_RKNE)./post.d_RFEPRTIO);
    
    post.qLFEPLFEO = axang2quat([post.LFTy' post.theta_LSHK]);
    post.qRFEPRFEO = axang2quat([post.RFTy' post.theta_RSHK]);
    outDebug.LFEO = outDebug.LTIO + quatrotate(quatconj(post.qLFEPLFEO), ...
                            post.v_LFEPLTIO./post.d0_LFEPLTIO.*body.LS_d);
    outDebug.RFEO = outDebug.RTIO + quatrotate(quatconj(post.qRFEPRFEO), ...
                            post.v_RFEPRTIO./post.d0_RFEPRTIO.*body.RS_d);
    
    post.W_R_LT = zeros(3,3,n2);
    post.W_R_RT = zeros(3,3,n2);
    post.LFEM_z = (outDebug.LFEP-outDebug.LFEO)';
    post.LFEM_z = post.LFEM_z ./ vecnorm(post.LFEM_z, 2, 1);
    post.W_R_LT(:,3,:) = post.LFEM_z;
    post.W_R_LT(:,2,:) = post.LFTy;
    post.W_R_LT(:,1,:) = cross(post.LFTy, post.LFEM_z); 
    post.RFEM_z = (outDebug.RFEP-outDebug.RFEO)';
    post.RFEM_z = post.RFEM_z ./ vecnorm(post.RFEM_z, 2, 1);
    post.W_R_RT(:,3,:) = post.RFEM_z;
    post.W_R_RT(:,2,:) = post.RFTy;
    post.W_R_RT(:,1,:) = cross(post.RFTy, post.RFEM_z); 
    outDebug.qLTH = rotm2quat(post.W_R_LT);
    outDebug.qRTH = rotm2quat(post.W_R_RT);
    
    post.W_R_LS = zeros(3,3,n2);
    post.W_R_RS = zeros(3,3,n2);  
    
    post.LSHK_z = (outDebug.LFEO-outDebug.LTIO)';
    post.LSHK_z = post.LSHK_z ./ vecnorm(post.LSHK_z, 2, 1);
    post.W_R_LS(:,3,:) = post.LSHK_z;
    post.W_R_LS(:,2,:) = post.LFTy;
    post.W_R_LS(:,1,:) = cross(post.LFTy, post.LSHK_z); 
    post.RSHK_z = (outDebug.RFEO-outDebug.RTIO)';
    post.RSHK_z = post.RSHK_z ./ vecnorm(post.RSHK_z, 2, 1);
    post.W_R_RS(:,3,:) = post.RSHK_z;
    post.W_R_RS(:,2,:) = post.RFTy;
    post.W_R_RS(:,1,:) = cross(post.RFTy, post.RSHK_z);
    outDebug.qLSK = rotm2quat(post.W_R_LS);
    outDebug.qRSK = rotm2quat(post.W_R_RS);
end