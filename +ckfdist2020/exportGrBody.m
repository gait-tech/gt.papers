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
    idxEndIdx = find(any(isnan(xtilde), 2), 1);
    if isempty(idxEndIdx)
        idx = 1:length(xtilde(:,1));
    elseif idxEndIdx == 0
        idx = []; 
    else
        idx = 1:(idxEndIdx-1); 
    end
else
    idx = 1:length(xtilde(:,1));
end

out = pelib.grBody('name', 'est', 'posUnit', 'm', 'oriUnit', 'deg', ...
           'lnSymbol', '--', 'ptSymbol', 'o', 'frame', 'world', ...
           'xyzColor', {'r', 'g', 'b'}, 'fs', debugData.fs, ...
           'MIDPEL', xtilde(idx, 1:3), ...
           'LFEP', debugData.LFEP(idx, :), ...
           'LFEO', debugData.LFEO(idx, :), ...
           'LTIO', xtilde(idx, 7:9), ...
           'RFEP', debugData.RFEP(idx, :), ...
           'RFEO', debugData.RFEO(idx, :), ...
           'RTIO', xtilde(idx, 13:15), ...
           'qRPV', debugData.qRPV(idx, :), ...
           'qLTH', debugData.qLTH(idx, :), ...
           'qRTH', debugData.qRTH(idx, :), ...
           'qLSK', debugData.qLSK(idx, :), ...
           'qRSK', debugData.qRSK(idx, :));

% calculate outDebug only when asked
if nargout > 1
    outDebug = out.repelem(3);
    n2 = outDebug.nSamples;
    
    outDebug.MIDPEL(1:3:n2, :) = debugData.predState(:,1:3);
    outDebug.MIDPEL(2:3:n2, :) = debugData.zuptState(:,1:3);
    outDebug.LTIO(1:3:n2, :) = debugData.predState(:,7:9);
    outDebug.LTIO(2:3:n2, :) = debugData.zuptState(:,7:9);
    outDebug.RTIO(1:3:n2, :) = debugData.predState(:,13:15);
    outDebug.RTIO(2:3:n2, :) = debugData.zuptState(:,13:15);
            
    v = quat2rotm(outDebug.qLSK); v = squeeze(v(:,3,:))';
    outDebug.LFEO = outDebug.LTIO + outDebug.calcLShankLength(3)*v;
    v = quat2rotm(outDebug.qRSK); v = squeeze(v(:,3,:))';
    outDebug.RFEO = outDebug.RTIO + outDebug.calcRShankLength(3)*v;
    v = quat2rotm(outDebug.qRPV); v = squeeze(v(:,2,:))';
    outDebug.LFEP = outDebug.MIDPEL + outDebug.calcPelvisLength(3)/2*v;
    outDebug.RFEP = outDebug.MIDPEL - outDebug.calcPelvisLength(3)/2*v;

    v = zeros(3, 3, n2);
    z = (outDebug.LFEP-outDebug.LFEO)';
    z = z ./ vecnorm(z, 2, 1);
    v(:, 3, :) = reshape(z, 3, 1, []);
    y =  quat2rotm(outDebug.qLSK);
    v(:, 2, :) = y(:, 2, :);
    x = cross(v(:, 2, :), v(:, 3, :));
    x = x ./ vecnorm(x, 2, 1);
    v(:, 1, :) =  reshape(x, 3, 1, []);
    outDebug.qLTH = rotm2quat(v);

    v = zeros(3, 3, n2);
    z = (outDebug.RFEP-outDebug.RFEO)';
    z = z ./ vecnorm(z, 2, 1);
    v(:, 3, :) = reshape(z, 3, 1, []);
    y =  quat2rotm(outDebug.qRSK);
    v(:, 2, :) = y(:, 2, :);
    x = cross(v(:, 2, :), v(:, 3, :));
    x = x ./ vecnorm(x, 2, 1);
    v(:, 1, :) =  reshape(x, 3, 1, []);
    outDebug.qRTH = rotm2quat(v);
        
    outDebug.fs = outDebug.fs*3;
end