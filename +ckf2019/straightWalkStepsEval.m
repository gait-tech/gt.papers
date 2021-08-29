dir = 'data/neura-sparse01/';
instruction = readtable('+papers/+ckf2019/neura-eval-step.csv');
expression = '^(?<dir>.+)/(?<ns>[^-]+)-(?<subj>[^-]+)-(?<act>.+)-(?<ns2>.+)-(?<algo>.+)\+Aw.+.mat$';
targetSeg = {'qLTH', 'qRTH'};

report = table();
for i = 1:size(instruction, 1)
    n = table2struct(instruction(i, :));
    n.cmd = lower(n.cmd);
    if strcmp(n.cmd, 'load')
        tokens = regexp(n.arg0, expression, 'names');

        name = sprintf("%s-%s-%s", tokens.ns, tokens.subj, tokens.act);
        data = mocapdb.loadNeuRaTrialData(dir, tokens.subj, tokens.act, false);
        actBody = data.dataV.togrBody(1:data.dataV.nSamples, ...
            {'name', name, 'oriUnit', 'deg', 'lnSymbol', '-', 'ptSymbol', '*', ...
             'fs', data.dataV.fs, 'xyzColor', {'m', 'y', 'c'}});
        load(n.arg0);
        fprintf("Loaded %s\n", n.arg0);
    elseif strcmp(n.cmd, 'eval')
        actBody2 = actBody.getSubset(n.startIndex:n.endIndex);
        
        estBody2a = estBody.getSubset(n.startIndex:n.endIndex);
        estBodyRel = estBody2a.changeRefFrame('MIDPEL');
        estBody2 = estBodyRel.toWorldFrame(actBody2.MIDPEL, estBody2a.qRPV);
        
        % calculate kinematic gait parameters against actBody (vicon)
        r_buf = estBody2.diffRMSEandMean(actBody2, true, targetSeg);
        r_buf2 = struct('name', name, 'side', string(n.arg0), ...
                        'sIdx', n.startIndex, 'eIdx', n.endIndex);
                    
        paramFields = fieldnames(r_buf);
        for paramI = 1:length(paramFields)
            sname = paramFields{paramI};
            if strncmp(sname, n.arg0, 1)
                r_buf2.(extractAfter(sname, 1)) = r_buf.(sname);
            elseif strncmp(sname, sprintf('q%s', n.arg0), 2) && ~strncmp(sname, 'qRPV', 4)
                r_buf2.(sprintf('q%s', extractAfter(sname, 2))) = r_buf.(sname);
            end 
        end
        
        report = [report; struct2table(r_buf2)];
    end
end 