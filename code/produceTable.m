%% produceTable
% This script loads a blink data set into a MATLAB table variable. When
% run, it will aggregate data for a given subject(s) and parameter(s). It 
% will then produce a table which describes the mean R2 by puff pressure
% across subjects, the R of the test-retest correlation for slope and for 
% offset. The R2 of the habituation effect across trials, and the R the 
% test retest of the slope of that habituation effect.

%% set up parameters

% load file path
dataPath = fileparts(fileparts(mfilename('fullpath')));
spreadsheet1 ='UPenn Ipsi Summary_25ms_02062022.csv';
spreadsheet2 ='Upenn_Ipsilateral Afiles_clean_full.csv';

% choose subject and parameters
% only run subjects 14590, 14589, and 14588 for highest 3 PSI levels
subList = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591, ...
    14590, 14589, 14588, 14587, 14586};
varNamesToPlot1 = {'aucI', 'latencyI', 'timeUnderI', 'openTimeI', 'initVelocityI', ...
     'closeTimeI', 'maxClosingVelocityI', 'maxOpeningVelocityI', 'blinkRate'};
varNamesToPlot2 = {'auc', 'latency', 'timeUnder20', 'openTime', 'initialVelocity', ...
     'closeTime', 'maxClosingVelocity', 'maxOpeningVelocity', 'blinkRate'};
highestOnly = true;

xFit = linspace(log10(3),log10(70),50);

% create MATLAB table variable puff pressure
T1 = readtable(fullfile(dataPath,'data',spreadsheet1));
allVarNames1 = T1.Properties.VariableNames;

% create MATLAB table variable habituation
T2 = readtable(fullfile(dataPath,'data',spreadsheet2));
allVarNames2 = T2.Properties.VariableNames;

%% get mean R2s

R2s = [];

for vv = 1:length(varNamesToPlot1)
    
    varR2 = [];
    
    for ss = 1:length(subList)

        % find scans for desired subject
        scans = T1(ismember(T1.subjectID,subList{ss}),:);
        scans = scans(ismember(scans.valid,'TRUE'),:);
        if highestOnly
           A = scans(ismember(scans.intendedPSI, 15),:);
           B = scans(ismember(scans.intendedPSI, 30),:);
           C = scans(ismember(scans.intendedPSI, 60),:);
           scans = vertcat(A, B, C);
        end
        ii = find(strcmp(varNamesToPlot1{vv},allVarNames1));

        % across session data
        y = scans.(allVarNames1{ii});
        goodPoints = ~isnan(y);
        x = log10(scans.PSI);
        x = x(goodPoints);
        y = y(goodPoints);
        [x,idxX]=sort(x);
        y = y(idxX);
        weights = scans.numIpsi;
        mSize = weights*20;

        % get R2
        fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weights);
        rsquare = fitObj.Rsquared.Ordinary;
        if rsquare > 1 || rsquare < 0
            rsquare = nan;
        end
        varR2(end+1) = rsquare;
        
    end
    
    R2s(end+1) = mean(varR2,'omitnan');
    
end

%% get r of test retest for slope and offset

roffset = [];
rslope = [];

for vv = 1:length(varNamesToPlot1)
    
    oX = [];
    oY = [];
    pX = [];
    pY = [];
    
    for ss = 1:length(subList)

        % find scans for desired subject
        scans = T1(ismember(T1.subjectID,subList{ss}),:);
        scans = scans(ismember(scans.valid,'TRUE'),:);

        % separate scans into a table for each of the sessions
        dates = unique(scans.scanDate);
        if highestOnly
           A = scans(ismember(scans.intendedPSI, 15),:);
           B = scans(ismember(scans.intendedPSI, 30),:);
           C = scans(ismember(scans.intendedPSI, 60),:);
           scans = vertcat(A, B, C);
        end
        sessOne = scans(ismember(scans.scanDate,dates(1,1)),:);
        sessTwo = scans(ismember(scans.scanDate,dates(2,1)),:);
        ii = find(strcmp(varNamesToPlot1{vv},allVarNames1));

        % session one data
        y = sessOne.(allVarNames1{ii});
        goodPoints = ~isnan(y);
        x = log10(sessOne.PSI);
        x = x(goodPoints);
        y = y(goodPoints);
        [x,idxX]=sort(x);
        y = y(idxX);
        weights = sessOne.numIpsi;
        mSize = weights*20;
        fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weights);
        oX(end+1) = fitObj.Coefficients.Estimate(1);
        pX(end+1) = fitObj.Coefficients.Estimate(2);

        % session two data
        y = sessTwo.(allVarNames1{ii});
        goodPoints = ~isnan(y);
        x = log10(sessTwo.PSI);
        x = x(goodPoints);
        y = y(goodPoints);
        [x,idxX]=sort(x);
        y = y(idxX);
        weights = sessTwo.numIpsi;
        mSize = weights*20;
        fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weights);
        oY(end+1) = fitObj.Coefficients.Estimate(1);
        pY(end+1) = fitObj.Coefficients.Estimate(2);
    end
    
    co = corrcoef(oX,oY);
    roffset(end+1) = co(1,2);
    co = corrcoef(pX,pY);
    rslope(end+1) = co(1,2);
    
end

%% get R2 of the habituation effect across trials

r2Hab = [];

for vv = 1:length(varNamesToPlot2)
    
    varR2 = [];
    
    for ss = 1:length(subList)
        
        acqMeans = NaN(1,25);
        psi = [30 3.75 15 60 60 15 3.75 7.5 30 7.5 7.5 3.75 60 30 15 15 30 60 7.5 15 7.5 60 3.75 3.75 30];
        ii = find(strcmp(varNamesToPlot2{vv},allVarNames2));

        % find scans for desired subject
        scans = T2(ismember(T2.subjectID,subList{ss}),:);
        scans = scans(ismember(scans.valid,'TRUE'),:);

        % get mean values for each scan
        for zz = 1:25
            temp = scans(ismember(scans.scanNumber, zz+1),:);
            if ~isempty(temp) 
               acqMeans(zz) = mean(temp.(varNamesToPlot2{vv}), 'omitnan');
             end
        end

        fitObj = fitlm(log10(psi),acqMeans,'RobustOpts', 'on');
        modelY = fitObj.Fitted;
        
        % calculate residuals as a function of trial number
        resByTrial = NaN(25, 8);
        resMeansByTrial = NaN(1,8);
        for zz = 1:8
            temp = scans(ismember(scans.stimIndex, zz),:);
            if ~isempty(temp)
                for yy = 1:25
                    tt = temp(ismember(temp.scanNumber, yy+1),:);
                    if isempty(tt)
                        residual = NaN;
                    elseif length(tt.(varNamesToPlot2{vv})) == 1
                        residual = tt.(varNamesToPlot2{vv})(1) - modelY(yy);
                    else
                        res1 = tt.(varNamesToPlot2{vv})(1) - modelY(yy);
                        res2 = tt.(varNamesToPlot2{vv})(1) - modelY(yy);
                        residual = mean([res1 res2]);
                    end
                    resByTrial(yy,zz) = residual;
                end
            end
            resMeansByTrial(1,zz) = mean(resByTrial(:,zz), 'omitnan');
        end
        
        % get R2
        fitObj = fitlm((1:8),resMeansByTrial,'RobustOpts', 'on');
        rsquare = fitObj.Rsquared.Ordinary;
        if rsquare > 1 || rsquare < 0
            rsquare = nan;
        end
        varR2(end+1) = rsquare;

    end
        
    r2Hab(end+1) = mean(varR2,'omitnan');
    
end

%% get r of the test retetest of habituation across trials

rHab = [];

for vv = 1:length(varNamesToPlot2)
    
    subjectTrialMeans = [];
    sessOneSlopes = [];
    sessTwoSlopes = [];
    
    for ss = 1:length(subList)
        
        psi = [30 3.75 15 60 60 15 3.75 7.5 30 7.5 7.5 3.75 60 30 15 15 30 60 7.5 15 7.5 60 3.75 3.75 30];
        ii = find(strcmp(varNamesToPlot2{vv},allVarNames2));

        % find scans for desired subject
        scans = T2(ismember(T2.subjectID,subList{ss}),:);
        scans = scans(ismember(scans.valid,'TRUE'),:);
        dates = unique(scans.scanDate);
        sessOne = scans(ismember(scans.scanDate,dates(1,1)),:);
        sessTwo = scans(ismember(scans.scanDate,dates(2,1)),:);

        % calculate residuals as a function of trial number session 1
        acqMeans = NaN(1,25);
        for zz = 1:25
            temp = sessOne(ismember(sessOne.scanNumber, zz+1),:);
            if ~isempty(temp) 
               acqMeans(zz) = mean(temp.(varNamesToPlot2{vv}), 'omitnan');
            end
        end

        fitObj = fitlm(log10(psi),acqMeans,'RobustOpts', 'on');
        modelY = fitObj.Fitted;
        
        resByTrial = NaN(25, 8);
        resMeansByTrial = NaN(1,8);
        for zz = 1:8
            temp = sessOne(ismember(sessOne.stimIndex, zz),:);
            if ~isempty(temp)
                for yy = 1:25
                    tt = temp(ismember(temp.scanNumber, yy+1),:);
                    if isempty(tt)
                        residual = NaN;
                    elseif length(tt.(varNamesToPlot2{vv})) == 1
                        residual = tt.(varNamesToPlot2{vv})(1) - modelY(yy);
                    else
                        res1 = tt.(varNamesToPlot2{vv})(1) - modelY(yy);
                        res2 = tt.(varNamesToPlot2{vv})(1) - modelY(yy);
                        residual = mean([res1 res2]);
                    end
                    resByTrial(yy,zz) = residual;
                end
            end
            resMeansByTrial(1,zz) = mean(resByTrial(:,zz), 'omitnan');
        end
        
        % get session one slope
        fitObj = fitlm((1:8),resMeansByTrial,'RobustOpts', 'on');
        sessOneSlope = fitObj.Coefficients.Estimate(2);
        sessOneSlopes(end+1) = sessOneSlope;
        
        % calculate residuals as a function of trial number session 2
        acqMeans = NaN(1,25);
        for zz = 1:25
            temp = sessTwo(ismember(sessTwo.scanNumber, zz+1),:);
            if ~isempty(temp) 
               acqMeans(zz) = mean(temp.(varNamesToPlot2{vv}), 'omitnan');
            end
        end

        fitObj = fitlm(log10(psi),acqMeans,'RobustOpts', 'on');
        modelY = fitObj.Fitted;
        
        resByTrial = NaN(25, 8);
        resMeansByTrial = NaN(1,8);
        for zz = 1:8
            temp = sessTwo(ismember(sessTwo.stimIndex, zz),:);
            if ~isempty(temp)
                for yy = 1:25
                    tt = temp(ismember(temp.scanNumber, yy+1),:);
                    if isempty(tt)
                        residual = NaN;
                    elseif length(tt.(varNamesToPlot2{vv})) == 1
                        residual = tt.(varNamesToPlot2{vv})(1) - modelY(yy);
                    else
                        res1 = tt.(varNamesToPlot2{vv})(1) - modelY(yy);
                        res2 = tt.(varNamesToPlot2{vv})(1) - modelY(yy);
                        residual = mean([res1 res2]);
                    end
                    resByTrial(yy,zz) = residual;
                end
            end
            resMeansByTrial(1,zz) = mean(resByTrial(:,zz), 'omitnan');
        end
        
        % get session two slope
        fitObj = fitlm((1:8),resMeansByTrial,'RobustOpts', 'on');
        sessTwoSlope = fitObj.Coefficients.Estimate(2);
        sessTwoSlopes(end+1) = sessTwoSlope;
        
        % get mean slope across sessions
        meanSlope = mean([sessOneSlope sessTwoSlope]);
        subjectTrialMeans(end+1) = meanSlope;

    end
    
    co = corrcoef(sessOneSlopes,sessTwoSlopes);
    rHab(end+1) = co(1,2);
    
end

%% write table

col = {'R2 puff pressure model', 'TR r slope', 'TR r offset', 'R2 habituation', 'TR r habituation'};
row = {'AUC', 'Latency', 'Time under 20', 'Time to open', 'Initial velocity', ...
     'Time to close', 'Max closing velocity', 'Max opening velocity', 'Blink rate'};
T = table(R2s', rslope', roffset', r2Hab', rHab', 'VariableNames', col, 'RowNames', row);
writetable(T, fullfile(dataPath,'data','correlationTable.csv'), 'WriteRowNames', 1);