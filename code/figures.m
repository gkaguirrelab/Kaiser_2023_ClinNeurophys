%% figures
%
%       Scan PSI index (out of 26 scans):
%          3.5 PSI: [3 8 13 24 25]
%          7.5 PSI: [9 11 12 20 22]
%          15 PSI: [4 7 16 17 21]
%          30 PSI: [1 2 10 15 18 26]
%          60 PSI: [5 6 14 19 23]
%

%% set up parameters

% load file path
dataPath = fileparts(fileparts(mfilename('fullpath')));
spreadsheet ='UPenn Ipsi Summary_25ms_02062022.csv';

% choose subject and parameters
highestOnly = true;

varNamesToPlot = {'aucI', 'latencyI', 'timeUnderI', 'openTimeI', 'initVelocityI', ...
     'closeTimeI', 'maxClosingVelocityI', 'maxOpeningVelocityI', 'blinkRateI'};
 
varNamesToPlot2 = {'auc', 'latency', 'timeUnder20', 'openTime', 'initialVelocity', ...
      'closeTime', 'maxClosingVelocity', 'maxOpeningVelocity', 'blinkRate'};

if highestOnly
    subList = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591, ...
    14590, 14589, 14588, 14587, 14586};
else
    subList = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591, ...
    14587, 14586};
end

xFit = linspace(log10(3),log10(70),50);

% create MATLAB table variable
T = readtable(fullfile(dataPath,'data',spreadsheet));
allVarNames = T.Properties.VariableNames;

% create MATLAB table variable for habituation data
spreadsheet ='Upenn_Ipsilateral Afiles_clean_full.csv';
T2 = readtable(fullfile(dataPath,'data',spreadsheet));
allVarNames2 = T2.Properties.VariableNames;

%% Figure 1b: plot puff intensity as a function of scan number to show counterbalanced sequence

pressures = [30 30 3.75 15 60 60 15 3.75 7.5 30 7.5 7.5 3.75 60 30 15 15 30 60 7.5 15 7.5 60 3.75 3.75 30];

figure();
x = (1:26);
y = log10(pressures);

scatter(x,y,30,'filled');
hold on;
stairs(x,y);
title(['Counterbalanced pressure sequence'], 'FontSize', 14)
ylabel(['Puff pressure [log psi]'], 'FontSize', 14)
xlabel('Acquisition number', 'FontSize', 14)

%% Figure 2a: example combined session blink feature across puff intensities

subject = 14591;
feature = 'maxClosingVelocityI';

figure();

% find scans for desired subject
scans = T(ismember(T.subjectID,subject),:);
scans = scans(ismember(scans.valid,'TRUE'),:);
dates = unique(scans.scanDate);
sessOne = scans(ismember(scans.scanDate,dates(1,1)),:);
sessTwo = scans(ismember(scans.scanDate,dates(2,1)),:);
ii = find(strcmp(feature,allVarNames));

% data
y = scans.(allVarNames{ii});
goodPoints = ~isnan(y);
x = log10(scans.PSI);
x = x(goodPoints);
y = y(goodPoints);
[x,idxX]=sort(x);
y = y(idxX);
weights = scans.numIpsi;
mSize = weights*20;

% means
sm = NaN(1,5);
aa = scans(ismember(scans.intendedPSI, 3.5),:);
bb = scans(ismember(scans.intendedPSI, 7.5),:);
cc = scans(ismember(scans.intendedPSI, 15),:);
dd = scans(ismember(scans.intendedPSI, 30),:);
ee = scans(ismember(scans.intendedPSI, 60),:);
sm(1) = mean(aa.(allVarNames{ii}), 'omitnan');
sm(2) = mean(bb.(allVarNames{ii}), 'omitnan');
sm(3) = mean(cc.(allVarNames{ii}), 'omitnan');
sm(4) = mean(dd.(allVarNames{ii}), 'omitnan');
sm(5) = mean(ee.(allVarNames{ii}), 'omitnan');

% make plot
scatter(x,y,mSize);
fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weights);
hold on
plot(x,fitObj.Fitted,'-r');
scatter(log10([3.5 7.5 15 30 60]),sm,300,'r');
rsquare = fitObj.Rsquared.Ordinary;
if rsquare > 1 || rsquare < 0
    rsquare = nan;
end
title(['Subject ' num2str(subject)], 'FontSize', 14)
ylabel(feature, 'FontSize', 14)
xlabel('puff pressure [log psi]', 'FontSize', 14);
xlim([0.3 2]);

%% Figure 2b: puff pressure R2 histogram
figure();
rsquares = [];
    
for ss = 1:length(subList)

    % find scans for desired subject
    scans = T(ismember(T.subjectID,subList{ss}),:);
    scans = scans(ismember(scans.valid,'TRUE'),:);

    % separate scans into a table for each of the sessions
    dates = unique(scans.scanDate);
    if highestOnly
       A = scans(ismember(scans.intendedPSI, 15),:);
       B = scans(ismember(scans.intendedPSI, 30),:);
       C = scans(ismember(scans.intendedPSI, 60),:);
       scans = vertcat(A, B, C);
    end
    ii = find(strcmp(feature,allVarNames));

    % session one data
    y = scans.(allVarNames{ii});
    goodPoints = ~isnan(y);
    x = log10(scans.PSI);
    x = x(goodPoints);
    y = y(goodPoints);
    [x,idxX]=sort(x);
    y = y(idxX);
    weights = scans.numIpsi;
    mSize = weights*20;

    % get rsquare
    fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weights);
    rsquare = fitObj.Rsquared.Ordinary;
    if rsquare > 1 || rsquare < 0
        rsquare = nan;
    end
    rsquares(end+1) = rsquare;

end

% make histogram
rsquares = sort(rsquares);
inc = 0.1;
count = 0;
max = 0;
while ~isempty(rsquares)
    if rsquares(1) < inc
        count = count + 1;
        plot(inc, count, 'k.', 'LineWidth', 2, 'MarkerSize', 100);
        hold on
        if count > max
            max = count;
        end
        rsquares(1) = [];
    else
        inc = inc + 0.1;
        count = 0;
    end
end
title(['R2 for max closing velocity fit across puff pressure'], 'FontSize', 14)
xlabel('R2 value', 'FontSize', 14)
ylabel('Number of subjects', 'FontSize', 14);
ylim([0 max]);
xlim([0 1]);
xticks(0:0.1:1);
hold off

%% Figure 2a: subject example session 1 and 2

% session one data
y = sessOne.(allVarNames{ii});
goodPoints = ~isnan(y);
x = log10(sessOne.PSI);
x = x(goodPoints);
y = y(goodPoints);
[x,idxX]=sort(x);
y = y(idxX);
weights = sessOne.numIpsi;
mSize = weights*20;

% make plot
subplot(1,2,1);
scatter(x,y,mSize);
%         hold on
%         scatter(log10([3.5 7.5 15 30 60]),sm,300);
fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weights);
hold on
plot(x,fitObj.Fitted,'-r');
xlim(log10([2 100]));
rsquare = fitObj.Rsquared.Ordinary;
if rsquare > 1 || rsquare < 0
    rsquare = nan;
end
title(['Subject ' num2str(subject)], 'FontSize', 14)
ylabel(['Session one ' feature], 'FontSize', 14)
xlabel('puff pressure [log psi]', 'FontSize', 14)

% session two data
y = sessTwo.(allVarNames{ii});
goodPoints = ~isnan(y);
x = log10(sessTwo.PSI);
x = x(goodPoints);
y = y(goodPoints);
[x,idxX]=sort(x);
y = y(idxX);
weights = sessTwo.numIpsi;
mSize = weights*20;

% make plot
subplot(1,2,2);
scatter(x,y,mSize);
fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weights);
hold on
plot(x,fitObj.Fitted,'-r')
rsquare = fitObj.Rsquared.Ordinary;
if rsquare > 1 || rsquare < 0
    rsquare = nan;
end
ylabel(['Session two ' feature], 'FontSize', 14)

%% Figure 3: slope and offset test retest example
feature = 'maxClosingVelocityI';

pX = [];
pY = [];
oX = [];
oY = [];

for ss = 1:length(subList)

    % find scans for desired subject
    scans = T(ismember(T.subjectID,subList{ss}),:);
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
    ii = find(strcmp(feature,allVarNames));

    % session one data
    y = sessOne.(allVarNames{ii});
    goodPoints = ~isnan(y);
    x = log10(sessOne.PSI);
    x = x(goodPoints);
    y = y(goodPoints);
    [x,idxX]=sort(x);
    y = y(idxX);
    weights = sessOne.numIpsi;

    % get session one slope and offset values
    fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weights);
    pX(end+1) = fitObj.Coefficients.Estimate(2);
    oX(end+1) = fitObj.Coefficients.Estimate(1);
    rsquare = fitObj.Rsquared.Ordinary;

    % session two data
    y = sessTwo.(allVarNames{ii});
    goodPoints = ~isnan(y);
    x = log10(sessTwo.PSI);
    x = x(goodPoints);
    y = y(goodPoints);
    [x,idxX]=sort(x);
    y = y(idxX);
    weights = sessTwo.numIpsi;
    mSize = weights*20;

    % get session two slope and offset values
    fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weights);
    pY(end+1) = fitObj.Coefficients.Estimate(2);
    oY(end+1) = fitObj.Coefficients.Estimate(1);
end

% plot parameter test retest values across subjects
figure();
pl = subplot(1,1,1);
plot(pX, pY, 'ob', 'MarkerSize', 10);
fitObj = fitlm(pX,pY,'RobustOpts', 'on');
hold on
plot(fitObj);
plot((-2:5),(-2:5),'k');
pl.Box = 'off';
rsquare = fitObj.Rsquared.Ordinary;
if rsquare > 1 || rsquare < 0
    rsquare = nan;
end
title([feature ' slope by session - ' sprintf(' R^2=%2.2f',rsquare)], 'FontSize', 16)
xlabel(['Slope'], 'FontSize', 16)
ylabel(['Slope'], 'FontSize', 16)
ylim(xlim);
axis(pl, 'square');

% plot offset test retest values across subjects
figure();
pl = subplot(1,1,1);
plot(oX, oY, 'ob', 'MarkerSize', 10);
fitObj = fitlm(oX,oY,'RobustOpts', 'on');
hold on
plot(fitObj);
plot((-2:14),(-2:14),'k');
pl.Box = 'off';
rsquare = fitObj.Rsquared.Ordinary;
if rsquare > 1 || rsquare < 0
    rsquare = nan;
end
title([feature ' offset by session - ' sprintf(' R^2=%2.2f',rsquare)], 'FontSize', 16)
xlabel(['Offset'], 'FontSize', 16)
ylabel(['Offset'], 'FontSize', 16)
ylim(xlim);
axis(pl, 'square');

%% Figure 5a: example residuals for a subject across trials

subject = 14591;
feature = 'maxClosingVelocity';

figure();

acqMeans = NaN(1,25);
resMeansByAcq = NaN(1,25);
psi = [30 3.75 15 60 60 15 3.75 7.5 30 7.5 7.5 3.75 60 30 15 15 30 60 7.5 15 7.5 60 3.75 3.75 30];
ii = find(strcmp(feature,allVarNames2));

% find scans for desired subject
scans = T2(ismember(T2.subjectID,subject),:);
scans = scans(ismember(scans.valid,'TRUE'),:);

% get mean values for each scan
for zz = 1:25
    temp = scans(ismember(scans.scanNumber, zz+1),:);
    if ~isempty(temp) 
       acqMeans(zz) = mean(temp.(feature), 'omitnan');
     end
end

fitObj = fitlm(log10(psi),acqMeans,'RobustOpts', 'on');
modelY = fitObj.Fitted;

% calculate residuals as a function of acquisition number
for zz = 1:25            
    temp = scans(ismember(scans.scanNumber, zz+1),:);            
    if ~isempty(temp)
        residuals = temp.(feature)' - modelY(zz);
        resMeansByAcq(zz) = mean(residuals,'omitnan');
    end            
end

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
            elseif length(tt.(feature)) == 1
                residual = tt.(feature)(1) - modelY(yy);
            else
                res1 = tt.(feature)(1) - modelY(yy);
                res2 = tt.(feature)(1) - modelY(yy);
                residual = mean([res1 res2]);
            end
            resByTrial(yy,zz) = residual;
        end
    end
    resMeansByTrial(1,zz) = mean(resByTrial(:,zz), 'omitnan');
end

% plot mean residual as a function of trial number
scatter((1:8),resMeansByTrial);
fitObj = fitlm((1:8),resMeansByTrial,'RobustOpts', 'on');
hold on
plot((1:8),fitObj.Fitted,'-r');
rsquare = fitObj.Rsquared.Ordinary;
if rsquare > 1 || rsquare < 0
    rsquare = nan;
end
title({['Subject ' num2str(subject)],sprintf(' R^2=%2.2f',rsquare)}, 'FontSize', 14)
ylabel([feature ' residual'], 'FontSize', 14)
xlabel('Trial number', 'FontSize', 14)

%% Figure 5b: histogram of residual R2 values across subjects

rsquares = [];
feature = 'maxClosingVelocity';
    
for ss = 1:length(subList)

    acqMeans = NaN(1,25);
    resMeansByAcq = NaN(1,25);
    psi = [30 3.75 15 60 60 15 3.75 7.5 30 7.5 7.5 3.75 60 30 15 15 30 60 7.5 15 7.5 60 3.75 3.75 30];
    ii = find(strcmp(feature,allVarNames2));

    % find scans for desired subject
    scans = T2(ismember(T2.subjectID,subList{ss}),:);
    scans = scans(ismember(scans.valid,'TRUE'),:);

    % get mean values for each scan
    for zz = 1:25
        temp = scans(ismember(scans.scanNumber, zz+1),:);
        if ~isempty(temp) 
           acqMeans(zz) = mean(temp.(feature), 'omitnan');
         end
    end

    fitObj = fitlm(log10(psi),acqMeans,'RobustOpts', 'on');
    modelY = fitObj.Fitted;

    % calculate residuals as a function of acquisition number
    for zz = 1:25            
        temp = scans(ismember(scans.scanNumber, zz+1),:);            
        if ~isempty(temp)
            residuals = temp.(feature)' - modelY(zz);
            resMeansByAcq(zz) = mean(residuals,'omitnan');
        end            
    end

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
                elseif length(tt.(feature)) == 1
                    residual = tt.(feature)(1) - modelY(yy);
                else
                    res1 = tt.(feature)(1) - modelY(yy);
                    res2 = tt.(feature)(1) - modelY(yy);
                    residual = mean([res1 res2]);
                end
                resByTrial(yy,zz) = residual;
            end
        end
        resMeansByTrial(1,zz) = mean(resByTrial(:,zz), 'omitnan');
    end

    % get rsquares
    fitObj = fitlm((1:8),resMeansByTrial,'RobustOpts', 'on');
    rsquare = fitObj.Rsquared.Ordinary;
    if rsquare > 1 || rsquare < 0
        rsquare = nan;
    end
    rsquares(end+1) = rsquare;
end

% make histogram
rsquares = sort(rsquares);
inc = 0.1;
count = 0;
max = 0;
while ~isempty(rsquares)
    if rsquares(1) < inc
        count = count + 1;
        plot(inc, count, 'k.', 'LineWidth', 2, 'MarkerSize', 100);
        hold on
        if count > max
            max = count;
        end
        rsquares(1) = [];
    else
        inc = inc + 0.1;
        count = 0;
    end
end
title(['Max closing velocity residual fit across trials'], 'FontSize', 14)
xlabel('R2 value', 'FontSize', 14)
ylabel('Number of subjects', 'FontSize', 14);
ylim([0 max]);
xlim([0 1]);
xticks(0:0.1:1);
hold off