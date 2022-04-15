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
     'closeTimeI', 'maxClosingVelocityI', 'maxOpeningVelocityI', 'excursionI', 'closuresI'};

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

%% Figure 1c left: example combined session blink feature across puff intensities

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

%% Figure 1c right: puff pressure R2 histogram
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

    % session one means
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

    % get rsquare
    fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weights);
    rsquare = fitObj.Rsquared.Ordinary;
    if rsquare > 1 || rsquare < 0
        rsquare = nan;
    end
    rsquares(end+1) = rsquare;

end

histogram(rsquares,10);
title(['R2 for max closing velocity fit across puff pressure'], 'FontSize', 14)
xlabel('R2 value', 'FontSize', 14)
ylabel('Number of subjects', 'FontSize', 14);

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

% session one means
sm = NaN(1,5);
aa = sessOne(ismember(sessOne.intendedPSI, 3.5),:);
bb = sessOne(ismember(sessOne.intendedPSI, 7.5),:);
cc = sessOne(ismember(sessOne.intendedPSI, 15),:);
dd = sessOne(ismember(sessOne.intendedPSI, 30),:);
ee = sessOne(ismember(sessOne.intendedPSI, 60),:);
sm(1) = mean(aa.(allVarNames{ii}), 'omitnan');
sm(2) = mean(bb.(allVarNames{ii}), 'omitnan');
sm(3) = mean(cc.(allVarNames{ii}), 'omitnan');
sm(4) = mean(dd.(allVarNames{ii}), 'omitnan');
sm(5) = mean(ee.(allVarNames{ii}), 'omitnan');

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