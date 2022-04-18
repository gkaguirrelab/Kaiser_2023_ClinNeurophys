%% timeSeriesAnalysis
% Loads I-Files and conducts an analysis of time series data

%% set up parameters

% load file path
subject = 14591;
location = strcat('/Users/brianahaggerty/Documents/MATLAB/projects/blinkCNSAnalysis/data/iFiles/', num2str(subject),'/');
dir(location);
ds = tabularTextDatastore(location,'FileExtensions',[".csv"]);
ds.ReadSize = 'file';


%% plot stimulus arrival time as a function of puff pressure
x = [];
y = [];

% loop through scan files
for ii = 1:52
   % make table from file
   [data,info] = read(ds);
   T = data;
   subStr = eraseBetween(info.Filename,1,102);
   subStr = extractBetween(subStr,"_",".csv");
   scan = str2num(subStr{1});
   allVarNames = T.Properties.VariableNames;
   
   % find stimulus arrival
   start = vertcat(T(ismember(T.Stimulus,'TA-OD'),:), T(ismember(T.Stimulus,'TA-OS'),:));
   arrivals = vertcat(T(ismember(T.Stimulus,'MC-OD'),:), T(ismember(T.Stimulus,'MC-OS'),:));
   tEnd = table2array(arrivals(:,1));
   tStart = table2array(start(:,1));
   
   if ismember(scan, [3 8 13 24 25 29 34 38 50 51])
       % 3.5 PSI
       for jj = 1:length(tEnd)
           x(end+1) = 3.5;
           y(end+1) = tEnd(jj) - tStart(jj);
       end
   elseif ismember(scan, [9 11 12 20 22 35 37 38 46 48])
       % 7.5 PSI       
       for jj = 1:length(tEnd)
           x(end+1) = 7.5;
           y(end+1) = tEnd(jj) - tStart(jj);
       end
   elseif ismember(scan, [4 7 16 17 21 30 33 42 43 47])
       % 15 PSI
       for jj = 1:length(tEnd)
           x(end+1) = 15;
           y(end+1) = tEnd(jj) - tStart(jj);
       end
   elseif ismember(scan, [2 10 15 18 26 28 36 41 44 52])
       % 30 PSI
       for jj = 1:length(tEnd)
           x(end+1) = 30;
           y(end+1) = tEnd(jj) - tStart(jj);
       end
   elseif ismember(scan, [5 6 14 19 23 31 32 40 45 49])
       % 60 PSI
       for jj = 1:length(tEnd)
           x(end+1) = 60;
           y(end+1) = tEnd(jj) - tStart(jj);
       end
   end
end

figure();
x = log10(x);
scatter(x,y);
hold on
fitObj = fitlm(x,y,'RobustOpts', 'on');
plot(x,fitObj.Fitted,'-r')
xlim(log10([2 100]));
title('Stimulus arrival time across puff pressure', 'FontSize', 16);
xlabel(['Puff pressure [log psi]'], 'FontSize', 16);
ylabel(['Time [msec]'], 'FontSize', 16);
