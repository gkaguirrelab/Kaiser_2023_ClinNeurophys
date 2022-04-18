%% timeSeriesAnalysis
% Loads I-Files and conducts an analysis of time series data

%% set up parameters

% load file path
location = '/Users/brianahaggerty/Documents/MATLAB/projects/blinkCNSAnalysis/data/iFiles/14591/';
dir(location);
ds = spreadsheetDatastore(location);
ds.ReadSize = 'file';


%% plot stimulus arrival time as a function of puff pressure
x = [];
y = [];

% loop through scan files
for ii = 1:52
   % make table from file
   temp = read(ds);
   T = readtable(temp);
   allVarNames = T.Properties.VariableNames;
   
   % find stimulus arrival
   r = T(ismember(T.Stimulus,'MC-OD'),:);
   t = r.Time;
   
   if ismember(ii, [3 8 13 24 25 29 34 38 50 51])
       % 3.5 PSI
       x(end+1) = 3.5;
       y(end+1) = t;       
   elseif ismember(ii, [9 11 12 20 22 35 37 38 46 48])
       % 7.5 PSI
       x(end+1) = 7.5;
       y(end+1) = t;
   elseif ismember(ii, [4 7 16 17 21 30 33 42 43 47])
       % 15 PSI
       x(end+1) = 15;
       y(end+1) = t;
   elseif ismember(ii, [2 10 15 18 26 28 36 41 44 52])
       % 30 PSI
       x(end+1) = 30;
       y(end+1) = t;
   else
       % 60 PSI
       x(end+1) = 60;
       y(end+1) = t;
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
