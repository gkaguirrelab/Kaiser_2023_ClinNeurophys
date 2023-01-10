%% preRegAnalysis
% Calculates the deviation for PSI.

%% set up parameters

% load file path
dataPath = fileparts(fileparts(mfilename('fullpath')));
spreadsheet1 ='UPenn Ipsi Summary_25ms_02062022.csv';

% create MATLAB table
T = readtable(fullfile(dataPath,'data',spreadsheet1));
allVarNames = T.Properties.VariableNames;

%% calculate deviation

% find valid scans
scans = T(ismember(T.valid,'TRUE'),:);
scans = scans(ismember(scans.numIpsi,(3:8)),:);

% get PSI values
intendedPSI = scans.intendedPSI;
actualPSI = scans.PSI;

% calculate deviation
deviation = actualPSI - intendedPSI;
absDev = abs(deviation);
meanDev = mean(absDev)