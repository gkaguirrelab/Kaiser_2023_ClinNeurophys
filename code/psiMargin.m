%% psiMargin
% This script loads a blink data set into a MATLAB table variable. When
% run, it will calculate the margin of error between the intended and
% actual PSI.
%% set up parameters

% load file path
dataPath = fileparts(fileparts(mfilename('fullpath')));
spreadsheet ='UPenn Ipsi Summary_25ms_02062022.csv';

% choose subject and parameters
highestOnly = true;
if highestOnly
    subList = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591, ...
    14590, 14589, 14588, 14587, 14586};
else
    subList = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591, ...
    14587, 14586};
end

xFit = linspace(log10(3),log10(70),50);
% ylims = {[30 65]};

% create MATLAB table variable
T = readtable(fullfile(dataPath,'data',spreadsheet));
allVarNames = T.Properties.VariableNames;

%% calculate

scans = T(ismember(T.valid,'TRUE'),:);
intended = scans.intendedPSI;
actual = scans.PSI;
difference = abs(intended - actual);
percent = (difference ./ intended)*100;
maxNum = max(percent)

col = {'intended psi', 'actual psi', 'difference', 'percent'};
T = table(intended, actual, difference, percent, 'VariableNames', col);
writetable(T, fullfile(dataPath,'data','psiMargin.csv'), 'WriteRowNames', 1);