%% compareParams
% This script loads a blink data set into a MATLAB table variable. When
% run, it will aggregate data for a given subject for the given parameters 
% across sessions. It will create a matrix of slopes for each parameter.
%%

% load file path
dataPath = fileparts(fileparts(mfilename('fullpath')));
spreadsheet ='2_2022.csv';

% choose subject and parameters
subList = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591};
varNamesToPlot = {'aucI', 'latencyI', 'timeUnderI', 'openTimeI', 'initVelocityI', ...
    'closeTimeI', 'maxClosingVelocityI', 'maxOpeningVelocityI', 'excursionI', 'closuresI'};

% create MATLAB table variable
T = readtable(fullfile(dataPath,'data',spreadsheet));
allVarNames = T.Properties.VariableNames;
slopes = zeros(length(subList), length(varNamesToPlot));

% create slopes matrix containing the slope values for each var and subject
for vv = 1:length(varNamesToPlot)

    for ss = 1:length(subList)

        % find scans for desired subject
        scans = T(ismember(T.subjectID,subList{ss}),:);
        scans = scans(ismember(scans.valid,'TRUE'),:);
        ii = find(strcmp(varNamesToPlot{vv},allVarNames));
        weights = scans.numIpsi;

        % subject parameter data
        y = scans.(allVarNames{ii});
        goodPoints = ~isnan(y);
        x = log10(scans.PSI);
        x = x(goodPoints);
        y = y(goodPoints);
        [x,idxX]=sort(x);
        y = y(idxX);

        % get slope
        fitObj = fitlm(x,y,'RobustOpts', 'on', 'Weight', weights);
        slopes(ss, vv) = fitObj.Coefficients.Estimate(2);

    end

end

comp = corr(slopes);
