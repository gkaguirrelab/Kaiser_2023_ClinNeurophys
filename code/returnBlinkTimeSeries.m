function [blinkVector,temporalSupport, nTrials, blinkVectorRaw, trialIndices] = returnBlinkTimeSeries( subjectID,targetPSI,sessionID,ipsiOrContra,minValidIpsiBlinksPerAcq,minValidAcq,nSamplesBeforeStim,nSamplesAfterStim )
% Loads I-Files and conducts an analysis of time series data
%
% Syntax:
%   [blinkVector,temporalSupport] = returnBlinkTimeSeries( subjectID, targetPSI )
%
% Description:
%   Briana: Describe how the iFiles are structured and how we load them
%
% Inputs:
%   subjectID             - Scalar. 5 digit integer that identifies the
%                           subject.
%   targetPSI             - Scalar. One of the valid PSI targets:
%                               {3.5, 7.5, 15, 30, 60}
%                           If set empty, all PSI levels are used
%   sessionID             - Scalar. Valid values are:
%                               {1, 2}
%                           If set empty, both sessions are used
%
% Outputs:
%   blinkVector           - Vector
%   temporalSupport       - Vector
%   nTrials               - Scalar
%   blinkVectorRaw        - nTrials x time matrix of raw responses
%
% Examples:
%{
    targetPSI = 7.5;
    subjectID = 15513;
    [blinkVector,temporalSupport] = returnBlinkTimeSeries( subjectID, targetPSI, 1, 'ipsi' );
    figure
    plot(temporalSupport,blinkVector,'-r')
%}


arguments
    subjectID (1,1) {mustBeNumeric}
    targetPSI = [];
    sessionID = [];
    ipsiOrContra = 'ipsi';
    minValidIpsiBlinksPerAcq (1,1) {mustBeNumeric} = 0;
    minValidAcq = 0;
    nSamplesBeforeStim = 10;
    nSamplesAfterStim = 150;
end

% Set a counter to return
nTrials = 0;

% Initialize vectors for return
blinkVectorRaw = [];
trialIndices = [];

% Define the location of the i-files
dataDirPath = fileparts(fileparts(mfilename('fullpath')));

% Define the location of the summary spreadsheet that we will use to
% determine if a given trial is valid
spreadsheet ='UPENN Summary with IPSI Responses_02072022_SquintCheck.csv';

% Turn off a warning during readtable
warnState = warning();
warning('off','MATLAB:table:ModifiedAndSavedVarnames');

% Read the table
T = readtable(fullfile(dataDirPath,'data',spreadsheet));

% Restore the warning state
warning(warnState);

% Find the scans for this subject
scanTable = T(ismember(T.subjectID,subjectID),:);

% Store the scan dates
scanDates = unique(scanTable.scanDate);

% Now cull the table to remove invalid scans
scanTable = scanTable(ismember(scanTable.notSquint,'TRUE'),:);
scanTable = scanTable(scanTable.numIpsi>=minValidIpsiBlinksPerAcq,:);

% If we have a targetPSI, filter the table to include just those
if ~isempty(targetPSI)
    scanTable = scanTable(scanTable.intendedPSI==targetPSI,:);
end

% Filter out "invalid" blinks (per BlinkCNS processing) if the setting
% minValidIpsiBlinksPerAcq is greater than zero
if minValidIpsiBlinksPerAcq > 0
    scanTable = scanTable(ismember(scanTable.valid,'TRUE'),:);
end

% Handle working with session 1, 2, or both
if nargin>2 && ~isempty(scanTable)
    dateMatchIdx = scanTable.scanDate == scanDates(sessionID);
    scanTable = scanTable(dateMatchIdx,:);
end

% Check if we have an empty scanTable
if isempty(scanTable)
    blinkVector = nan(1,nSamplesBeforeStim+nSamplesAfterStim+1);
    temporalSupport = nan(1,nSamplesBeforeStim+nSamplesAfterStim+1);
    return
end

% Check if we have too few acquisitions
if size(scanTable,1) < minValidAcq
    blinkVector = nan(1,nSamplesBeforeStim+nSamplesAfterStim+1);
    temporalSupport = nan(1,nSamplesBeforeStim+nSamplesAfterStim+1);
    return
end

% Turn off a warning during readtable
warnState = warning();
warning('off','MATLAB:table:ModifiedAndSavedVarnames');

% Loop over the acquisitions
for ii = 1:size(scanTable,1)

    % In principle, the iFile name should be fully defined by the
    % information in the data table provided by Andy from BlinkTBI. In
    % practice, the scanID and scanNumber in the filenames do not agree
    % with that in the table. It seems that the scanNumber value from the
    % table can be ignored, and the scanID used as the unique identifier of
    % the file. So, we find the matching scanID file, and check to make
    % sure it is the only match
    iFileName = ['l-file_' num2str(scanTable.subjectID(ii)) '_' num2str(scanTable.scanID(ii)) '_' '*' '.csv'];
    fullFilePath = fullfile(dataDirPath,'data','iFiles',num2str(subjectID),iFileName);
    fileList = dir(fullFilePath);
    if length(fileList)>1
        error('Too many files');
    else
        fullFilePath = fullfile(fileList.folder,fileList.name);
    end

    % Load the iFile into a table
    T = readtable(fullFilePath);

    % Get the deltaT for this measure
    timeVarName = T.Properties.VariableNames{1};
    deltaT(ii) = mean(diff(T.(timeVarName)));

    % find stimulus arrivals
    stimVarName = T.Properties.VariableNames{2};
    rights = find(strcmp('MC-OD',T.(stimVarName)(:,1)));
    lefts = find(strcmp('MC-OS',T.(stimVarName)(:,1)));
    all = sort(cat(1,rights,lefts));

    % get times
    starts = all - nSamplesBeforeStim;
    ends = all + nSamplesAfterStim;

    % Issue a warning if there are fewer than 8 trials
    if length(starts)<8
        fprintf(['Only %d trials: ' fullfile(num2str(subjectID),iFileName) '\n'],length(starts));
    end

    % get means across trials
    pos = nan(length(starts),nSamplesBeforeStim+nSamplesAfterStim+1);
    for jj = 1:length(starts)

        % Handle the edge case of the time-series starting after the
        % desired "numBefore" window
        offset = max([1 -starts(jj)+2]);

        % Handle the laterality of the stimulus and the choice of returning
        % the ipsi or contra response
        if ismember(all(jj),rights)
            switch ipsiOrContra
                case 'ipsi'
                    columnIdx = 3;
                case 'contra'
                    columnIdx = 4;
            end
        else
            switch ipsiOrContra
                case 'ipsi'
                    columnIdx = 4;
                case 'contra'
                    columnIdx = 3;
            end
        end

        % Get this timeseries
        temp = table2array(T(max([1 starts(jj)]):ends(jj),columnIdx));

        % Add it to the matrix
        pos(jj,offset:end) = temp';

        % Store the trial index
        trialIndices = [trialIndices jj];

        % Increment the trial count
        nTrials = nTrials + 1;

    end

    % center pre-stimulus around zero
    posAvg = nanmean(pos);
    posAvgPreStim = mean(posAvg(1:nSamplesBeforeStim));
    posAvg = posAvg - posAvgPreStim;

    % Store the response
    respByAcq(ii,:) = posAvg;

    % Create a concatenated, raw vector of responses
    blinkVectorRaw = [blinkVectorRaw; pos-mean(pos(:,1:nSamplesBeforeStim),2)];

end

% Restore the warning state
warning(warnState);

% Get the mean deltaT and assemble the temporal support
deltaT = round(mean(deltaT),4);
temporalSupport = -nSamplesBeforeStim*deltaT:deltaT:nSamplesAfterStim*deltaT;

% Return the blinkVector
blinkVector = mean(respByAcq,1);

end
