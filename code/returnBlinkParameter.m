function [meanParameter,nAcq,parameterVector] = returnBlinkParameter( parameterName,subjectID,targetPSI,sessionID,minValidIpsiBlinksPerAcq,minValidAcq )
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
    targetPSI = 60;
    subjectID = 15513;
    meanParameter = returnBlinkParameter( 'excursionI', subjectID, targetPSI, 1 );
%}


arguments
    parameterName {isscalar,mustBeMember(parameterName,{...
        'excursionI','excursionC','latencyI','latencyC'...
        'maxClosingVelocityI','maxClosingVelocityC','aucI','aucC'
        })}
    subjectID (1,1) {mustBeNumeric}
    targetPSI = [];
    sessionID = [];
    minValidIpsiBlinksPerAcq (1,1) {mustBeNumeric} = 0;
    minValidAcq = 0;
end

% Set a counter to return
nAcq = 0;

% Initialize vectors for return
    meanParameter = nan;
    parameterVector = nan;

% Define the location of the i-files
dataDirPath = fileparts(fileparts(mfilename('fullpath')));

% Define the location of the summary spreadsheet that we will use to
% determine if a given trial is valid
spreadsheet ='UPENN Summary with IPSI Responses_02072022_SquintCheck.csv';

% Define the sub-directory in which these data live
dataSubdir = 'Kaiser2023_17PatientTrial';

% Turn off a warning during readtable
warnState = warning();
warning('off','MATLAB:table:ModifiedAndSavedVarnames');

% Read the table
T = readtable(fullfile(dataDirPath,'data',dataSubdir,spreadsheet));

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
if ~isempty(sessionID) && ~isempty(scanTable)
    dateMatchIdx = scanTable.scanDate == scanDates(sessionID);
    scanTable = scanTable(dateMatchIdx,:);
end

% Check if we have an empty scanTable
if isempty(scanTable)
    return
end

% Check if we have too few acquisitions
if size(scanTable,1) < minValidAcq
    return
end

% Extract the parameter
nAcq = length(scanTable.(parameterName));
parameterVector = scanTable.(parameterName);
meanParameter = mean(parameterVector);

end
