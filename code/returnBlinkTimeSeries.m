function [blinkVector,temporalSupport] = returnBlinkTimeSeries( subjectID, targetPSI )
% Loads I-Files and conducts an analysis of time series data
%
% Syntax:
%   output = myFunc(input)
%
% Description:
%   Lorem ipsum dolor sit amet, consectetur adipiscing elit. Aenean euismod
%   nulla a tempor scelerisque. Maecenas et lobortis est. Donec et turpis
%   sem. Sed fringilla in metus ut malesuada. Pellentesque nec eros
%   efficitur, pellentesque nisl vel, dapibus felis. Morbi eu gravida enim.
%   Sed sodales ipsum eget finibus dapibus. Fusce sagittis felis id orci
%   egestas, non convallis neque porttitor. Proin ut mi augue. Cras posuere
%   diam at purus dignissim, vel vestibulum tellus ultrices
%
% Inputs:
%   none
%   foo                   - Scalar. Foo foo foo foo foo foo foo foo foo foo
%                           foo foo foo foo foo foo foo foo foo foo foo foo
%                           foo foo foo
%
% Optional key/value pairs:
%   none
%  'bar'                  - Scalar. Bar bar bar bar bar bar bar bar bar bar
%                           bar bar bar bar bar bar bar bar bar bar bar bar
%                           bar bar bar bar bar bar
%
% Outputs:
%   none
%   baz                   - Cell. Baz baz baz baz baz baz baz baz baz baz
%                           baz baz baz baz baz baz baz baz baz baz baz baz
%                           baz baz baz
%
% Examples:
%{
	foo = 1;
    bar = myFunc(foo);
	fprintf('Bar = %d \n',bar);   
%}

minValidIpsiBlinksPerAcq = 3;

% Define the location of the i-files
dataDirPath = fileparts(fileparts(mfilename('fullpath')));

% Define the location of the summary spreadsheet that we will use to
% determine if a given trial is valid
spreadsheet ='UPENN Summary with IPSI Responses_02072022_SquintCheck.csv';
T = readtable(fullfile(dataDirPath,'data',spreadsheet));

% Find valid scans for this subject and PSI
scanTable = T(ismember(T.subjectID,subjectID),:);
scanTable = scanTable(ismember(scanTable.valid,'TRUE'),:);
scanTable = scanTable(ismember(scanTable.notSquint,'TRUE'),:);
scanTable = scanTable(scanTable.numIpsi>=minValidIpsiBlinksPerAcq,:);
scanTable = scanTable(scanTable.intendedPSI==targetPSI,:);

% decide row bounds
numBefore = 10;
numAfter = 150;

% Check if we have an empty scanTable
if isempty(scanTable)
    blinkVector = nan(1,numBefore+numAfter+1);
    return
end

temporalSupport = -numBefore:1:numAfter;

% loop through scan files
scanNumberOffset = 0;

for ii = 1:size(scanTable,1)


    iFileName = ['l-file_' num2str(scanTable.subjectID(ii)) '_' num2str(scanTable.scanID(ii)) '_' num2str(scanTable.scanNumber(ii)) '.csv'];
    fullFilePath = fullfile(dataDirPath,'data','iFiles',num2str(subjectID),iFileName);

    % Load the iFile into a table
    T = readtable(fullFilePath);

    % Get the deltaT for this measure
    deltaT(ii) = mean(diff(T.Time_msec_));

    % find stimulus arrivals
    [rights,col,v] = find(strcmp('MC-OD',T.Stimulus(:,1)));
    [lefts,col,v] = find(strcmp('MC-OS',T.Stimulus(:,1)));
    all = sort(cat(1,rights,lefts));

    % get times
    starts = all - numBefore;
    ends = all + numAfter;

    % get means across trials
    pos = nan(length(starts),numBefore+numAfter+1);
    for jj = 1:length(starts)

        % Handle the edge case of the time-series starting after the
        % desired "numBefore" window
        offset = max([1 -starts(jj)+2]);

        % Handle the laterality of the stimulus and
        if ismember(all(jj),rights)
            columnIdx = 3;
        else
            columnIdx = 4;
        end

        % Get this timeseries
        temp = table2array(T(max([1 starts(jj)]):ends(jj),columnIdx));

        % Add it to the matrix
        pos(jj,offset:end) = temp';

    end

    pos = nanmean(pos);

    % center pre-stimulus around zero and convert to % change
    pre = pos(1:numBefore);
    mm = mean(pre);
    pos = (pos - mm);%./mm;

    % Store the response
    respByAcq(ii,:) = pos;
end

deltaT = round(mean(deltaT),4);
temporalSupport = -numBefore*deltaT:deltaT:numAfter*deltaT;
blinkVector = mean(respByAcq);

end
