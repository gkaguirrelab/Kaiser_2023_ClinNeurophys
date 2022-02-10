dataPath = fileparts(fileparts(mfilename('fullpath')));
spreadsheet ='UPenn Ipsi Summary_25ms_02062022.csv';

T = readtable(fullfile(dataPath,'data',spreadsheet));
allVarNames = T.Properties.VariableNames;

scans = T(ismember(T.subjectID,15514),:);

% separate scans into a table for each of the sessions
dates = unique(scans.scanDate);
sessOne = scans(ismember(scans.scanDate,dates(1,1)),:);
sessTwo = scans(ismember(scans.scanDate,dates(2,1)),:);

aOne = sessOne(ismember(sessOne.intendedPSI,3.5),:);
aTwo = sessTwo(ismember(sessTwo.intendedPSI,3.5),:);

bOne = sessOne(ismember(sessOne.intendedPSI,7.5),:);
bTwo = sessTwo(ismember(sessTwo.intendedPSI,7.5),:);

cOne = sessOne(ismember(sessOne.intendedPSI,15),:);
cTwo = sessTwo(ismember(sessTwo.intendedPSI,15),:);

dOne = sessOne(ismember(sessOne.intendedPSI,30),:);
dTwo = sessTwo(ismember(sessTwo.intendedPSI,30),:);

eOne = sessOne(ismember(sessOne.intendedPSI,60),:);
eTwo = sessTwo(ismember(sessTwo.intendedPSI,60),:);

numValid = {height(aOne(ismember(aOne.valid,'TRUE'),:)), height(aTwo(ismember(aTwo.valid,'TRUE'),:)),...
    height(bOne(ismember(bOne.valid,'TRUE'),:)), height(bTwo(ismember(bTwo.valid,'TRUE'),:)),...
    height(cOne(ismember(cOne.valid,'TRUE'),:)), height(cTwo(ismember(cTwo.valid,'TRUE'),:)),...
    height(dOne(ismember(dOne.valid,'TRUE'),:)), height(dTwo(ismember(dTwo.valid,'TRUE'),:)),...
    height(eOne(ismember(eOne.valid,'TRUE'),:)), height(eTwo(ismember(eTwo.valid,'TRUE'),:))}