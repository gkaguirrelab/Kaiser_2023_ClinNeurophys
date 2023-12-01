function Kaiser_2023_PsychophysiologyLocalHook
%  Kaiser_2023_PsychophysiologyLocalHook
%
% For use with the ToolboxToolbox.
%
% As part of the setup process, ToolboxToolbox will copy this file to your
% ToolboxToolbox localToolboxHooks directory (minus the "Template" suffix).
% The defalt location for this would be
%   ~/localToolboxHooks/Kaiser_2023_PsychophysiologyLocalHook.m
%
% Each time you run tbUseProject('Kaiser_2023_Psychophysiology'), ToolboxToolbox will
% execute your local copy of this file to do setup for Kaiser_2023_Psychophysiology.
%
% You should edit your local copy with values that are correct for your
% local machine, for example the output directory location.
%


projectName = 'Kaiser_2023_Psychophysiology';

%% Delete any old prefs
if (ispref(projectName))
    rmpref(projectName);
end

% Get user name
if ismac
    [~, userName] = system('whoami');
    userName = strtrim(userName);
elseif isunix
    userName = getenv('USER');
elseif ispc
    userName = getenv('username');
else
    disp('What are you using?')
end

% Get the DropBox path
if ismac
    dbJsonConfigFile = '~/.dropbox/info.json';
    fid = fopen(dbJsonConfigFile);
    raw = fread(fid,inf);
    str = char(raw');
    fclose(fid);
    val = jsondecode(str);
    dropboxBaseDir = val.business.path;
else
    error('Need to set up DropBox path finding for non-Mac machine')
end

% Path to data and analysis directories
switch userName
    case 'aguirre'
        plotSaveDir = fullfile(dropboxBaseDir,'_Papers/xCompleted/2023/Kaiser_2023_psychometricBlink/Figures/matlabFigures');
        analysisDir = fullfile(dropboxBaseDir,'BLNK_analysis','noise_cancellation');
        dataDir = fullfile(dropboxBaseDir,'BLNK_data','noise_cancellation');
    otherwise
        plotSaveDir = fullfile(dropboxBaseDir,'_Papers/xCompleted/2023/Kaiser_2023_psychometricBlink/Figures/matlabFigures');
        analysisDir = fullfile(dropboxBaseDir,'BLNK_analysis','expt01_summer2023');
        dataDir = fullfile(dropboxBaseDir,'BLNK_data','expt01_summer2023');
end

% Set the prefs
setpref(projectName,'plotSaveDir',plotSaveDir);
setpref(projectName,'analysisDir',analysisDir);
setpref(projectName,'dataDir',dataDir);


