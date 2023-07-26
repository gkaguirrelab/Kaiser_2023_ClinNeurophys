function Kaiser_2023_ClinNeurophysLocalHook
%  Kaiser_2023_ClinNeurophysLocalHook
%
% For use with the ToolboxToolbox.
%
% As part of the setup process, ToolboxToolbox will copy this file to your
% ToolboxToolbox localToolboxHooks directory (minus the "Template" suffix).
% The defalt location for this would be
%   ~/localToolboxHooks/Kaiser_2023_ClinNeurophysLocalHook.m
%
% Each time you run tbUseProject('Kaiser_2023_ClinNeurophys'), ToolboxToolbox will
% execute your local copy of this file to do setup for Kaiser_2023_ClinNeurophys.
%
% You should edit your local copy with values that are correct for your
% local machine, for example the output directory location.
%


projectName = 'Kaiser_2023_ClinNeurophys';

%% Delete any old prefs
if (ispref(projectName))
    rmpref(projectName);
end

% Define a save location for plots
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


% DropBox
dropboxBaseDir = ...
    fullfile('/Users', userName, ...
    'Dropbox (Aguirre-Brainard Lab)');

% Path to paper directory
switch userName
    case 'aguirre'
        plotSaveDir = fullfile(dropboxBaseDir,'_Papers/xCompleted/2023/Kaiser_2023_psychometricBlink/Figures/matlabFigures');
    otherwise
        plotSaveDir = '';
end

setpref(projectName,'plotSaveDir',plotSaveDir);
