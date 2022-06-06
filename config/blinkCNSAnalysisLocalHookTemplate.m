function blinkCNSAnalysisLocalHook
%  blinkCNSAnalysisLocalHook
%
% For use with the ToolboxToolbox.
%
% As part of the setup process, ToolboxToolbox will copy this file to your
% ToolboxToolbox localToolboxHooks directory (minus the "Template" suffix).
% The defalt location for this would be
%   ~/localToolboxHooks/blinkCNSAnalysisLocalHook.m
%
% Each time you run tbUseProject('blinkCNSAnalysis'), ToolboxToolbox will
% execute your local copy of this file to do setup for blinkCNSAnalysis.
%
% You should edit your local copy with values that are correct for your
% local machine, for example the output directory location.
%


projectName = 'blinkCNSAnalysis';

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
        plotSaveDir = fullfile(dropboxBaseDir,'_Papers/Kaiser_2022_psychometricBlink/Figures/matlabFigures');
    otherwise
        plotSaveDir = '';
end

setpref(projectName,'plotSaveDir',plotSaveDir);
