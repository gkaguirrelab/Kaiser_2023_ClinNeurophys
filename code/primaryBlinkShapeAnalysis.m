%% primaryBlinkShapeAnalysis
% This routine conducts an analysis upon the time-series measurements of
% lid position evoked by air puffs, and as measured using the BlinkTBI
% device. Each of 17 subjects were studied in two separate sessions. Each
% session had multiple acquisitions. Each acquisition used one of 5 air
% puff pressures, and collected lid position data in response to 8 air
% puffs directed randomly to the left and right eye. The time-series data
% of lid position is stored in "iFiles" provided by the company.
% Separately, the company provides data files that measure multiple
% parameters of the blink response, and further mark individual trials as
% having produced a "valid" or "invalid" blink (effectively, a blink above
% a certain criterion size).
%
% The analysis pursued here attempts to use all data from all subjects, in
% a departure from our pre-registered plan. The average blink response for
% each subject, pressure level, and session is obtained, including those
% trials that BlinkTBI considers as "invalid" (due to having an
% insufficient blink response). The time-series data is loaded using the
% function "returnBlinkTimeSeries.m". We do, however, remove some trials
% during which subjects engaged in a "squint" behavior, in that they closed
% their eyes and never subsequently opened them.
%
% The data are combined across sessions, and the 17x5 (subject x pressure)
% blink responses are subjected to an analysis using a model derived from
% the amplitude, velocity, and remaining shape variation of the response.

% housekeeping
close all
clear

% There is a large blink response to the first of the 8 trials in each
% acquisition. This flag controls whether this first trial is discarded in
% the analysis.
discardFirstTrialFlag = true;
discardSquintScansFlag = false;
nBoots = 1000;

% Get the location to save plots
plotSaveDir = getpref('blinkCNSAnalysis','plotSaveDir');

% List of subject IDs
subjectIDs = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591, ...
    14590, 14589, 14588, 14587, 14586, 15513, 15514};
nSubs = length(subjectIDs);

% The set of intended PSI values
targetPSISet = [3.5,7.5,15,30,60];
nPSIs = length(targetPSISet);

% A log-transformed version of the PSI values to use for fitting later
xVals = log10(targetPSISet);
xValMid = xVals(3);
xVals = xVals - xValMid;

% The number of time-points
nTimePoints = 161;

% Number of blinks per acquisition
nBlinksPerAcq = 8;

% Number of session
nSessions = 2;

% Define variables to hold the time-series data
X = zeros(nSubs,nPSIs,nTimePoints);
XIpsi = zeros(nSubs,nPSIs,nTimePoints);
XSessIpsi = zeros(nSessions,nSubs,nPSIs,nTimePoints);
XSessContra = zeros(nSessions,nSubs,nPSIs,nTimePoints);
nTrials = zeros(nSubs,5);

% Loop through the subjects and identify the largest excursion seen for
% each of the sessions. This may happen in the 30 or 60 PSI stimulus
for ss=1:nSubs
    for tt=1:nSessions
        maxExcursionIpsi(ss,tt) = -min([...
            min(returnBlinkTimeSeries( subjectIDs{ss}, 30, tt, 'ipsi', discardFirstTrialFlag, discardSquintScansFlag )),...
            min(returnBlinkTimeSeries( subjectIDs{ss}, 60, tt, 'ipsi', discardFirstTrialFlag, discardSquintScansFlag )),...
            ]);
        maxExcursionContra(ss,tt) = -min([...
            min(returnBlinkTimeSeries( subjectIDs{ss}, 30, tt, 'contra', discardFirstTrialFlag, discardSquintScansFlag )),...
            min(returnBlinkTimeSeries( subjectIDs{ss}, 60, tt, 'contra', discardFirstTrialFlag, discardSquintScansFlag )),...
            ]);
    end
end

% Loop through subjects and pressure levels. Load the full and by-session
% data separately. As each time-series is loaded, scale it by the max
% excursion so that the response is in units of proportion of maximal blink
% for that session and eye.
for ss=1:nSubs
    for pp=1:nPSIs
        for tt=1:nSessions
            [XSessIpsi(tt,ss,pp,:),~,~,~,~,XessIpsiBoots] = ...
                returnBlinkTimeSeries( subjectIDs{ss}, targetPSISet(pp), tt, 'ipsi',...
                discardFirstTrialFlag, discardSquintScansFlag, nBoots);
            XSessIpsi(tt,ss,pp,:) = XSessIpsi(tt,ss,pp,:)./ maxExcursionIpsi(ss,tt);
            XSessIpsiBoots{tt,ss,pp} = XessIpsiBoots./ maxExcursionIpsi(ss,tt);

            [XSessContra(tt,ss,pp,:),~,~,~,~,XessContraBoots] = ...
                returnBlinkTimeSeries( subjectIDs{ss}, targetPSISet(pp), tt, 'contra',...
                discardFirstTrialFlag, discardSquintScansFlag, nBoots);
            XSessContra(tt,ss,pp,:) = XSessContra(tt,ss,pp,:)./ maxExcursionContra(ss,tt);
            XSessContraBoots{tt,ss,pp} = XessContraBoots./ maxExcursionContra(ss,tt);
        end
        XIpsi(ss,pp,:) = (XSessIpsi(1,ss,pp,:)+XSessIpsi(2,ss,pp,:))./2;
        XContra(ss,pp,:) = (XSessContra(1,ss,pp,:)+XSessContra(2,ss,pp,:))./2;
        XBoth(ss,pp,:) = ( XIpsi(ss,pp,:) +  XContra(ss,pp,:) )./2;
        XBothBoot{ss,pp} = ( XSessIpsiBoots{1,ss,pp} + XSessIpsiBoots{2,ss,pp} + ...
            XSessContraBoots{1,ss,pp} + XSessContraBoots{2,ss,pp} )./4;
    end

    % Get the raw vector and obtain the blinks averaged by trial index for
    % the analysis of habituation effects. Do not discard the first trial.
    trialBySess = [];
    for tt=1:nSessions
        [~,~,~,blinkVectorRaw,trialIndices] = returnBlinkTimeSeries( subjectIDs{ss}, [], tt, 'ipsi', false );
        for aa=1:nBlinksPerAcq
            trialBySess(tt,aa,:) = nanmean(blinkVectorRaw(trialIndices==aa,:)) ./ maxExcursionIpsi(ss,tt);
        end
    end
    trialX(ss,:,:) = squeeze(mean(trialBySess));

end

% Call the function once more to grab the temporal support
[~,temporalSupport]=returnBlinkTimeSeries( subjectIDs{ss}, targetPSISet(pp), 2);
[~,zeroIdx]=min(abs(temporalSupport));

% Reshape into a matrix
XAll_mat = [reshape(XIpsi,nSubs*nPSIs,nTimePoints); reshape(XContra,nSubs*nPSIs,nTimePoints)];

%% Create a regression matrix

% We will have three covariates
q=3;

% This is the average response and its first derivative
components = [];
components(:,1)=mean(XAll_mat);
components(:,2)=[0 diff(mean(XAll_mat))];

% mean center
components = components - mean(components);

% Obtain the residual of the X_mat matrix after removing these components.
for ii=1:size(XAll_mat,1)
    y = XAll_mat(ii,:);
    offset = mean(y);
    y = y - offset;
    y = y - (regress(y',components)'*components');
    y = y + offset;
    X_resid(ii,:) = y;
end

% Reset the random seed to ensure a consistent PCA result
rng default;

% Conduct a PCA to find the most informative third component
pcaCoeff = pca(X_resid);
components = [components, -pcaCoeff(:,1)];

% Scale these to have unit excursion and mean center
components = components - mean(components);
components = components ./ range(components);

% Conduct the regression and obtain the coefficients.
for ss=1:nSubs
    for pp=1:nPSIs

        y = squeeze(XIpsi(ss,pp,:))';
        offset = mean(y);
        y = y - offset;
        b = regress(y',components)';
        XIpsi_coeff(ss,pp,:) = b;
        XIpsi_fit(ss,pp,:) = b*components' + offset;

        y = squeeze(XContra(ss,pp,:))';
        offset = mean(y);
        y = y - offset;
        b = regress(y',components)';
        XContra_coeff(ss,pp,:) = b;
        XContra_fit(ss,pp,:) = b*components' + offset;

        y = squeeze(XBoth(ss,pp,:))';
        offset = mean(y);
        y = y - offset;
        b = regress(y',components)';
        XBoth_coeff(ss,pp,:) = b;
        XBoth_fit(ss,pp,:) = b*components' + offset;
    
        % Across bootstraps 
        boothMat = XBothBoot{ss,pp};
        for bb=1:nBoots
            y = squeeze(boothMat(bb,:));
            offset = mean(y);
            y = y - offset;
            b = regress(y',components)';
            XBoothBoot_coeff(bb,ss,pp,:) = b;
        end

    end
end

% Reshape
XAll_fit_mat = [reshape(XIpsi_fit,nSubs*nPSIs,nTimePoints); reshape(XContra_fit,nSubs*nPSIs,nTimePoints)];

% Now do the regression across trials
for ii=1:nSubs
    for aa=1:nBlinksPerAcq
        y = squeeze(trialX(ii,aa,:))';
        offset = mean(y);
        y = y - offset;
        b = regress(y',components)';
        trialX_coeff(ii,aa,:) = b;
    end
end
