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
% blink responses are subjected to an Independent Component Analysis,
% intialized with the mean chained derivatives of the average blink
% response.

% housekeeping
close all
clear

% Flag to control if the ICA is intialized with a derivative set
initializeICA = false;

% List of subject IDs
subjectIDs = {15512, 15507, 15506, 15505, 14596, 14595, 14594, 14593, 14592, 14591, ...
    14590, 14589, 14588, 14587, 14586, 15513, 15514};
nSubs = length(subjectIDs);

% The set of intended PSI values
targetPSISet = [3.5,7.5,15,30,60];
nPSIs = length(targetPSISet);

% The number of time-points
nTimePoints = 161;

% Load the time-series data
X = zeros(nSubs,nPSIs,nTimePoints);
X1 = zeros(nSubs,nPSIs,nTimePoints);
X2 = zeros(nSubs,nPSIs,nTimePoints);
nTrials = zeros(nSubs,5);

% Loop through subjects and pressure levels. Load the full and by-session
% data separately for ease of coding below
for ss=1:nSubs
    for pp=1:5
        [X(ss,pp,:),~,nTrials(ss,pp)]=returnBlinkTimeSeries( subjectIDs{ss}, targetPSISet(pp));
        X1(ss,pp,:)=returnBlinkTimeSeries( subjectIDs{ss}, targetPSISet(pp), 1);
        X2(ss,pp,:)=returnBlinkTimeSeries( subjectIDs{ss}, targetPSISet(pp), 2);
    end
end

% Call the function once more to grab a temporal support
[~,temporalSupport]=returnBlinkTimeSeries( subjectIDs{ss}, targetPSISet(pp), 2);

% Reshape into a matrix
X_ICA = reshape(X,nSubs*nPSIs,nTimePoints);
X1_ICA = reshape(X1,nSubs*nPSIs,nTimePoints);
X2_ICA = reshape(X2,nSubs*nPSIs,nTimePoints);

% Remove any "bad" blink averages that are all nans
goodIdx = ~any(isnan(X_ICA'));
X_ICA = X_ICA(goodIdx,:);
X1_ICA = X1_ICA(goodIdx,:);
X2_ICA = X2_ICA(goodIdx,:);

%% Conduct the ICA
% After some trial-and-error, I find that 4 dimensions fits the data
% very well, and supports the creation of independent "amplitude" and
% "timing" components.
rng default % For reproducibility
q = 4; % four dimensions

% Create some initial weights which are smoothed derivatives of the average
% blink response
if initializeICA
initialWeights = zeros(nTimePoints,q);
initialWeights(:,1)=mean(X_ICA);
options = fitoptions('Method','Smooth','SmoothingParam',0.2);
for ii=2:q
    fitObj = fit([1:nTimePoints]',initialWeights(:,ii-1),'smooth',options);
    initialWeights(:,ii) = differentiate(fitObj, 1:1:nTimePoints);
end
initialWeights=initialWeights-mean(initialWeights);
initialWeights=initialWeights./vecnorm(initialWeights);

% ICA time
Mdl = rica(X_ICA,q,'InitialTransformWeights',initialWeights);

else

    Mdl = rica(X_ICA,q);

end

% Derive the coefficients
X_ICAcoeff = Mdl.transform(X_ICA);
X1_ICAcoeff = Mdl.transform(X1_ICA);
X2_ICAcoeff = Mdl.transform(X2_ICA);

% Extract the components
components = Mdl.TransformWeights;

% Flip some signs if we used the initialize weights, which has an effect
% upon the (arbitrary) signs of the ICA components
if initializeICA
    X_ICAcoeff(:,4) = -X_ICAcoeff(:,4);
    X1_ICAcoeff(:,4) = -X1_ICAcoeff(:,4);
    X2_ICAcoeff(:,4) = -X2_ICAcoeff(:,4);
    components(:,4) = -components(:,4);
end

% Generate the fit
X_ICAfit = components*X_ICAcoeff';

% Reshape the results
Xcoeff = nan(nSubs*nPSIs,q);
Xcoeff(goodIdx,:)=X_ICAcoeff;
Xcoeff = reshape(Xcoeff,nSubs,nPSIs,q);
X1coeff = nan(nSubs*nPSIs,q);
X1coeff(goodIdx,:)=X1_ICAcoeff;
X1coeff = reshape(X1coeff,nSubs,nPSIs,q);
X2coeff = nan(nSubs*nPSIs,q);
X2coeff(goodIdx,:)=X2_ICAcoeff;
X2coeff = reshape(X2coeff,nSubs,nPSIs,q);
Xfit = nan(nSubs*nPSIs,nTimePoints);
Xfit(goodIdx,:) = X_ICAfit';
Xfit = reshape(Xfit,nSubs,nPSIs,nTimePoints);

% Fit a slope to the first and fourth component coefficients
for ii=1:nSubs
    ampPuffCoeff(ii,:)=polyfit(-2:2,Xcoeff(ii,:,1),1);
    ampPuffCoeff1(ii,:)=polyfit(-2:2,X1coeff(ii,:,1),1);
    ampPuffCoeff2(ii,:)=polyfit(-2:2,X2coeff(ii,:,1),1);
    speedPuffCoeff(ii,:)=polyfit(-2:2,Xcoeff(ii,:,4),1);
    speedPuffCoeff1(ii,:)=polyfit(-2:2,X1coeff(ii,:,4),1);
    speedPuffCoeff2(ii,:)=polyfit(-2:2,X2coeff(ii,:,4),1);
end

%% Create some plots

% Define a gray-to-red color set for puff-pressure
psiColors = [0.5:0.125:1.0; 0.5:-0.125:0; 0.5:-0.125:0]';

% Average blink response by puff pressure
figure
tmp = squeeze(mean(X,1));
for pp = 1:nPSIs
    plot(temporalSupport,tmp(pp,:),'-','Color',psiColors(pp,:),'LineWidth',1.5)
    hold on
end
xlabel('time [msecs]');
ylabel('blink depth [pixels]');

% Illustration of the ICA components
figure
componentNames = {'amplitude','shape1','shape2','speed'};
componentColors = [0 0 0; 0.85 0.85 0.85; 0.65 0.65 0.65; 0 0 1];
componentWidths = [1.5, 1, 1, 1.5];
plotOrder = [1 4 2 3];
for cc = plotOrder
    plot(temporalSupport,components(:,cc),'-','Color',componentColors(cc,:),'LineWidth',componentWidths(cc))
    hold on
end
legend(componentNames(plotOrder))
xlabel('time [msecs]');
ylabel('componnt value [a.u.]');

% Plot of the coefficients by puff pressure
figure
meanCoeff = squeeze(mean(Xcoeff,1));
semCoeff = squeeze(std(Xcoeff,1))./sqrt(nSubs);
for cc=1:4
    subplot(2,2,cc)
    for pp = 1:nPSIs
        plot([log10(targetPSISet(pp)) log10(targetPSISet(pp))],[meanCoeff(pp,cc)+semCoeff(:,cc),meanCoeff(pp,cc)-semCoeff(:,cc)],'-k');
        hold on
        plot(log10(targetPSISet(pp)),meanCoeff(pp,cc),'o','Color',componentColors(cc,:));
    end
end

% Calculate the correlation of the fit with each average blink response
for ss=1:nSubs
    for pp=1:5
        varExplained(ss,pp) = corr(squeeze(X(ss,pp,:)),squeeze(Xfit(ss,pp,:)))^2';
    end
end

% Show scatter plots of test / retest of overall amplitude and speed
figure
subplot(2,2,1);
plot(ampPuffCoeff1(:,2),ampPuffCoeff2(:,2),'ok'); axis square
refline(1,0);
title(sprintf('amplitude offset, r=%2.2f',corr(ampPuffCoeff1(:,2),ampPuffCoeff2(:,2))));
subplot(2,2,2);
plot(speedPuffCoeff1(:,2),speedPuffCoeff2(:,2),'ob'); xlim([-50 150]); ylim([-50 150]); axis square
refline(1,0);
title(sprintf('speed offset, r=%2.2f',corr(speedPuffCoeff1(:,2),speedPuffCoeff2(:,2))));
subplot(2,2,3);
plot(ampPuffCoeff1(:,1),ampPuffCoeff2(:,1),'ok'); xlim([0 200]); ylim([0 200]);axis square
refline(1,0);
title(sprintf('amplitude slope, r=%2.2f',corr(ampPuffCoeff1(:,1),ampPuffCoeff2(:,1))));
subplot(2,2,4);
plot(speedPuffCoeff1(:,1),speedPuffCoeff2(:,1),'ob'); xlim([-20 80]); ylim([-20 80]);axis square
refline(1,0);
title(sprintf('speed slope, r=%2.2f',corr(speedPuffCoeff1(:,1),speedPuffCoeff2(:,1))));

