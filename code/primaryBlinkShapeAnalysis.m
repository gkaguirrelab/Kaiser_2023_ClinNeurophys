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

% Load the time-series data
X = zeros(nSubs,nPSIs,nTimePoints);
X1 = zeros(nSubs,nPSIs,nTimePoints);
X2 = zeros(nSubs,nPSIs,nTimePoints);
nTrials = zeros(nSubs,5);

% Loop through subjects and pressure levels. Load the full and by-session
% data separately for ease of coding below
for ss=1:nSubs
    for pp=1:nPSIs
        [X(ss,pp,:),~,nTrials(ss,pp)]=returnBlinkTimeSeries( subjectIDs{ss}, targetPSISet(pp));
        X1(ss,pp,:)=returnBlinkTimeSeries( subjectIDs{ss}, targetPSISet(pp), 1);
        X2(ss,pp,:)=returnBlinkTimeSeries( subjectIDs{ss}, targetPSISet(pp), 2);
    end

    % Get the raw vector and obtain the blinks averaged by trial index for
    % the analysis of habituation effects
    [~,~,~,blinkVectorRaw,trialIndices] = returnBlinkTimeSeries( subjectIDs{ss} );
    for tt=1:nBlinksPerAcq
        trialX(ss,tt,:) = nanmean(blinkVectorRaw(trialIndices==tt,:));
    end

end

% Call the function once more to grab the temporal support
[~,temporalSupport]=returnBlinkTimeSeries( subjectIDs{ss}, targetPSISet(pp), 2);
[~,zeroIdx]=min(abs(temporalSupport));

% Reshape into a matrix
X_ICA = reshape(X,nSubs*nPSIs,nTimePoints);
X1_ICA = reshape(X1,nSubs*nPSIs,nTimePoints);
X2_ICA = reshape(X2,nSubs*nPSIs,nTimePoints);
trialX_ICA = reshape(trialX,nSubs*nBlinksPerAcq,nTimePoints);

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

% ICA time
Mdl = rica(X_ICA,q);

% Derive the coefficients
X_ICAcoeff = Mdl.transform(X_ICA);
X1_ICAcoeff = Mdl.transform(X1_ICA);
X2_ICAcoeff = Mdl.transform(X2_ICA);
trialX_ICAcoeff = Mdl.transform(trialX_ICA);

% Extract the components
components = Mdl.TransformWeights;

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
trialXcoeff = reshape(trialX_ICAcoeff,nSubs,nBlinksPerAcq,q);

Xfit = nan(nSubs*nPSIs,nTimePoints);
Xfit(goodIdx,:) = X_ICAfit';
Xfit = reshape(Xfit,nSubs,nPSIs,nTimePoints);

% Fit a slope to the first and fourth component coefficients
for ii=1:nSubs
    ampPuffCoeff(ii,:)=polyfit(xVals,Xcoeff(ii,:,1),1);
    ampPuffCoeff1(ii,:)=polyfit(xVals,X1coeff(ii,:,1),1);
    ampPuffCoeff2(ii,:)=polyfit(xVals,X2coeff(ii,:,1),1);
    speedPuffCoeff(ii,:)=polyfit(xVals(2:end),Xcoeff(ii,2:end,4),1);
    speedPuffCoeff1(ii,:)=polyfit(xVals(2:end),X1coeff(ii,2:end,4),1);
    speedPuffCoeff2(ii,:)=polyfit(xVals(2:end),X2coeff(ii,2:end,4),1);
end

% Fit a weibullCDF to the amplitude data

% Define an increasing weibull CDF with 4 parameters. Assume that there is
% a zero amplitude response at some zero stimulus.
figure
xShift = log10(0.01);
x = log10(targetPSISet)-xShift;
deltaX = x(2)-x(1);
xFit = linspace(xShift,log10(100))-xShift;
deltaXFit = xFit(2)-xFit(1);
lb = [0 0 0 0];
ub = [0 nan 5 12];
weibullCDF = @(x,p) p(1) + p(2) - p(2)*exp( - (x./p(3)).^p(4) ) ;
options = optimoptions('fmincon','Display','off');
for ii=1:nSubs
    y=Xcoeff(ii,:,1);
    % Determine the max amplitude from y
    if y(nPSIs)<y(nPSIs-1)
        maxY(ii) = mean(y(nPSIs-1:nPSIs));
    else
        maxY(ii) = y(nPSIs);
    end
    ub(2) = maxY(ii);
    myObj = @(p) norm(y-weibullCDF(x,p));
    p(ii,:) = fmincon(myObj,[0 1000 1 1],[],[],[],[],lb,ub,[],options);
    yFit = weibullCDF(xFit,p(ii,:))./p(ii,2);
    subplot(2,9,ii);
    plot(x+xShift,y./p(ii,2),'ok'); hold on; plot(xFit+xShift,yFit,'-r');
    [~,idx] = min(abs(yFit-0.5));
    plot([xFit(idx)+xShift xFit(idx)+xShift],[0 0.5],'-m');
    plot([xShift xFit(idx)+xShift],[0.5 0.5],'-m');
    x50(ii) = 10^(xFit(idx)+xShift);
    xlim([-1.1 2.1]);
    ylim([-0.1 1.1]);
    axHandle = gca;
    axHandle.XTickLabel = cellstr(string([0.1 1 10 100]));
end

% Repeat now and obtain the coefficients from the first and second session
for ii=1:nSubs
    for ss = 1:2
        switch ss
            case 1
                y=X1coeff(ii,:,1);
            case 2
                y=X2coeff(ii,:,1);
        end
    ub(2) = maxY(ii);
    myObj = @(p) norm(y-weibullCDF(x,p));
    pSess(ss,ii,:) = fmincon(myObj,[0 1000 1 1],[],[],[],[],lb,ub,[],options);
    yFit = weibullCDF(xFit,p(ii,:))./p(ii,2);
    [~,idx] = min(abs(yFit-0.5));
    x50Sess(ss,ii) = 10^(xFit(idx)+xShift);
    end
end
