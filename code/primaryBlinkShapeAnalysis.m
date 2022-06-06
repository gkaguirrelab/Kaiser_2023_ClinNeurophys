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
    for pp=1:5
        [X(ss,pp,:),~,nTrials(ss,pp)]=returnBlinkTimeSeries( subjectIDs{ss}, targetPSISet(pp));
        X1(ss,pp,:)=returnBlinkTimeSeries( subjectIDs{ss}, targetPSISet(pp), 1);
        X2(ss,pp,:)=returnBlinkTimeSeries( subjectIDs{ss}, targetPSISet(pp), 2);
    end
end

% Call the function once more to grab the temporal support
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

% ICA time
Mdl = rica(X_ICA,q);

% Derive the coefficients
X_ICAcoeff = Mdl.transform(X_ICA);
X1_ICAcoeff = Mdl.transform(X1_ICA);
X2_ICAcoeff = Mdl.transform(X2_ICA);

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


%%%%%%%%%%%%%%%%%%%%%
%% FIGURES
%%%%%%%%%%%%%%%%%%%%%

% Define a gray-to-red color set for puff-pressure
psiColors = [0.5:0.125:1.0; 0.5:-0.125:0; 0.5:-0.125:0]';


%% Acquisition order and example raw set of blinks
psiAcqOrder = [4 4 1 3 5 5 3 1 2 4 2 2 1 5 4 3 3 4 5 2 3 2 5 1 1 4];
figure
set(gcf, 'Position',  [100, 100, 560, 200])
subplot(2,1,1)
plotIdx = 1;
dash = 10;
dot = 2;
for aa=1:length(psiAcqOrder)
  plot( plotIdx:plotIdx+dash-1, repmat(psiAcqOrder(aa),dash,1), '-', 'Color', psiColors(psiAcqOrder(aa),:),'LineWidth',2)
  if aa==1; hold on; end
  plotIdx = plotIdx+dash+dot;
end
xlim([1 plotIdx])
axis off

subplot(2,1,2)
subjectID = 14591;
[~,temporalSupport,nTrials,blinkVectorRaw] = returnBlinkTimeSeries( 14591, [], 1 );
spacing = 10;
plotIdx = (nTimePoints)*nBlinksPerAcq;
blinkIdx = 1;
for aa=2:length(psiAcqOrder)
    for bb = 1:nBlinksPerAcq
        plot( plotIdx:plotIdx+nTimePoints-1, blinkVectorRaw(blinkIdx,:), '-', 'Color', psiColors(psiAcqOrder(aa),:))
        if aa==2 && bb==1; hold on; end
        plotIdx = plotIdx + nTimePoints;
        blinkIdx = blinkIdx+1;
    end
end
xlim([1 plotIdx])
axis off

% Average blink response by puff pressure
figure
tmpX = squeeze(mean(X,1));
tmpXfit = squeeze(mean(Xfit,1));
for pp = 1:nPSIs
    plot(temporalSupport,tmpX(pp,:),'-','Color',psiColors(pp,:),'LineWidth',1.5)
    hold on
    plot(temporalSupport,tmpXfit(pp,:),':','Color',psiColors(pp,:),'LineWidth',1.5)
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
ylabel('component value [a.u.]');

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
    % Add a linear fit line
    if cc==4
        pp = polyfit(xVals(2:end),meanCoeff(2:end,cc),1);
        plot([xVals(2)+xValMid xVals(end)+xValMid],polyval(pp,[xVals(2) xVals(end)]),'-r')
        plot([xVals(1)+xValMid xVals(2)+xValMid],polyval(pp,[xVals(2) xVals(2)]),'-r')
    else
        pp = polyfit(xVals,meanCoeff(:,cc),1);
        plot([xVals(1)+xValMid xVals(end)+xValMid],polyval(pp,[xVals(1) xVals(end)]),'-r')
    end
end


% Calculate the correlation of the fit with each average blink response
for ss=1:nSubs
    for pp=1:5
        varExplained(ss,pp) = corr(squeeze(X(ss,pp,:)),squeeze(Xfit(ss,pp,:)))^2';
    end
end

% Show scatter plots of test / retest of overall amplitude and speed
limVals = {...
    [0 700],[-100 300];...
    [0 1000],[-150 150]};
symbolColors={'k','b'};
nameRow = {'slope','offset'};
nameColumn = {'amplitude','rapidity'};
figure
for cc=1:2
    for rr=1:2
        subplot(2,2,(2-rr)*2+cc);
        if cc==1
            vals1 = ampPuffCoeff1(:,rr);
            vals2 = ampPuffCoeff2(:,rr);
        else
            vals1 = speedPuffCoeff1(:,rr);
            vals2 = speedPuffCoeff2(:,rr);
        end
        scatter(vals1,vals2,'MarkerFaceColor',symbolColors{cc},'MarkerEdgeColor','none','MarkerFaceAlpha',0.5);
        xlim(limVals{rr,cc}); ylim(limVals{rr,cc});
        axis square; box off
        refline(1,0);
        titleStr = sprintf([nameColumn{cc} ' ' nameRow{rr} ' r=%2.2f'],corr(vals1,vals2));
        title(titleStr);
    end
end


%% Illustration of all blink responses and ICA model fit
figure
nPSIs = 5;
betweenSubGap = 3;
[~,zeroIdx]=min(abs(temporalSupport));

% Loop over the three display panels
for xx = 1:3
    switch xx
        case 1
            % Create uniC, which is permuted to order rows by subject then puff
            uniC = reshape(permute(X,[2 1 3]),nSubs*nPSIs,nTimePoints);
            titleStr = 'average blinks';
            C = ones(85+(betweenSubGap*16),161+10,3);
        case 2
            uniC = reshape(permute(Xfit,[2 1 3]),nSubs*nPSIs,nTimePoints);
            titleStr = 'model fit';
            C = ones(85+(betweenSubGap*16),161+10,3);
        case 3
            uniC = reshape(permute(X-Xfit,[2 1 3]),nSubs*nPSIs,nTimePoints);
            titleStr = 'residuals';
            C = ones(85+(betweenSubGap*16),161+10,3);
    end

    % Map uniC to the 0-1 range. We store the scaling factors to use them
    % for all three matrix displays.
    if xx == 1
        minX = min(uniC(:));
        uniC = uniC-minX;
        maxAbsX = max(abs(uniC(:)));
        uniC = uniC ./ maxAbsX;
    else
        uniC = uniC-minX;
        uniC = uniC ./ maxAbsX;
    end

    for ss=1:17
        XrowStart = (ss-1)*nPSIs+1;
        CrowStart = (ss-1)*(nPSIs+betweenSubGap)+1;

        % Place the blink vectors into the matrix
        C(CrowStart:CrowStart+nPSIs-1,11:end,:) = repmat(uniC(XrowStart:XrowStart+4,:),1,1,3);

        % Add a color bar
        C(CrowStart:CrowStart+nPSIs-1,1:10,:) = permute(repmat(psiColors,1,1,10),[1 3 2]);

        % Add a marker for stimulus onset
        C(CrowStart:CrowStart+nPSIs-1,10+zeroIdx,:) = repmat([0 0 1],nPSIs,1);

    end

    subplot(1,3,xx);
    image(imresize(C,[size(C,1)*4,size(C,2)],"nearest")); axis off
    axis equal
    title(titleStr);
end

