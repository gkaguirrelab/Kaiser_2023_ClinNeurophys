%% primaryBlinkShapeFigures
% Run after primaryBlinkShapeAnalysis.m

% Define a gray-to-red color set for puff-pressure
psiColors = [0.5:0.125:1.0; 0.5:-0.125:0; 0.5:-0.125:0]';

% Define a blue-to-gray color set for trial number
trialColors = [0:(0.5/7):0.5; 0:(0.5/7):0.5; 1.0:-(0.5/7):0.5]';


%% Report the variance explained by the model
tVar = var(X_mat(:));
fVar = var(X_matfit(:));
fprintf('Proportion variance explained by the model is: %2.2f \n',fVar/tVar);


%% Acquisition order and example raw set of blinks
psiAcqOrder = [4 4 1 3 5 5 3 1 2 4 2 2 1 5 4 3 3 4 5 2 3 2 5 1 1 4];
figure
set(gcf, 'Position',  [100, 100, 400, 200])
subplot(3,2,1:2)
plotIdx = 1;
dash = 10;
dot = 2;
for aa=1:length(psiAcqOrder)
    plot( plotIdx:plotIdx+dash-1, repmat(psiAcqOrder(aa),dash,1), '-', 'Color', psiColors(psiAcqOrder(aa),:),'LineWidth',4)
    if aa==1; hold on; end
    plotIdx = plotIdx+dash+dot;
end
xlim([1 plotIdx])
axis off

subplot(3,2,3:4)
subjectID = 14594;
[~,~,~,blinkVectorRaw] = returnBlinkTimeSeries( subjectID, [], 1 );
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
plot([0 0],[0 -250],'-','Color',[0.5 0.5 0.5],'LineWidth',2)
t=text(-600,-250,'250 pixels');
t.Rotation = 90;
plot([0 1000],[-250 -250],'-','Color',[0.5 0.5 0.5],'LineWidth',2)
t=text(0,-150,'1s');
xlim([1 plotIdx])
axis off
title(num2str(subjectID));

for ss=1:2
    subplot(3,2,ss+4)
    for pp=1:nPSIs
        plot(temporalSupport, returnBlinkTimeSeries( subjectID, targetPSISet(pp), ss ), '-', 'Color', psiColors(pp,:),'LineWidth',1.5);
        if pp==1;hold on; end
    end
    axis off
    plot([0 0],[-125 25],'-b')
    title(sprintf('Session %d',ss));
    if ss==1
        plot([-100 -100],[0 -100],'-','Color',[0.5 0.5 0.5],'LineWidth',2)
        plot([-100 0],[-125 -125],'-','Color',[0.5 0.5 0.5],'LineWidth',2)
        t=text(-100,-50,'100 msec');
        t=text(-175,-100,'100 pixels');
        t.Rotation = 90;
    end
end
saveas(gcf,fullfile(plotSaveDir,'acquisitionOrder.pdf'));


%% Average blink response by puff pressure
figure
set(gcf, 'Position',  [100, 100, 600, 200])
tmpX = squeeze(mean(X,1));
tmpXfit = squeeze(mean(Xfit,1));
tmpXfit = tmpXfit - tmpXfit(:,1);
for pp = 1:nPSIs
    subplot(1,2,1)
    plot(temporalSupport,tmpX(pp,:),'-','Color',psiColors(pp,:),'LineWidth',1.5)
    hold on
    subplot(1,2,2)
    plot(temporalSupport,tmpXfit(pp,:),'-','Color',psiColors(pp,:),'LineWidth',1.5)
    hold on
end
for mm=1:2
    subplot(1,2,mm)
    plot([0 0],[-1 0],'-b')
    plot([-100 -100],[-1 0],'-','Color',[0.5 0.5 0.5],'LineWidth',2)
    plot([-100 0],[-1 -1],'-','Color',[0.5 0.5 0.5],'LineWidth',2)
    xlabel('time [msecs]');
    ylabel('blink [proportion]');
end
saveas(gcf,fullfile(plotSaveDir,'averageBlnkResponseByPSI.pdf'));


%% Illustration of the model components
figure
componentNames = {'amplitude','velocity','shape'};
componentColors = [0 0 0; 0 0 1; 0.65 0.65 0.65];
componentWidths = [1.5, 1.5, 1];
plotOrder = [1 2 3];
for cc = plotOrder
    plot(temporalSupport,components(:,cc)-components(1,cc),'-','Color',componentColors(cc,:),'LineWidth',componentWidths(cc))
    hold on
end
legend(componentNames(plotOrder))
xlabel('time [msecs]');
ylabel('component value');
saveas(gcf,fullfile(plotSaveDir,'ModelComponents.pdf'));


%% Plot of the coefficients by puff pressure

% Calculate the latency shift in msecs implied by one unit of the velocity
% model component
deltaVelocityComponent = 0.1;
scaleVals = -deltaVelocityComponent*2:deltaVelocityComponent:deltaVelocityComponent*2;
tFine = 0:0.1:max(temporalSupport);
for ii=1:length(scaleVals)
    vec = components(:,1)+components(:,2)*scaleVals(ii);
    vec = vec-vec(1);
    vecFine = interp1(temporalSupport,vec,tFine,'cubic');
    [~,idx] = min( abs(vecFine+0.5) );
    tPoint(ii)=tFine(idx);
end
velocityFactor = mean(diff(tPoint))/deltaVelocityComponent;

% Time to plot
figure
set(gcf, 'Position',  [100, 100, 200, 600])
meanCoeff = squeeze(mean(Xcoeff,1));
semCoeff = squeeze(std(Xcoeff,1))./sqrt(nSubs);
plotOrder = [1 2 3];
for cc=1:3
    subplot(3,1,plotOrder(cc))
    for pp = 1:nPSIs
        plot([log10(targetPSISet(pp)) log10(targetPSISet(pp))],[meanCoeff(pp,cc)+2.*semCoeff(:,cc),meanCoeff(pp,cc)-2.*semCoeff(:,cc)],'-k');
        hold on
        plot(log10(targetPSISet(pp)),meanCoeff(pp,cc),'o',...
            'MarkerFaceColor',componentColors(cc,:),'MarkerEdgeColor','none' );
    end
    xticks(log10(targetPSISet));
    xticklabels(arrayfun(@num2str, targetPSISet, 'UniformOutput', 0));
    xlabel('stimulus intensity [log PSI]')
    title(componentNames{cc})
    box off
    switch cc
        case 1
            ylabel('proportion blink');
        case 2
            % Convert the coefficient values to units of msecs. A
            % coefficient value equivalent to a 5 msec latency shift is
            coeffVal5msec = -5/velocityFactor;
            ylim([-coeffVal5msec*1.5 coeffVal5msec*1.75]);
            axHandle = gca;
            axHandle.YTick = [-coeffVal5msec 0 coeffVal5msec];
            axHandle.YTickLabel = {'5','0','-5'};
            ylabel('latency shift [msecs]');
        case 3
            ylim([-0.2 0.2]);
            ylabel('arbitrary units');
    end
end
saveas(gcf,fullfile(plotSaveDir,'coefficientsByPSI.pdf'));


%% Fit a weibullCDF to the amplitude data

% Define an increasing weibull CDF with 4 parameters. Assume that there is
% a zero amplitude response at some zero stimulus.
figure
set(gcf, 'Position',  [100, 100, 200, 800])
xShift = log10(0.01);
x = log10(targetPSISet)-xShift;
deltaX = x(2)-x(1);
xFit = linspace(xShift,log10(100),1000)-xShift;
deltaXFit = xFit(2)-xFit(1);
lb = [0 1 2 5];
ub = [0 1 5 15];
p0 = [0 1 2.5 6];
weibullCDF = @(x,p) p(1) + p(2) - p(2)*exp( - (x./p(3)).^p(4) ) ;
options = optimoptions('fmincon','Display','off');
for ii=1:nSubs

    % Ipsilateral response
    y = Xcoeff(ii,:,1);
    y = y ./ max(y);
    if y(nPSIs-1)>y(nPSIs)
        weights = [ones(1,nPSIs-1) 0.5];
    else
        weights = ones(1,nPSIs);
    end
    myObj = @(p) norm((y-weibullCDF(x,p)).*weights);
    p(ii,:) = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);
    yFit = weibullCDF(xFit,p(ii,:));
    [~,idx] = min(abs(yFit-0.5));
    x50(ii) = 10^(xFit(idx)+xShift);
    slopeAt50(ii) = max(diff(yFit))/deltaXFit;

    % Contralateral response
    y = XContraCoeff(ii,:,1);
    y = y ./ max(y);
    if y(nPSIs-1)>y(nPSIs)
        weights = [ones(1,nPSIs-1) 0.5];
    else
        weights = ones(1,nPSIs);
    end
    myObj = @(p) norm((y-weibullCDF(x,p)).*weights);
    pContra(ii,:) = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);
    yFit = weibullCDF(xFit,pContra(ii,:));
    [~,idx] = min(abs(yFit-0.5));
    x50Contra(ii) = 10^(xFit(idx)+xShift);
    slopeAt50Contra(ii) = max(diff(yFit))/deltaXFit;

    % Now fit by session
    for tt=1:nSessions

        y = squeeze(XSessCoeff(tt,ii,:,1))';
        y = y ./ max(y);
        if y(nPSIs-1)>y(nPSIs)
            weights = [ones(1,nPSIs-1) 0.5];
        else
            weights = ones(1,nPSIs);
        end

        myObj = @(p) norm((y-weibullCDF(x,p)).*weights);
        pSess(tt,ii,:) = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);
        yFit = weibullCDF(xFit,squeeze(pSess(tt,ii,:))');
        [~,idx] = min(abs(yFit-0.5));
        x50Sess(tt,ii) = 10^(xFit(idx)+xShift);
        slopeAt50Sess(tt,ii) = max(diff(yFit))/deltaXFit;

    end

end

% Plot the subject fits in order of sensitivity
[~,sortOrder]=sort(x50,'descend');
for ss=1:nSubs
    ii = sortOrder(ss);
    y=Xcoeff(ii,:,1);
    y = y ./ max(y);
    yFit = weibullCDF(xFit,p(ii,:));
    subplot(9,2,ss)
    plot(xFit+xShift,yFit,'-r');
    hold on;
    plot(x+xShift,y,'ok');
    [~,idx] = min(abs(yFit-0.5));
    plot([xFit(idx)+xShift, xFit(idx)+xShift],[0 0.5],':b')
    plot([xShift, xFit(idx)+xShift],[0.5 0.5],':b')
    box off
    xlim([-1 2.1]);
    ylim([0 1.2]);
    axHandle = gca;
    axHandle.XTickLabel = cellstr(string(10.^axHandle.XTick));
    if mod(ss-1,6)>0
        set(gca,'ytick',[])
    end
    if ss==1
        xlabel('stimulus [PSI]');
        ylabel('proportion');
    end
end
saveas(gcf,fullfile(plotSaveDir,'weibullFitToAmplitude.pdf'));



%% Illustrate individual differences in amplitude parameters
figure
subIdxToShow = [sortOrder(2),sortOrder(end-1)];
for ss=1:length(subIdxToShow)
    subplot(1,2,ss);
    for pp=1:nPSIs
        plot(temporalSupport, squeeze(X( subIdxToShow(ss), pp, :)), 'Color', psiColors(pp,:),'LineWidth',1);
        hold on
    end
    ylim([-1.1 0.2]);
    box off
end
saveas(gcf,fullfile(plotSaveDir,'blinkResponseSubjectExtremeX50.pdf'));


%% Scatter plot of test / retest of x50 on amplitude
figure
for pp=1:2
    subplot(1,2,pp)
    switch pp
        case 1
            vals = log10(x50Sess);
        case 2
            vals = slopeAt50Sess;
    end
    scatter(vals(1,:),vals(2,:),100,'r');
    [Rval,Pval] = corrcoef(vals(1,:),vals(2,:));
    titleStr = sprintf('corr param r=%2.2f, p=%2.5f',Rval(2),Pval(2));
    title(titleStr);
    axis square; box off
    switch pp
        case 1
            xlabel('S1 half-response stimulus [PSI]');
            ylabel('S2 half-response stimulus [PSI]');
            a=gca;
            xlim([0 1.7]);
            ylim([0 1.7]);
            axis square
            a.XTick = log10([1 5 10 50]);
            a.YTick = log10([1 5 10 50]);
            a.XTickLabel = {'1','5','10','50'};
            a.YTickLabel = {'1','5','10','50'};
        case 2
            xlabel({'S1 max slope','[proportion / log PSI]'});
            ylabel({'S2 max slope','[proportion / log PSI]'});
    end
    refline(1,0);
end
saveas(gcf,fullfile(plotSaveDir,'testRetestCoefficients.pdf'));


%% Illustration of ipsi vs. contra sensitivity
faceColors = {'k','none'};
edgeColors = {'none','k'};
lineColors = {'-r','--r'};
figure
set(gcf, 'Position',  [100, 100, 840, 420])

subplot(1,3,1)
for cc = 1:2
    switch cc
        case 1
            y = squeeze(mean(Xcoeff(:,:,1),1));
            ySEM = squeeze(std(Xcoeff,1))./sqrt(nSubs);
        case 2
            y = squeeze(mean(XContraCoeff(:,:,1),1));
            ySEM = squeeze(std(XContraCoeff,1))./sqrt(nSubs);
    end
    y = y./max(y);
    if y(nPSIs-1)>y(nPSIs)
        weights = [ones(1,nPSIs-1) 0.5];
    else
        weights = ones(1,nPSIs);
    end
    myObj = @(p) norm((y-weibullCDF(x,p)).*weights);
    pFit = fmincon(myObj,p0,[],[],[],[],lb,ub,[],options);
    yFit = weibullCDF(xFit,pFit);
    for pp = 1:nPSIs
        plot(x(pp)+xShift,y(pp),'o',...
            'MarkerFaceColor',faceColors{cc},'MarkerEdgeColor',edgeColors{cc} );
        hold on
    end
    plot(xFit+xShift,yFit,lineColors{cc});
end
xticks(log10(targetPSISet));
xticklabels(arrayfun(@num2str, targetPSISet, 'UniformOutput', 0));
xlabel('stimulus intensity [log PSI]')
ylabel('proportion blink');
ylim([0 1.1]);
xlim([0.1 2]);
axis square

subplot(1,3,2)
scatter(log10(x50),log10(x50Contra),100,'r');
axis square; box off
xlabel('ipsi half-response stimulus [PSI]');
ylabel('contra half-response stimulus [PSI]');
a=gca;
xlim([0 1.7]);
ylim([0 1.7]);
axis square
a.XTick = log10([1 5 10 50]);
a.YTick = log10([1 5 10 50]);
a.XTickLabel = {'1','5','10','50'};
a.YTickLabel = {'1','5','10','50'};
axis square
refline(1,0);
saveas(gcf,fullfile(plotSaveDir,'ipsiVsContraSensitivity.pdf'));

subplot(1,3,3)
scatter(slopeAt50,slopeAt50Contra,100,'r');
axis square; box off
xlabel({'ipsi max slope','[proportion / log PSI]'});
ylabel({'contra max slope','[proportion / log PSI]'});
a=gca;
xlim([0.5 2]);
ylim([0.5 2]);
axis square
refline(1,0);
saveas(gcf,fullfile(plotSaveDir,'ipsiVsContraSensitivity.pdf'));

% Report the ttest ipsi vs. contra
[~,Tpval,~,Tstats] = ttest(log10(x50),log10(x50Contra));
fprintf('T-test ipsi vs. contral x50 vals: means = [%2.1f, %2.1f], t(df)=%2.2f (%d), p=%2.9f \n',...
    10^mean(log10(x50)),10^mean(log10(x50Contra)),Tstats.tstat,Tstats.df,Tpval);
[~,Tpval,~,Tstats] = ttest(slopeAt50,slopeAt50Contra);
fprintf('T-test ipsi vs. contral max slope vals: means = [%2.1f, %2.1f], t(df)=%2.2f (%d), p=%2.9f \n',...
    mean(slopeAt50),mean(slopeAt50Contra),Tstats.tstat,Tstats.df,Tpval);



%% Plot of the coefficients by trial number
figure
set(gcf, 'Position',  [100, 100, 200, 600])
meanCoeff = squeeze(mean(trialX_coeff,1));
semCoeff = squeeze(std(trialX_coeff,1))./sqrt(nSubs);
plotOrder = [1 2 3];
for cc=1:length(plotOrder)
    subplot(3,1,plotOrder(cc))
    for pp = 1:nBlinksPerAcq
        plot([pp pp],[meanCoeff(pp,cc)+2.*semCoeff(:,cc),meanCoeff(pp,cc)-2.*semCoeff(:,cc)],'-k');
        hold on
        plot(pp,meanCoeff(pp,cc),'o',...
            'MarkerFaceColor',componentColors(cc,:),'MarkerEdgeColor','none' );
    end
    switch cc
        case 1
            ylabel('proportion blink');
            myExpFunc = @(p,x) p(1).*exp(-( x./p(2) )) + p(3);
            myObj = @(p) norm(meanCoeff(:,1)' - myExpFunc(p,1:nBlinksPerAcq));
            myFitP = fmincon(myObj,[1 1 1],[],[],[],[],[],[],[],options);
            plot(1:0.1:nBlinksPerAcq,myExpFunc(myFitP,1:0.1:nBlinksPerAcq),'-r')
        case 2
            myExpFunc = @(p,x) p(3) - p(1).*exp(-( x./p(2) ));
            myObj = @(p) norm(meanCoeff(:,2)' - myExpFunc(p,1:nBlinksPerAcq));
            myFitP = fmincon(myObj,[1 1 1],[],[],[],[],[],[],[],options);
            plot(1:0.1:nBlinksPerAcq,myExpFunc(myFitP,1:0.1:nBlinksPerAcq),'-r')
            % Convert the coefficient values to units of msecs. A
            % coefficient value equivalent to a 5 msec latency shift is
            coeffVal5msec = -5/velocityFactor;
            ylim([-coeffVal5msec*1.5 coeffVal5msec*1.5]);
            axHandle = gca;
            axHandle.YTick = [-coeffVal5msec 0 coeffVal5msec];
            axHandle.YTickLabel = {'5','0','-5'};
            ylabel('latency shift [msecs]');
        case 3
            myExpFunc = @(p,x) p(1).*exp(-( x./p(2) )) + p(3);
            myObj = @(p) norm(meanCoeff(:,3)' - myExpFunc(p,1:nBlinksPerAcq));
            myFitP = fmincon(myObj,[1 1 1],[],[],[],[],[],[],[],options);
            plot(1:0.1:nBlinksPerAcq,myExpFunc(myFitP,1:0.1:nBlinksPerAcq),'-r')
            ylim([-0.2 0.2]);
            ylabel('arbitrary units');
    end
    xticks(1:8);
    xlabel('trial number')
    title(componentNames{cc})
    box off
end
saveas(gcf,fullfile(plotSaveDir,'coefficientsByTrialNumber.pdf'));


%% Plot of the mean response by trial number
figure
meanResponseByTrial=squeeze(mean(trialX));
for ii=1:8
    plot(temporalSupport,meanResponseByTrial(ii,:),'-','Color',trialColors(ii,:),'LineWidth',2);
    hold on
end
saveas(gcf,fullfile(plotSaveDir,'meanResponseByTrialNumber.pdf'));


%% Illustration of all blink responses and model fit
betweenSubGap = 1;
grayScaleRangePixels = [0.25 -1];

% Loop over the three display panels
for xx = 1:3
    figure

    % Create a blank image
    C = ones(85+(betweenSubGap*16),161+10,3);

    % Depending upon xx, pick a variable to display
    switch xx
        case 1
            % Create uniC, which is permuted to order rows by subject then puff
            uniC = reshape(permute(X,[2 1 3]),nSubs*nPSIs,nTimePoints);
            titleStr = 'average blinks';
        case 2
            uniC = reshape(permute(Xfit,[2 1 3]),nSubs*nPSIs,nTimePoints);
            titleStr = 'model fit';
        case 3
            uniC = reshape(permute(X-Xfit,[2 1 3]),nSubs*nPSIs,nTimePoints);
            titleStr = 'residuals';
    end

    % Map uniC to the 0-1 range. We store the scaling factors to use them
    % for all three matrix displays.
    uniC(uniC>grayScaleRangePixels(1))=grayScaleRangePixels(1);
    uniC(uniC<grayScaleRangePixels(2))=grayScaleRangePixels(2);
    uniC = (uniC-grayScaleRangePixels(1));
    uniC = 1+(uniC ./ range(grayScaleRangePixels));
    if xx==1
        grayAtZeroPixels = mean(mean(uniC(:,1:zeroIdx)));
    end

    for ss=1:17
        XrowStart = (ss-1)*nPSIs+1;
        CrowStart = (ss-1)*(nPSIs+betweenSubGap)+1;

        % Place the blink vectors into the matrix
        C(CrowStart:CrowStart+nPSIs-1,11:end,:) = repmat(uniC(XrowStart:XrowStart+4,:),1,1,3);

        % Add a color bar
        C(CrowStart:CrowStart+nPSIs-1,1:7,:) = permute(repmat(psiColors,1,1,7),[1 3 2]);

        % Add a marker for stimulus onset
        C(CrowStart:CrowStart+nPSIs-1,10+zeroIdx,:) = repmat([0 0 1],nPSIs,1);

    end

    image(imresize(C,[size(C,1)*4,size(C,2)],"nearest")); axis off
    axis equal
    exportgraphics(gca,fullfile(plotSaveDir,sprintf('blinkAndFitAllSubjects_%d.png',xx)));
end
