%% tableBlinkRawVals
% Run after primaryBlinkPCnalysis.m
% This routine derives a set of summary measures of the blink response by
% refernce to the "raw" blink response time series. This was added to
% satisfy reviewers who wished to have an account of the data that had not
% passed through the PCA approach.

% Define a gray-to-red color set for puff-pressure
psiColors = [0.5:0.125:1.0; 0.5:-0.125:0; 0.5:-0.125:0]';

% Get some data property values
nSubs = size(XBoth,1);
nLevels = size(XBoth,2);
nSamples = size(XBoth,3);

% Define some variables to hold the values across bootstraps
bootAmpByLevel = [];
bootLatencyByLevel = [];
bootFWHMByLevel = [];

% Define the blink threhsold for the latency measure
thresh = -0.1;

% Bootstrap across subjects and derive some blink response metrics
nBoots = 100;
for bb=1:nBoots

    % Boot-strap resample the data across subjects
    bootIdx = datasample(1:nSubs,nSubs);
    tmpX = squeeze(mean(XBoth(bootIdx,:,:),1));

    % Loop over the levels
    for ss = 1:nLevels

        % Get the response vector at this level
        vec = tmpX(ss,:);

        % Get the amplitude
        [bootAmpByLevel(bb,ss),minIdx] = min(vec);

        % Calculate the latency to threshold
        idx = find(vec<thresh,1);        
        bootLatencyByLevel(bb,ss) = temporalSupport(idx);

        % Calculate temporal FWHM
        halfMax = min(vec)/2;

        % The idx at which the response initially passes the half max on
        % the way down and the way up
        idx1 = find(vec<halfMax,1);        
        vec2 = vec; vec2(1:minIdx) = -1;
        idx2 = find(vec2>halfMax,1);        
        bootFWHMByLevel(bb,ss) = temporalSupport(idx2)-temporalSupport(idx1);

    end

end

% Print a table with the results
bootAmpByLevelSEM = std(bootAmpByLevel);
bootAmpByLevel = mean(bootAmpByLevel);
bootLatencyByLevelSEM = std(bootLatencyByLevel);
bootLatencyByLevel = mean(bootLatencyByLevel);
bootFWHMByLevelSEM = std(bootFWHMByLevel);
bootFWHMByLevel = mean(bootFWHMByLevel);


str = '\nstim level [PSI]\tproportion blink\tlatency [ms]\t FWHM [ms]\n';
fprintf(str);
for ss = 1:nLevels
    v(1) = targetPSISet(ss);
    v(2) = -bootAmpByLevel(ss); v(3) = bootAmpByLevelSEM(ss);
    v(4) = bootLatencyByLevel(ss); v(5) = bootLatencyByLevelSEM(ss);
    v(6) = bootFWHMByLevel(ss); v(7) = bootFWHMByLevelSEM(ss);
    str = sprintf('%2.1f\t%2.2f ± %2.1f\t%2.1f ± %2.1f\t%2.1f ± %2.1f\n',v);
    fprintf(str);
end

str = ['\nSummary measures of average blink response. Errors are SEM by\n' ...
    'bootstrap resampling across subjects. The latency is time in ms to\n' ...
    'reach 10 percent closure. FWHM is the temporal width of the blink at \n' ...
    'half the maxmimum blink for that puff pressure.\n'];
fprintf(str);