

%% timeSeriesAnalysis
% Loads I-Files and conducts an analysis of time series data

%% set up parameters

% load file path

dataPath = fileparts(fileparts(mfilename('fullpath')));

subject = 14591;
location = fullfile(dataPath,'data','iFiles',num2str(subject));

% define scan numbers at each PSI level
psiArray = [];
psiArray(:,1) = [3 8 13 24 25 29 34 39 50 51];
psiArray(:,2) = [9 11 12 20 22 35 37 38 46 48];
psiArray(:,3) = [4 7 16 17 21 30 33 42 43 47];
psiArray(:,4) = [2 10 15 18 26 28 36 41 44 52];
psiArray(:,5) = [5 6 14 19 23 31 32 40 45 49];
psi = ["3.5" "7.5" "15" "30" "60"];

% decide row bounds
numBefore = 60;
numAfter = 700;
       
%% plot stimulus arrival time as a function of puff pressure
ds = tabularTextDatastore(location,'FileExtensions',[".csv"]);
ds.ReadSize = 'file';
x = [];
y = [];

% loop through scan files
for ii = 1:52
   % make table from file
   [data,info] = read(ds);
   T = data;
   [~,subStr] = fileparts(info.Filename);
   subStr = split(subStr,'_'); subStr = subStr{end};
   scan = str2num(subStr);
   allVarNames = T.Properties.VariableNames;
   
   % find stimulus arrival
   onset = vertcat(T(ismember(T.Stimulus,'TA-OD'),:), T(ismember(T.Stimulus,'TA-OS'),:));
   arrival = vertcat(T(ismember(T.Stimulus,'MC-OD'),:), T(ismember(T.Stimulus,'MC-OS'),:));
   tEnd = table2array(arrival(:,1));
   tStart = table2array(onset(:,1));
   
   if ismember(scan, [3 8 13 24 25 29 34 39 50 51])
       % 3.5 PSI
       for jj = 1:length(tEnd)
           x(end+1) = 3.5;
           y(end+1) = tEnd(jj) - tStart(jj);
       end
   elseif ismember(scan, [9 11 12 20 22 35 37 38 46 48])
       % 7.5 PSI       
       for jj = 1:length(tEnd)
           x(end+1) = 7.5;
           y(end+1) = tEnd(jj) - tStart(jj);
       end
   elseif ismember(scan, [4 7 16 17 21 30 33 42 43 47])
       % 15 PSI
       for jj = 1:length(tEnd)
           x(end+1) = 15;
           y(end+1) = tEnd(jj) - tStart(jj);
       end
   elseif ismember(scan, [2 10 15 18 26 28 36 41 44 52])
       % 30 PSI
       for jj = 1:length(tEnd)
           x(end+1) = 30;
           y(end+1) = tEnd(jj) - tStart(jj);
       end
   elseif ismember(scan, [5 6 14 19 23 31 32 40 45 49])
       % 60 PSI
       for jj = 1:length(tEnd)
           x(end+1) = 60;
           y(end+1) = tEnd(jj) - tStart(jj);
       end
   end
end

figure();
x = log10(x);
scatter(x,y);
hold on
fitObj = fitlm(x,y,'RobustOpts', 'on');
plot(x,fitObj.Fitted,'-r')
xlim(log10([2 100]));
title('Stimulus arrival time across puff pressure', 'FontSize', 16);
xlabel(['Puff pressure [log psi]'], 'FontSize', 16);
ylabel(['Time [msec]'], 'FontSize', 16);
hold off

%% plot time series example from single scan
ds = tabularTextDatastore(location,'FileExtensions',[".csv"]);
ds.ReadSize = 'file';
scanNum = 5;
x = [];
y = [];

% loop through scan files
for ii = 1:52
   % make table from file
   [data,info] = read(ds);
   T = data;
   [~,subStr] = fileparts(info.Filename);
   subStr = split(subStr,'_'); subStr = subStr{end};
   scan = str2num(subStr);
   allVarNames = T.Properties.VariableNames;
   
   if scan == scanNum
       % find stimulus arrivals
       [rights,col,v] = find(strcmp('MC-OD',T.Stimulus(:,1)));
       [lefts,col,v] = find(strcmp('MC-OS',T.Stimulus(:,1)));
       all = sort(cat(1,rights,lefts));
       
       % get times
       starts = all - numBefore;
       ends = all + numAfter;
       time = table2array(T(starts(1):ends(1),1));
       time = time - time(numBefore + 1);
       
       % get means across trials
       pos = [];
       for jj = 1:length(starts)
           temp = T(starts(jj):ends(jj),:);
           if ismember(all(jj),rights)
               pos(:,end+1) = table2array(temp(:,3))';
           else
               pos(:,end+1) = table2array(temp(:,4))';
           end
       end
       pos = mean(pos,2);
       
       % center pre-stimulus around zero
       pre = pos(1:numBefore);
       mm = mean(pre);
       pos = pos - mm;
   end
end

% calculate velocity plot
dt = time;
dt(1) = [];
vel = diff(pos);

% center pre-stimulus velocity around zero
pre = vel(1:numBefore);
mm = mean(pre);
vel = vel - mm;

figure();
plot(time,pos,'-k');
hold on
xl = xline(time(numBefore + 1),'-k','Stimulus arrival');
xl.LabelVerticalAlignment = 'bottom';
xl.LabelHorizontalAlignment = 'left';
title('Average eyelid position across trials in an acquisition', 'FontSize', 16);
xlabel(['Time [msec]'], 'FontSize', 16);
ylabel(['Eyelid position [px]'], 'FontSize', 16);
ylim([-300 100]);
yyaxis right
plot(dt,vel,'-r');
ylabel(['Velocity [px/msec]'], 'FontSize', 16);
ylim([-30 10]);
xmin = time(1);
xmax = time(1) + 1250;
xticks([0 250 500 750 1000])
xlim([xmin xmax]);
hold off

%% plot time series for each PSI
figure();

for pp = 1:5
    ds = tabularTextDatastore(location,'FileExtensions',[".csv"]);
    ds.ReadSize = 'file';
    x = [];
    y = [];
    pos = [];
    
    % loop through scan files
    for ii = 1:52
       % make table from file
       [data,info] = read(ds);
       T = data;
   [~,subStr] = fileparts(info.Filename);
   subStr = split(subStr,'_'); subStr = subStr{end};
   scan = str2num(subStr);
       allVarNames = T.Properties.VariableNames;

       if ismember(scan,psiArray(:,pp))
           % find stimulus arrivals
           [rights,col,v] = find(strcmp('MC-OD',T.Stimulus(:,1)));
           [lefts,col,v] = find(strcmp('MC-OS',T.Stimulus(:,1)));
           all = sort(cat(1,rights,lefts));

           % get starting and ending rows
           starts = all - numBefore;
           ends = all + numAfter;
           if starts(1) < 1
               all(1) = [];
               starts(1) = [];
               ends(1) = [];
           end
           time = table2array(T(starts(1):ends(1),1));
           time = time - time(numBefore + 1);

           % get means across trials
           for jj = 1:length(starts)
               temp = T(starts(jj):ends(jj),:);
               if ismember(all(jj),rights)
                   pos(:,end+1) = table2array(temp(:,3))';
               else
                   pos(:,end+1) = table2array(temp(:,4))';
               end
           end
       end
    end
    
    % get mean across acquisitions
    pos = mean(pos,2);
    
    % center pre-stimulus around zero
    pre = pos(1:numBefore);
    mm = mean(pre);
    pos = pos - mm;
    
    % plot time series
    plot(time,pos,'-k');
    hold on
    if pp == 1
        xl = xline(time(numBefore + 1),'-k','Stimulus arrival');
        xl.LabelVerticalAlignment = 'bottom';
        xl.LabelHorizontalAlignment = 'left';        
    end
    title(['Average eyelid position across trials'], 'FontSize', 16);
    xlabel(['Time [msec]'], 'FontSize', 16);
    ylabel(['Eyelid position [px]'], 'FontSize', 16);
    xmin = time(1);
    xmax = time(1) + 1250;
    xlim([xmin xmax]);
    xticks([0 250 500 750 1000])
    ylim([-300 100]);
    legend('3.5 PSI', '', '7.5 PSI', '15 PSI', '30 PSI', '60 PSI', 'Location', 'southeast');
    colormap(gray);
    
end