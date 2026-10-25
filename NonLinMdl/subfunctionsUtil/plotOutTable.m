function axPlotAll = plotOutTable(OutTable,OutTableTest1,OutTableTest2,testCaseStr, testCaseCell, figNo2,figDir,strFig,timeEval,OutTableTest3,plotTest3Lines)
% timeEval (optional): time window [s] of the plots and the statistics
% (Default: plots from 20 s, statistics from 15 s to the end)
% OutTableTest3 (optional): fourth case, label testCaseCell{4}; additional
% line (color 4 of the color map) and last column of the statistics tables
% plotTest3Lines (optional, default true): false writes the fourth case only
% to the statistics tables, not to the plot


%% Set parameters and missing inputs
% yAxCell = {'Wind V [m/s]', 'GenTq T_g [kNm]', 'Pitch \beta [°]', 'RotSpd \omega_r [rpm]',...
%     'GenPwr P_g [MW]','Twr_{FA} y_t [m/s^2]'};
yAxCell = {'Wind (m/s)', 'GenTq (kNm)', 'Pitch (°)', 'RotSpd (rpm)',...
    'GenPwr (MW)','TwrAcc_{FA} (m/s^2)'};

if nargin < 5 || isempty(testCaseCell)
    testCaseCell = {'Mdl1: Rot+Twr ', 'Mdl2: Rot,Twr,Bld+Act ','FAST'};
end

if nargin < 6 || isempty(figNo2)
    figNo2 = 1;
end

if nargin<7 || isempty(figDir)
    % Set path  of default figure directory
    workDir = fileparts(mfilename('fullpath'));
    mainDir = fileparts(workDir);
    figDir = fullfile(mainDir,'figDir');
end
if ~isfolder(figDir)
    mkdir(figDir)
end

if nargin<8 || isempty(strFig)
    strFig = 'testConstr';
end

if nargin<9 || isempty(timeEval)
    timeEval = [];
end

if nargin<10
    OutTableTest3 = [];
end
hasTest3 = ~isempty(OutTableTest3);
if nargin < 11 || isempty(plotTest3Lines)
    plotTest3Lines = true;
end
showTest3 = hasTest3 && plotTest3Lines; % fourth case in the plot
% Assign data from FAST to variables

%% Handle different variable names
idxTime =1: min([height(OutTable), height(OutTableTest1),height(OutTableTest2)]);
if hasTest3
    idxTime = 1:min(numel(idxTime), height(OutTableTest3));
end
OutTable = OutTable(idxTime,:);
OutTable.BldPitch1 = OutTable.BlPitch1;
OutTable.Torque = OutTable.GenTq;

% Calculate and assign wind amplitude from comonents
idxWind = contains(OutTable.Properties.VariableNames, 'Wind');
vectWind = OutTable{:,idxWind};
vectAmpWind = sqrt(sum((vectWind.^2),2));

% Get color map for 'title legends'
cl = colormap('lines');
% titleStr2 = [testCaseStr,': ',testCaseCell{3},' {\color[rgb]{',num2str(cl(1,:)),'}',testCaseCell{3},' ',...
%     '\color[rgb]{',num2str(cl(2,:)),'}',testCaseCell{1}];
% titleStr = [testCaseStr,': Mdl1: Rot+Twr {\color[rgb]{',num2str(cl(1,:)),'}Mdl2: Rot,Twr,Bld+Act ',...
%     '\color[rgb]{',num2str(cl(2,:)),'}FAST} '];
titleStr = [testCaseStr,': ',testCaseCell{1},'{\color[rgb]{',num2str(cl(1,:)),'}',testCaseCell{2},...
    '\color[rgb]{',num2str(cl(2,:)),'}',testCaseCell{3},'}'];
if showTest3
    titleStr = [titleStr(1:end-1), ' \color[rgb]{',num2str(cl(4,:)),'}',testCaseCell{4},'}'];
end

% Create output plot time reference
time = OutTable.Time(idxTime);
if isempty(timeEval)
    idxTime = time > 20;
else
    idxTime = time >= timeEval(1) & time <= timeEval(2);
end
time = time(idxTime);

% Check if generator torque has meaningful variation
plotGenTq = std(OutTableTest2.GenTq(idxTime)/10^3) > 1 || ...
            std(OutTable.GenTq(idxTime))            > 1 || ...
            std(OutTableTest1.GenTq(idxTime)/10^3)  > 1;

% Determine number of subplots
nPlots = 5 + plotGenTq;

figure(figNo2)
t = tiledlayout(nPlots, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

% --- Wind speed
axPlotAll(1) = nexttile;
plot(time, vectAmpWind(idxTime), time, OutTableTest1.Wind(idxTime), '--');
axis tight; grid on;
ylabel(yAxCell{1});
title(titleStr);

% --- Generator torque (conditional)
if plotGenTq
    axPlotAll(2) = nexttile;
    plot(time, OutTableTest2.GenTq(idxTime)/10^3, ...
         time, OutTable.GenTq(idxTime), ...
         time, OutTableTest1.GenTq(idxTime)/10^3, 'k--');
    plotTest3(showTest3, time, OutTableTest3, 'GenTq', idxTime, 1/10^3, cl(4,:));
    axis tight;
    posAxis = axis;
    axis([posAxis(1:2), min(43, posAxis(3)), 44]);
    grid on;
    ylabel(yAxCell{2});
end

% --- Blade pitch
axPlotAll(3) = nexttile;
plot(time, OutTableTest2.BlPitch1(idxTime), ...
     time, OutTable.BlPitch1(idxTime), ...
     time, OutTableTest1.BlPitch1(idxTime), 'k--');
plotTest3(showTest3, time, OutTableTest3, 'BlPitch1', idxTime, 1, cl(4,:));
axis tight; grid on;
ylabel(yAxCell{3});

% --- Rotor speed
axPlotAll(4) = nexttile;
plot(time, OutTableTest2.RotSpeed(idxTime), ...
     time, OutTable.RotSpeed(idxTime), ...
     time, OutTableTest1.RotSpeed(idxTime), 'k--');
plotTest3(showTest3, time, OutTableTest3, 'RotSpeed', idxTime, 1, cl(4,:));
axis tight; grid on;
ylabel(yAxCell{4});

% --- Tower fore-aft acceleration
axPlotAll(5) = nexttile;
plot(time, OutTableTest2.NcIMUTAxs(idxTime), ...
     time, OutTable.NcIMUTAxs(idxTime), ...
     time, OutTableTest1.NcIMUTAxs(idxTime), 'k--');
plotTest3(showTest3, time, OutTableTest3, 'NcIMUTAxs', idxTime, 1, cl(4,:));
axis tight; grid on;
ylabel(yAxCell{6});

% --- Generator power
axPlotAll(6) = nexttile;
plot(time, OutTableTest2.GenPwr(idxTime)/1000, ...
     time, OutTable.GenPwr(idxTime)/1000, ...
     time, OutTableTest1.GenPwr(idxTime)/1000, 'k--');
plotTest3(showTest3, time, OutTableTest3, 'GenPwr', idxTime, 1/1000, cl(4,:));
axis tight; grid on;
ylabel(yAxCell{5});
xlabel('Time (s)');

linkaxes(axPlotAll(isgraphics(axPlotAll, 'axes')), 'x');

set(gcf,'Name',['cmpTimeDomain_All',strFig])
posDefault = [520   378   560   420]; %get(gcf, 'position');
set(gcf, 'position', [posDefault(1:3),posDefault(4)*1.7]);
% set(groot,'defaultLineLineWidth',defaultLineWidth);

set(findall(gcf,'-property','FontSize'),'FontSize',11.5)
set(findall(gcf,'-property','LineWidth'),'LineWidth',0.75)

% 
print(fullfile(figDir,['cmpTimeDomain_All',strFig]), '-dpng');
print(fullfile(figDir,['cmpTimeDomain_All',strFig]), '-depsc');


%% Calculate variance and standard deviation for key signals
% Use the same time index as in plotOutTable (time > 15s)
time    = OutTable.Time;
if isempty(timeEval)
    idxTime = time > 15;
else
    idxTime = time >= timeEval(1) & time <= timeEval(2);
end

% ---- Extract signals (matching unit conversions used in plotOutTable) ----
dt = mean(diff(OutTable.Time(idxTime)));   % [s]

% ---- Delta Pitch: pitch rate in °/s ----
dPitch_1   = diff(OutTable.BlPitch1(idxTime))     / dt;
dPitch_2   = diff(OutTableTest2.BlPitch1(idxTime)) / dt;
dPitch_Ref = diff(OutTableTest1.BlPitch1(idxTime)) / dt;

% Rotor speed in rpm (the Simulink test models log rpm already)
rotSpd_1   = toRpm(OutTable.RotSpeed(idxTime));
rotSpd_2   = toRpm(OutTableTest2.RotSpeed(idxTime));
rotSpd_Ref = toRpm(OutTableTest1.RotSpeed(idxTime));

twrFA_1    = OutTable.NcIMUTAxs(idxTime);
twrFA_2    = OutTableTest2.NcIMUTAxs(idxTime);
twrFA_Ref  = OutTableTest1.NcIMUTAxs(idxTime);

genPwr_1   = OutTable.GenPwr(idxTime);
genPwr_2   = OutTableTest2.GenPwr(idxTime);
genPwr_Ref = OutTableTest1.GenPwr(idxTime);

% ---- Signal labels and data: order matches desired table output ----
criteriaNames = {'TwrFA (m/s²)', 'GenPwr (kW)', 'RotSpeed (rpm)', 'Pitch rate (°/s)'};
modelNames    = testCaseCell;

data = {
    twrFA_1,    twrFA_2,    twrFA_Ref;
    genPwr_1,   genPwr_2,   genPwr_Ref;
    rotSpd_1,   rotSpd_2,   rotSpd_Ref;
    dPitch_1,   dPitch_2,   dPitch_Ref;
};
if hasTest3 % fourth case as last column
    data(:,4) = {OutTableTest3.NcIMUTAxs(idxTime); OutTableTest3.GenPwr(idxTime); ...
        toRpm(OutTableTest3.RotSpeed(idxTime)); diff(OutTableTest3.BlPitch1(idxTime)) / dt};
end
nModels = size(data, 2);
modelNames = modelNames(1:nModels);

% ---- Compute variance and std dev matrices ----
nCriteria = numel(criteriaNames);

varMat = zeros(nCriteria, nModels);
stdMat = zeros(nCriteria, nModels);

for c = 1:nCriteria
    for m = 1:nModels
        varMat(c,m) = var(data{c,m});
        stdMat(c,m) = std(data{c,m});
    end
end

% ---- Build tables ----
varianceTable = [table(criteriaNames(:), 'VariableNames', {'Criterion'}), ...
    array2table(varMat, 'VariableNames', modelNames)];

stdDevTable = [table(criteriaNames(:), 'VariableNames', {'Criterion'}), ...
    array2table(stdMat, 'VariableNames', modelNames)];

% ---- Swap columns (Ref=Test1 is col3, swap to col1 position) ----
tmpNames = varianceTable.Properties.VariableNames;
tmpNames([2,4]) = tmpNames([4,2]);
varianceTable.Properties.VariableNames = tmpNames;

tmpNames = stdDevTable.Properties.VariableNames;
tmpNames([2,4]) = tmpNames([4,2]);
stdDevTable.Properties.VariableNames = tmpNames;

disp('=== Variance ===');       disp(varianceTable);
disp('=== Standard Deviation ==='); disp(stdDevTable);

% ---- Write LaTeX ----
fid = fopen(fullfile(figDir, ['statsTable',strFig,'.tex']), 'w');
writeLatexTable(fid, varianceTable, ...
    'Variance of key signals for different pitch rate constraints', ...
    'tab:variance');
writeLatexTable(fid, stdDevTable, ...
    'Standard deviation of key signals for different pitch rate constraints', ...
    'tab:stddev');
fclose(fid);
fprintf('LaTeX tables written to %s\n', fullfile(figDir, ['statsTable',strFig,'.tex']));
% Unrounded standard deviations for the ratios of merge_stats_tables
save(fullfile(figDir, ['statsTable',strFig,'.mat']), 'stdDevTable');
end

% Helper to write one table
function writeLatexTable(fid, tbl, caption, label)
    colNames = tbl.Properties.VariableNames;
    fprintf(fid, '\\begin{table}[htbp]\n');
    fprintf(fid, '\\centering\n');
    fprintf(fid, '\\caption{%s}\n', caption);
    fprintf(fid, '\\label{%s}\n', label);
    nCol = numel(colNames) - 1;
    fprintf(fid, '\\begin{tabular}{l%s}\n', repmat('r', 1, nCol));
    fprintf(fid, '\\toprule\n');
    fprintf(fid, '%s \\\\\n', strjoin(colNames, ' & '));
    fprintf(fid, '\\midrule\n');
    % Row labels matching new order: TwrFA, GenPwr, RotSpeed, Pitch rate.
    % Above rated the power reference is constant, so the std of the
    % generator power equals the std of the power tracking error (runCompareCtrl).
    % Number formats as in runCompareCtrl: tower acceleration and rotor
    % speed with four digits, the other signals with two digits
    criteria = {
        'Tower fore-aft acc.\ (m/s$^2$)',   '%.4f';
        'Power tracking error (kW)',        '%.2f';
        'Rotor speed (rpm)',                '%.4f';
        'Blade pitch rate (\textdegree/s)', '%.2f'
    };
    for r = 1:height(tbl)
        fmt = criteria{r,2};
        fprintf(fid, ['%s', repmat([' & ', fmt], 1, nCol), ' \\\\\n'], ...
            criteria{r,1}, tbl{r,2:end});
    end
    fprintf(fid, '\\bottomrule\n');
    fprintf(fid, '\\end{tabular}\n');
    fprintf(fid, '\\end{table}\n\n');
end

% Helper: rotor speed in rpm (signals in rad/s are converted)
function rpm = toRpm(rotSpd)
    if mean(rotSpd) < 3 % rad/s (rated 1.27 rad/s = 12.1 rpm)
        rpm = rotSpd * 60/(2*pi);
    else
        rpm = rotSpd;
    end
end

% Helper: line of the fourth case (if available) in the current axes
function plotTest3(hasTest3, time, OutTableTest3, sigName, idxTime, scale, color)
    if ~hasTest3
        return
    end
    hold on;
    plot(time, OutTableTest3.(sigName)(idxTime)*scale, 'Color', color);
    hold off;
end
