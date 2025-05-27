function runCompareCtrl(strWindType,loadData,figNo1,useFASTForComparison,figDirStr)
% runComparCtrl compares simulation results for baseline and qLPV MPC
% controller in closed loop with the Simulink WECS model.
%
% All inputs are optional:
% - strWindType: Two testcases: step sweep 4 to 25 ms and normal dist. with 
%   18 m/s mean (Default: Sweep)
% - loadData: load simulation output data if available instead of running
%   simulation (Default: 1)
% - figNo1: Number of figure (Default: 1)
% - useFASTForComparison: Uses FAST(1) or Simulink data for comparison

%% Handle optional inputs
% The default inputs are provided here.

if ~nargin || isempty(strWindType) 
    strWindType =  'NTW18';  % 'Sweep'; %
end

if nargin < 2 || isempty(loadData)
    loadData = 0; %load simulation output data if available;
end

if nargin < 3 || isempty(figNo1) 
    figNo1 = 1;
end

if nargin < 4 || isempty(useFASTForComparison)
    useFASTForComparison = 1;
end

if nargin <5 || isempty(figDirStr)
    figDirStr = 'figDir';
end

%% Set path to directories

% Set path to inputdata directory
workDir = fileparts(mfilename('fullpath'));
mainDir = fileparts(workDir);
dataInDir = fullfile( mainDir,'dataIn');
addpath(dataInDir);

% Set path to figure directory
figDir = fullfile(mainDir,figDirStr);
if ~isfolder(figDir)
    mkdir(figDir)
end

% Set path to data output directory
dataDirOut = fullfile(workDir, 'dataOut');
if ~isfolder(dataDirOut)
    mkdir(dataDirOut)
end

% Names of simulation output in correct order
varnames = {'Wind', 'RotSpeed', 'GenPwr', 'GenTq', 'BlPitch1', ...
    'NcIMUTAxs', 'NcIMUTAys'}; 

% Weight values (copied from MATLAB function)
q_ = [1 10^4 0 10^3 10^3 0 0 0];
r_ = [1 10^4]; 
p = 10^3;

% Switch between Sweep and NTW18
if strcmp(strWindType,'Sweep') == 1 % sweep from 4 to 25 in steps
    outDataSimulationMat = 'OutDataSweep.mat';
    strFig = ''; 
else % wind with average 18 m/s
    outDataSimulationMat = 'OutDataWind18NTW.mat';
    strFig = 'NTW18'; 
end

%% Load baseline and MPC data or run simulations

%For FASTtool simulation data
load(outDataSimulationMat ,'OutTable'); 

% Load data closed loop simulation PI
if useFASTForComparison % from FAST
    OutTableTest2 = OutTable;
else % from Simulink
    loadBaseline = 1; %always load if available
    matFileOutTableTest2 = fullfile(dataDirOut,['OutTableTest2',strFig,'.mat']);
    simMdlname2 = 'test_SimulinkMdl2_Baseline';
    OutTableTest2 = getSimulationOutputTable(matFileOutTableTest2,loadBaseline,OutTable,simMdlname2);    
end

% Check that SI units are used (can be removed, 
if mean(OutTableTest2.GenTq) <100 % protection against legacy data in kNm
    OutTableTest2.GenTq = OutTableTest2.GenTq*1000; % kNm -> Nm
end

if mean(OutTableTest2.RotSpeed) >1 
    OutTableTest2.RotSpeed = OutTableTest2.RotSpeed *pi/30; %RPM -> rad/s
end

% Run or load qLPV controller
simMdlname = 'test_SimulinkMdl2_qLPVMPCbeta.slx'; %
matFileOutTableTest1 = fullfile(dataDirOut,['OutTableMPC',strFig,'.mat']);
[OutTableMPC, tictoc_LPVMPC,GenPwrRef] = getSimulationOutputTable(matFileOutTableTest1,loadData,OutTable,simMdlname,varnames);

%% Prepare information for plots and plot data

% Information for plot
figStr = ['RefSimulinkClosedLoop',strFig,'_beta'];
meanOpt = strrep(sprintf(' mean cputime MPC: %1.2e s',mean(tictoc_LPVMPC.Data)),'e-0','e-');
timeForPlot = [15,400];
DT = 0.008;
r = GenPwrRef;

% Plot and compare closed-loop results
plotComparison(OutTableTest2,OutTableMPC,q_,r_,p,DT,figStr,figDir,r,timeForPlot,meanOpt,figNo1)

function plotComparison(OutTableTest2,OutTableMPC,q_,r_,p,DT,figStr,figDir,r,timeForPlot,meanOpt,figNo1)
%plotComparison plots the time series obtained with baseline and qLPVMPC

%% Handle optional inputs
maxTime = min(height(OutTableTest2),height(OutTableMPC));
if nargin <10
    timeForPlot(1) = 0;
    timeForPlot(2) = maxTime;
end

%% Vector for plot: Time index for plot
timeVec = 0:DT: maxTime*DT - DT;

idxPlot = timeVec >= timeForPlot(1) & timeVec <= timeForPlot(2);
timeVecPlot = timeVec(idxPlot);
OutTableTest2Plot = OutTableTest2(idxPlot,:);
OutTableMPCPlot = OutTableMPC(idxPlot,:);

PGRef = r.Data(idxPlot);

%% Title strings 
% Display optimization weights and calculate decrease in tower movement and
% power variance.

% String of optimization weights
qVec = deblank(sprintf('%d ',ceil(q_)));
rVec = deblank(sprintf('%d ', ceil(r_)));
pVec = deblank(sprintf('%d', p(1)/q_(1)));

% Improvement tower foreaft acceleration
varPI = var(OutTableTest2Plot.NcIMUTAxs);
varqLPV = var(OutTableMPCPlot.NcIMUTAxs);
titleStdPI = sprintf('%2.2e',varPI);
titleStdqLPV = sprintf('%2.2e',varqLPV);
ratioStd = sprintf('%2.1f ', varPI/varqLPV);

varPI_Pwr = var(abs(PGRef - OutTableTest2Plot.GenPwr));
varqLPV_Pwr = var(abs(PGRef - OutTableMPCPlot.GenPwr));
titleStdPI_Pwr = sprintf('%2.2e',varPI_Pwr);
titleStdqLPV_Pwr = sprintf('%2.2e',varqLPV_Pwr);
ratioStd_Pwr = sprintf('%2.1f ', varPI_Pwr/varqLPV_Pwr);

% axis for plot
yAxCell = {'wind V (m/s)', 'Twr_{FA} (m/s^2)','GenPwr P_g (kW)'};

%% Plot
% Plot input wind, tower fore-aft acceleration and generator power

figure(figNo1);

% Plot wind
ax1(1) = subplot(3,1,1);
plot(timeVecPlot, OutTableMPCPlot.Wind);
ylabel(yAxCell{1}); 
axis tight; grid on;
title({['Baseline (P-PI) vs. qLMPC; ', meanOpt],... %WT mdl: Model2
    ['q_ = [', qVec,'], r = [',rVec,'], p = ',pVec,'* q']})

% Tower fore-aft acceleration
ax1(2) = subplot(3,1,2);
plot(timeVecPlot, OutTableTest2Plot.NcIMUTAxs, timeVecPlot,OutTableMPCPlot.NcIMUTAxs,'-.');
title([yAxCell{2}, ': var_{PI}: ',titleStdPI , ', var_{qLMPC}: ' ,titleStdqLPV,', ratio: ' ,ratioStd])
ylabel(yAxCell{2}); %'y_t [m/s^2]')
axis tight; grid on;

% Generator power
ax1(3) = subplot(3,1,3);
plot(timeVecPlot, OutTableTest2Plot.GenPwr, timeVecPlot,OutTableMPCPlot.GenPwr,'-.',timeVecPlot,PGRef,'k:');
title(['|P_{g,ref}-P_g| [kW]: var_{PI}: ',titleStdPI_Pwr, ', var_{qLMPC}: ' ,...
    titleStdqLPV_Pwr,', ratio: ' ,ratioStd_Pwr])
ylabel(yAxCell{3}); %'P_g [kW]')
axis tight; grid on;
legend('P-PI','MPC','P_{g,ref}','Location','SouthEast')
xlabel('Time (s)');

% Link axes and set poisition
linkaxes(ax1,'x')
set(gcf,'Name', figStr)

posDefault = [520   378   560   420]; %get(gcf, 'position');
set(gcf, 'position', [posDefault(1:3),posDefault(4)*1.1]);

set(findall(gcf,'-property','FontSize'),'FontSize',11.5)
set(findall(gcf,'-property','LineWidth'),'LineWidth',0.75)


print(gcf,[fullfile(figDir,'cmpCtrlSimulink_PI_MPC'),'_',figStr,'_',num2str(timeVecPlot(end))], '-dpng');
print(gcf,[fullfile(figDir,'cmpCtrlSimulink_PI_MPC'),'_',figStr,'_',num2str(timeVecPlot(end))], '-depsc');


%% Figure for paper and dissertation
figure(figNo1+1)

% dc = [double('x'),776]; chX = char(dc);
%    yAxCell = {'Wind V (m/s)', 'GenTq T_g (kNm)', 'Pitch \beta (°)', 'RotSpd \omega_r (rpm)',...
%     'GenPwr P_g (MW)',['Twr_{FA}',chX,'_{t} (m/s^2)'], 'Twr_{SW} x_t (m/s^2)'};

   yAxCell = {'Wind V (m/s)', 'GenTq T_g (kNm)', 'Pitch \beta (°)', 'RotSpd \omega_r (rpm)',...
    'GenPwr P_g (MW)','TwrAcc_{FA} (m/s^2)', 'TwrAcc_{SW} (m/s^2)'};



titleStr = ['Baseline (P-PI) vs. qLPV MPC; ', meanOpt];
nAx = 7;
tiledlayout(nAx,1,'TileSpacing','Compact');

axPlotAll(1) = nexttile; %subplot(nAx,1,1);
plot(timeVecPlot,OutTableMPCPlot.Wind); 
axis tight; grid on;
ylabel(yAxCell{1}); %'wind [m/s]')
title(titleStr);
  
axPlotAll(2) = nexttile; %axPlotAll(2) = subplot(nAx,1,2);  
plot(timeVecPlot,OutTableTest2Plot.GenTq/1000,timeVecPlot,OutTableMPCPlot.GenTq/1000,'-.');
axis tight;
posAxis = axis;
axis([posAxis(1:2), min(43,posAxis(3)), 44]);
grid on; % axis tight;
%grid on;
ylabel(yAxCell{2}); %'GenTq T_g [Nm]') 

axPlotAll(3) = nexttile; %axPlotAll(3) = subplot(nAx,1,3);
plot(timeVecPlot,OutTableTest2Plot.BlPitch1,timeVecPlot,OutTableMPCPlot.BlPitch1,'-.'); 
axis tight; grid on;
ylabel(yAxCell{3}); %'Pitch \beta [°]')

axPlotAll(4) = nexttile; %axPlotAll(4) = subplot(nAx,1,4);
plot(timeVecPlot,OutTableTest2Plot.RotSpeed * 30/pi,timeVecPlot,OutTableMPCPlot.RotSpeed,'-.');
axis tight; grid on;
ylabel(yAxCell{4}); %'RotSpd \omega_r [rpm]')

axPlotAll(5) = nexttile; %axPlotAll(5) = subplot(nAx,1,5);
plot(timeVecPlot,OutTableTest2Plot.NcIMUTAxs,timeVecPlot,OutTableMPCPlot.NcIMUTAxs,'-.');
axis tight; grid on;
ylabel(yAxCell{6}); %'Twr_{FA} y_t[m/s^2]')

axPlotAll(6) = nexttile; %axPlotAll(6) = subplot(nAx,1,6);
plot(timeVecPlot,OutTableTest2Plot.NcIMUTAys,timeVecPlot,OutTableMPCPlot.NcIMUTAys,'-.');
axis tight; grid on;
ylabel(yAxCell{7}); %'Twr_{SW} x_t[m/s^2]')
set(axPlotAll(6),'YLim', get(axPlotAll(5) ,'YLim'))

axPlotAll(7) = nexttile; %axPlotAll(7) = subplot(nAx,1,7);
plot(timeVecPlot,OutTableTest2Plot.GenPwr/1000,timeVecPlot,OutTableMPCPlot.GenPwr/1000,'-.');
axis tight; grid on;
ylabel(yAxCell{5}); %'GenPwr P_g [MW]')
xlabel('Time (s)')
linkaxes(axPlotAll,'x');

set(gcf,'Name',['cmpCtrlTimeDomain_All_',figStr])
posDefault = get(0,'DefaultFigurePosition');
set(gcf, 'position', [posDefault(1),posDefault(2) - posDefault(4)*0.7,posDefault(3),posDefault(4)*2.1]);

set(findall(gcf,'-property','FontSize'),'FontSize',11.5)
set(findall(gcf,'-property','LineWidth'),'LineWidth',0.75)

%set(groot,'defaultLineLineWidth',defaultLineWidth);

print(fullfile(figDir,['cmpCtrlTimeDomain_All_',figStr]), '-dpng');
print(fullfile(figDir,['cmpCtrlTimeDomain_All_',figStr]), '-depsc');

