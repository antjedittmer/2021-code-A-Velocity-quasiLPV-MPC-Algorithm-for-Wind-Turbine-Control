function plotTimeseriesCtrlComparison(idxTime, outTableNRELorig, outTableNRELnew,outTableDLR, strNewCtrl,figNo1,figNo2,testCaseStr)
% plotTimeseriesCtrlComparison plots the timeseries of the NREL Ctrl shipped
% with FAST and OC3 against a new controller design
%
% Inputs:
% 

%% Prepare Inputs

% time selection
if nargin <1 || isempty(idxTime)
    idxTime = outTableNRELorig.Time >= 15;   
end

if nargin < 5 
    strNewCtrl = '';
end

if nargin < 6 
    figNo1 = 1; 
    figNo2 = 2;
end

if nargin < 67
    figNo1 = 1; 
    figNo2 = figNo1 +1;
end

if nargin <= 8
    testCaseStr = 'NTW 18';
end
%% Prepare data for plots


% Get time basis for plots
time = outTableNRELorig.Time(idxTime);

% Calculate wind amplitude
cellChannelWind = {'Wind'};
OutList = outTableNRELorig.Properties.VariableNames;
idxWind = contains(OutList, cellChannelWind);
vectWind = outTableNRELorig{:,idxWind};
vectAmpWind = sqrt(sum((vectWind.^2),2));

vectWind = outTableNRELnew{:,idxWind};
vectAmpWindNREL = sqrt(sum((vectWind.^2),2));

vectWind = outTableDLR{:,idxWind};
vectAmpWindDLR = sqrt(sum((vectWind.^2),2));

% Get color map for 'title legends'
cl = colormap('lines');
titleStr = ['Test Case ',testCaseStr,': {\color[rgb]{',num2str(cl(1,:)),'}PI-CPC ',...
    '\color[rgb]{',num2str(cl(2,:)),'}PI-IPC_{fminsearch} } ', strNewCtrl];

%% 
defaultLineWidth = get(groot,'defaultLineLineWidth');
set(groot,'defaultLineLineWidth',1);
pitchName = OutList{contains(OutList,'Pitch')};

figure(figNo1);
axPlot(1) = subplot(3,1,1);
plot(time,vectAmpWind(idxTime),time,vectAmpWindNREL(idxTime),'--',time,vectAmpWindNREL(idxTime),'k:'); 
axis tight; grid on;
ylabel('Wind Amplitude (m/s)')
title(titleStr);

axPlot(2) = subplot(3,1,2);
plot(time,outTableNRELorig.RotSpeed(idxTime),time,...
    outTableNRELnew.RotSpeed(idxTime),'--',time,outTableDLR.RotSpeed(idxTime),'k-.');
axis tight; grid on;
ylabel('RotSpeed (rpm)')

axPlot(3) = subplot(3,1,3);
plot(time,outTableNRELorig.GenPwr(idxTime)/1000, ....
    time,outTableNRELnew.GenPwr(idxTime)/1000,'--',time,outTableDLR.GenPwr(idxTime)/1000,'k-.');
axis tight; grid on;
ylabel('GenPwr (MW)')
xlabel('Time (s)')
linkaxes(axPlot,'x');


figure(figNo2)
axPlot2(1) = subplot(4,1,1);
plot(time,outTableNRELorig.(pitchName)(idxTime), time,outTableNRELnew.(pitchName)(idxTime),'--',...
    time,outTableDLR.(pitchName)(idxTime),'k-.'); 
axis tight; grid on;
ylabel('BldPitch1 (rad)')
title(titleStr);

axPlot2(2) = subplot(4,1,2);
plot(time,outTableNRELorig.GenTq(idxTime),time,outTableNRELnew.GenTq(idxTime),'--',...
    time,outTableNRELnew.GenTq(idxTime),'k-.' );
axis tight; grid on;
ylabel('GenTq (Nm)') 

axPlot2(3) = subplot(4,1,3);
plot(time,outTableNRELorig.NcIMUTAxs(idxTime),time,outTableNRELnew.NcIMUTAxs(idxTime),'--',...
    time,outTableDLR.NcIMUTAxs(idxTime),'k-.');
axis tight; grid on;
ylabel('Tower FA (m/s)')

axPlot2(4) = subplot(4,1,4);
plot(time,outTableNRELorig.NcIMUTAys(idxTime),time,outTableNRELnew.NcIMUTAys(idxTime),'--',...
time,outTableDLR.NcIMUTAxs(idxTime),'k-.');
axis tight; grid on;
ylabel('Tower SS (m/s)')
xlabel('Time (s)')
linkaxes(axPlot2,'x');
set(groot,'defaultLineLineWidth',defaultLineWidth);

