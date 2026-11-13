% Script that can be used to call bar plots

clc; clear; close all;
addpath(genpath(fileparts(pwd)));

dirFigAll ='dirFig';
if ~isdir(dirFigAll) %#ok<ISDIR>
    mkdir(dirFig);
end

dirFigDiss ='dirFigDiss';
if ~isdir(dirFigDiss) %#ok<ISDIR>
    mkdir(dirFigDiss);
end
onlyDissPics = 1;
fs = 12;

dirFig = dirFigDiss;
%% Parameter for DEL
% C                               = 2/sqrt(pi)*10;% [m/s] TC I
% k                               = 2;            % [-]
WoehlerExponentSt               = 4;            % [-]   for steel
WoehlerExponentCFP              = 10;           % [-]   for carbon fiber-reinforced polymers
% N_REF                           = 2e6/(20*8760);% [-]   fraction of 2e6 in 20 years for 1h
timePlot = 100;

bStrAll = 'MIP'; % in-plane moments (edgewise rotor coordinate system)
cStrAll = 'MOoP'; % out of plane moments (flapwise rotor coordinate system)
mzStrAll = 'Mzb';
tStrAll = 'TwrBsM';

pitchStrAll = 'BlPitch';
blStr1 = sprintf('%s%d',pitchStrAll,1);
blStr2 = sprintf('%s%d',pitchStrAll,2);
blStr3 = sprintf('%s%d',pitchStrAll,3);

allMatlabFiles = what;
allMatFiles = allMatlabFiles.mat;
idxIPCFiles = contains(allMatFiles,'IPC1');

% cellIPCNames = allMatFiles(idxIPCFiles);
cellIPCNames =  {'OutTable_IPC1_MaxRate4_Kp0dot00.mat',...
    'OutTable_IPC1_MaxRate8_Kp0dot04.mat', 'OutTable_IPC1_MaxRate8_Kp0dot71.mat',  ... %'OutTable_IPC1_MaxRate8_Kp0dot00.mat'... 'OutTable_IPC1_MaxRate8_Kp0dot35.mat', ... %
    'OutTable_IPC1_MaxRate13_Kp0dot71.mat','OutTable_IPC1_MaxRate13_Kp1dot78.mat'};

%vec1 = [1:2,4:5]; actLim = 5;
vec1 = [1,3,5]; actLim = 10;

cellIPCNames = cellIPCNames(vec1);

%% Load one table to get variable names and time index
load(cellIPCNames{1},'OutTable0'); % load first IPC file to get headers and time range

idxTime = OutTable0.Time > timePlot & OutTable0.Time < (timePlot+400);
outVarNames = OutTable0.Properties.VariableNames;
idxBlIP = contains(outVarNames,bStrAll);
idxBlOoP = contains(outVarNames,cStrAll);
idxBlMZb = contains(outVarNames,mzStrAll);
idxBlM =  idxBlIP | idxBlOoP | idxBlMZb;
idxPitch =  contains(outVarNames,pitchStrAll);
idxT = contains(outVarNames,tStrAll);

CPCName = strrep(cellIPCNames{1},'IPC1','IPC0'); %
load(CPCName,'OutTable0');
OutTableCPC = OutTable0(idxTime,:);
Mblade_all_CPC = OutTableCPC(:,idxBlM);
Mblade_pitch_CPC = OutTableCPC(:,idxBlMZb);
Mtower_CPC = OutTableCPC(:,idxT);
Blade_pitch_CPC = OutTableCPC(:,idxPitch);

%% Load IPC
DELvec = struct;
DELvec_RateLim = [4,8,13]';
for idx = 1: length(cellIPCNames)

    load(cellIPCNames{idx},'OutTable0');
    OutTableIPC = OutTable0(idxTime,:); % assign table as IPC table

    figNo = 1;
    strplot = sprintf('IPC limit %d °/s, NTM mean wind 16 m/s',DELvec_RateLim(idx));
    plotFASTresults(OutTableIPC,figNo,strplot,OutTableCPC); figNo = figNo + 1;
    filenamepng = ['fig00a_Timeseries',num2str(idx)];
    filepathpng = fullfile(dirFig,filenamepng);
    print(gcf,filepathpng, '-dpng');
    print(gcf,filepathpng, '-depsc');


    Mblade_all_IPC = OutTableIPC(:,idxBlM);
    % Mblade_OoP_IPC = OutTableIPC{:,idxBlOoP};

    Mblade_pitch_IPC = OutTableIPC(:,idxBlMZb);
    Mtower_IPC = OutTableIPC(:,idxT);
    Blade_pitch_IPC = OutTableIPC(:,idxPitch);

    CPCName = strrep(cellIPCNames{idx},'IPC1','IPC0');
    load(CPCName,'OutTable0');
    OutTableCPC1 = OutTable0(idxTime,:);
    Mblade_all_CPC1 = OutTableCPC1(:,idxBlM);
    Mblade_pitch_CPC1 = OutTableCPC1(:,idxBlMZb);
    Mtower_CPC1 = OutTableCPC1(:,idxT);
    Blade_pitch_CPC1 = OutTableCPC1(:,idxPitch);

    %% Calculate blade DEL ratio between IPC and CPC
    [damageBladesIPC, damageBladesVecIPC] = calculateDamageBladesNew(Mblade_all_IPC{:,:},WoehlerExponentCFP); % MOoP and MIP
    [damageBladesCPC, damageBladesVecCPC] = calculateDamageBladesNew(Mblade_all_CPC{:,:},WoehlerExponentCFP); % Ref CPC at 4 °/s
    [damageBladesCPC1, damageBladesVecCPC1] = calculateDamageBladesNew(Mblade_all_CPC1{:,:},WoehlerExponentCFP); % Ref CPC at rate of IPC
    damagesVec = damageBladesVecIPC./damageBladesVecCPC;
    DELratio_bladeIPC = sum(damageBladesVecIPC)/sum(damageBladesVecCPC);
    damagesVecCPC = damageBladesVecCPC1./damageBladesVecCPC;
    DELratio_bladeCPC = sum(damageBladesVecCPC1)/sum(damageBladesVecCPC);

    %% Calculate tower DEL ratio between IPC and CPC
    [damageTowerIPC, damageTowerVecIPC] = calculateDamageBladesNew(Mtower_IPC{:,:},WoehlerExponentSt);
    [damageTowerCPC, damageTowerVecCPC] = calculateDamageBladesNew(Mtower_CPC{:,:},WoehlerExponentSt);
    [damageTowerCPC1, damageTowerVecCPC1] = calculateDamageBladesNew(Mtower_CPC1{:,:},WoehlerExponentSt);
    DELratio_twrIPC = damageTowerIPC/damageTowerCPC;
    DELratio_twrCPC = damageTowerCPC1/damageTowerCPC;

    %% Calculate abs blade deviation
    absDev = [abs(OutTableIPC.(blStr1) - OutTableIPC.(blStr2)),...
        abs(OutTableIPC.(blStr1) - OutTableIPC.(blStr3)),...
        abs(OutTableIPC.(blStr2)- OutTableIPC.(blStr3))];

    absDiff = [abs(diff(OutTableIPC.(blStr1))),...
        abs(diff(OutTableIPC.(blStr2))),...
        abs(diff(OutTableIPC.(blStr3)))]/0.008;

    %% Calculate power mean and deviation
    meanGenPwrCPC = mean(OutTableCPC.GenPwr);
    stdGenPwrCPC = std(OutTableCPC.GenPwr);
    meanGenPwrCPC1 = mean(OutTableCPC1.GenPwr);
    stdGenPwrCPC1 = std(OutTableCPC1.GenPwr);
    meanGenPwr = mean(OutTableIPC.GenPwr);
    stdGenPwr = std(OutTableIPC.GenPwr);

    %% Calculate actuator activity
    absMomActCPC = abs( diff(Blade_pitch_CPC{:,:}) .* Mblade_pitch_CPC{1:end-1,:});
    meanActPwrCPC = mean( sum(absMomActCPC));

    absMomActCPC1 = abs( diff(Blade_pitch_CPC1{:,:}) .*  Mblade_pitch_CPC1{1:end-1,:});
    meanActPwrCPC1 = mean( sum(  absMomActCPC1));

    absMomActIPC =  abs( diff(Blade_pitch_IPC{:,:}) .* Mblade_pitch_IPC{1:end-1,:});
    meanActPwrIPC = mean( sum(absMomActIPC));

    meanActDiffCPC = mean( sum( abs( diff(Blade_pitch_CPC{:,:}))));
    meanActDiffCPC1 = mean( sum( abs( diff(Blade_pitch_CPC1{:,:}))));
    meanActDiffIPC = mean( sum( abs( diff(Blade_pitch_IPC{:,:}))));


    %% assignValues
    cellNameTC = regexp(cellIPCNames{idx},'_','split');

    DELvec(idx).RateLim = str2double(strrep(strrep(cellNameTC{3},'MaxRate',''),'.mat',''));
    DELvec(idx).KpName = strrep(strrep(cellNameTC{4},'Kp',''),'.mat','');

    DELvec(idx).bladeIPC = (1-DELratio_bladeIPC)*100 ;
    DELvec(idx).twrIPC = (1-DELratio_twrIPC)*100;
    DELvec(idx).actEne = meanActPwrIPC/meanActPwrCPC;
    DELvec(idx).meanGenPwr = meanGenPwr/meanGenPwrCPC;
    DELvec(idx).stdGenPwr = stdGenPwr/stdGenPwrCPC;
    DELvec(idx).bladeDev = max(absDev(:));

    DELvec(idx).bladeDevVec = absDev(:)';
    DELvec(idx).bladeDiffVec = absDiff(:)';

    % Blade moment vectors
    DELvec(idx).bladeMOoPVecIPC = Mblade_all_IPC.RootMOoP1'; %mean(OutTableIPC{:,idxBlOoP},2)';
    DELvec(idx).bladeMIPVecIPC = Mblade_all_IPC.RootMIP1';% mean(OutTableIPC{:,idxBlIP},2)';
    DELvec(idx).bladeMtorqueVecIPC = Mblade_all_IPC.RootMzb1'; %mean(OutTableIPC{:,idxBlMZb},2)';

    % Tower moment vectors
    DELvec(idx).twrBsMxtVecIPC = OutTableIPC.TwrBsMxt';
    DELvec(idx).twrBsMytVecIPC = OutTableIPC.TwrBsMyt';
    DELvec(idx).twrBsMztVecIPC = OutTableIPC.TwrBsMzt';

    % Power vector
    DELvec(idx).powerVecIPC = OutTableCPC.GenPwr';

    % Actuatormoment vector
    DELvec(idx).absMomAcVecIPC = mean(absMomActIPC,2)';

    DELvec(idx).bladeCPC = (1-DELratio_bladeCPC)*100 ;
    DELvec(idx).twrCPC = (1-DELratio_twrCPC)*100;
    DELvec(idx).actEneCPC = meanActPwrCPC1/meanActPwrCPC;
    DELvec(idx).meanGenPwrCPC = meanGenPwrCPC1/meanGenPwrCPC;
    DELvec(idx).stdGenPwrCPC = stdGenPwrCPC1/stdGenPwrCPC;


end

    dirFig = dirFigAll;
    %return;

    %% Write data into table for post-processing
    tableObj0 = struct2table(DELvec);
    tableObj = tableObj0(tableObj0.bladeIPC > 0 & tableObj0.actEne < actLim,:); % &
    tableObjAll = tableObj0(tableObj0.bladeIPC > 0 & tableObj0.actEne > actLim,:); % &

    vecRateLim = unique(tableObj.RateLim);

    %% Plot maximum blade deviation
    close all;
    
    % figNo = 1;
    figure(figNo); figNo = figNo +1;

    nameVec = {'bladeDiffVec','bladeDevVec'};
    ylabelVec = {'Pitch rate, °/s', 'Blade ang. dev., °'};
    tickVec ={{'','',''},{'C_{Bl}: 0.97; C_{Tow}: 1.00','C_{Bl}: 0.87; C_{Tow}: 0.98','C_{Bl}: 0.86; C_{Tow}: 0.96'}};
    %tickVec ={{'','',''},{'C_{Bl}: 0.97','C_{Bl}: 0.87','C_{Bl}: 0.86'}};

    vecBl = (100 - tableObj.bladeIPC)/100;
    vecTwr = (100 - tableObj.twrIPC)/100;

    cl = lines;

    for idx = 1: length(nameVec)
        data = tableObj0.(nameVec{idx})';

        yLimVec = max(data(:));

        ax1(idx) = subplot(2,1,idx); %#ok<SAGROW>
        hold on;
        for i = 1:size(data, 2)
            boxchart(ones(size(data, 1), 1) * i, data(:, i), ...
                'BoxFaceColor', cl(i, :), ...
                'MarkerStyle', '+');

            % Compute statistics
            medianValue = median(data(:, i));
            meanValue = mean(data(:, i));
            maxValue = max(data(:, i));
            maxDisp = max(maxValue,3.25);


            % Annotate maximum value near the max point
            text(i + 0.075, maxDisp, sprintf('Max: %.1f', maxValue), ...
                'HorizontalAlignment', 'left', 'FontSize', fs-1, 'Color', 'r');
            % Annotate maximum value near the max point
            text(i + 0.075, maxDisp - yLimVec/8, sprintf('Med: %.1f', medianValue), ...
                'HorizontalAlignment', 'left', 'FontSize', fs-1, 'Color', 'k');

        end
        hold off;
        axis tight; grid on;
        yLimVec = get(ax1(idx),"YLim");
        set(ax1(idx),"YLim",[yLimVec(1) yLimVec(2)*1.1]);
        ylabel(ylabelVec{idx});

        set(gca,'FontSize',fs-2);
        set(gca,'XTick',1:3,'XTickLabel',tickVec{idx})
        set(gca,'XTickLabelRotation',0)
    end

    legCell = regexp(sprintf('Rate lim. %d°/s,', vecRateLim),',','split');
    legCell0 = legCell(1:end-1);
    legend(legCell0,'Location','northwest')

    pos0 = get(0,'defaultFigurePosition');
    set(gcf,'Position', [pos0(1:3),pos0(4)*0.7])

    filenamepng = 'fig00_CmpBoxPlot';
    filepathpng = fullfile(dirFig,filenamepng);
    print(gcf,filepathpng, '-dpng');
    print(gcf,filepathpng, '-depsc');


    %% Plot the blade root moments
    nameVec = {'bladeMOoPVecIPC','bladeMIPVecIPC'};
    ylabelVec = {'Root MOoP (kNm)', 'Root MIP (kNm)'};
    tickVec ={{'','',''},{'C_{Bl}: 0.97; C_{Tow}: 1.00','C_{Bl}: 0.87; C_{Tow}: 0.98','C_{Bl}: 0.86; C_{Tow}: 0.96'}};

    figNo1 = figNo;
    fcnPlotBoxChart(tableObj,nameVec,vecRateLim,fs,ylabelVec,100 + figNo1);
    figNo1 = figNo1 + 1;

    filenamepng = 'fig00a_CmpBoxPlotBlRoots';
    filepathpng = fullfile(dirFig,filenamepng);
    print(gcf,filepathpng, '-dpng');
    print(gcf,filepathpng, '-depsc');

    %% Plot the tower base moments
    nameVec = {'twrBsMxtVecIPC','twrBsMytVecIPC'};
    ylabelVec = {'Twr Mx (kNm)', 'Twr My (kNm)'};

    fcnPlotBoxChart(tableObj,nameVec,vecRateLim,fs,ylabelVec,100 + figNo1);
    figNo1 = figNo1 + 1;

    filenamepng = 'fig00a_CmpBoxPlotTwrBs';
    filepathpng = fullfile(dirFig,filenamepng);
    print(gcf,filepathpng, '-dpng');
    print(gcf,filepathpng, '-depsc');

    %% Plot the power
    nameVec = {'powerVecIPC'};
    ylabelVec = {'Generator Power (kW)'};

    fcnPlotBoxChart(tableObj,nameVec,vecRateLim,fs,ylabelVec,100 + figNo1);
    figNo1 = figNo1 + 1;

    filenamepng = 'fig00a_CmpBoxPlotPwr';
    filepathpng = fullfile(dirFig,filenamepng);
    print(gcf,filepathpng, '-dpng');
    print(gcf,filepathpng, '-depsc');


    %% Plot the actuator rate and actuator limit
    nameVec = {'bladeDiffVec','absMomAcVecIPC'};
    ylabelVec = {'Pitch rate (°/s)', 'Actuator power (kNm/s)'};
    fcnPlotBoxChart(tableObj,nameVec,vecRateLim,fs,ylabelVec,100 + figNo1);


    filenamepng = 'fig00a_CmpBoxPlotAct';
    filepathpng = fullfile(dirFig,filenamepng);
    print(gcf,filepathpng, '-dpng');
    print(gcf,filepathpng, '-depsc');


    %% Boxplots
    fs = 12;
    figure(figNo); figNo = figNo + 1;
    cl = lines;
    objl = {'o','s','d','>'};
    ms = 12:-2:2;
    aTableBest = table;

 
    %'C_{DEL,Blades}','C_{DEL,Tower}'
    yLabBl = 'C_{Bl}'; %'C_{DEL,Blades}';
    yLabTwr = 'C_{Tow}'; %'C_{DEL,Tower}';

    for idx = 1: length(vecRateLim)
        aRateLim = vecRateLim(idx);
        structTableRateLim.(sprintf('RateLim%02d', aRateLim)) = tableObj(tableObj.RateLim == aRateLim,:);
        aTable = structTableRateLim.(sprintf('RateLim%02d', aRateLim));

        ax1(1) = subplot(2,1,1);
        plot(aTable.bladeDev, (100 - aTable.bladeIPC)/100 ,objl{idx},'Color',cl(idx,:),...
            'MarkerSize',ms(1), 'LineWidth',1.5);
        hold on; grid on;
        ylabel(yLabBl); %Blade DEL IPC/CPC %')

        ax1(2) = subplot(2,1,2);
        plot(aTable.bladeDev, (100 - aTable.twrIPC)/100,objl{idx},'Color',cl(idx,:),...
            'MarkerSize',ms(1), 'LineWidth',1.5)
        hold on; grid on;
        ylabel(yLabTwr); %'Tower DEL(IPC) / DEL(CPC), %')

        [~, idxmin] = min((1- aTable.twrIPC/100) .^2 + (1- aTable.bladeIPC/100)  .^2);
        aTableBest = [aTableBest; aTable(idxmin,:)]; %#ok<AGROW>

    end
    legCell = regexp(sprintf('Rate lim. %d °/s,', vecRateLim),',','split');
    legCell0 = legCell(1:end-1);
    legend(legCell0)
    xlabel('Maximum Pitch Deviation, °')

    vecRateLim1 = unique(tableObjAll.RateLim);
    for idx = 1: length(vecRateLim1)
        aRateLim = vecRateLim1(idx);

        structTableRateLim0.(sprintf('RateLim%02d', aRateLim)) = tableObjAll(tableObjAll.RateLim == aRateLim,:);
        aTableAll = structTableRateLim0.(sprintf('RateLim%02d', aRateLim));

        plot(ax1(1), tableObjAll.bladeDev, (100 - tableObjAll.bladeIPC)/100,objl{idx+2},'Color',0.5 *ones(1,3),...
            'MarkerSize',ms(1), 'LineWidth',1.5)


        plot(ax1(2), aTableAll.bladeDev, (100 - aTableAll.twrIPC)/100,objl{idx+2},'Color',0.5 *ones(1,3),...
            'MarkerSize',ms(1),'LineWidth',1.5)

    end

    legCella = regexp(sprintf('Rate lim. %d° (C_{Act} = 6),', vecRateLim1),',','split');

    legall = [legCell0, legCella(1:end-1)];
    legend(legall,"Location","South",'FontSize',fs);

    set(ax1(1), 'FontSize',fs);
    set(ax1(2), 'FontSize', fs);

    pos0 = get(0,'defaultFigurePosition');
    set(gcf,'Position', [pos0(1:3),pos0(4)*0.65])


    filenamepng = 'fig01_CmpIPCBlTwr';
    filepathpng = fullfile(dirFig,filenamepng);
    print(gcf,filepathpng, '-dpng');
    print(gcf,filepathpng, '-depsc');

    dirFig = dirFigDiss;



%% Plot MPC criteria
%return;
%load('matSubPlotBar.mat', 'K_matrix','strK','summaryTickLabel','printTitle');

load('matSubPlotBar_4to13_16.mat', 'K_matrix','strK','summaryTickLabel','printTitle');
selVec = [1,3:size(K_matrix,1)];
selVecMPC = selVec(2:end);

K_matrix = K_matrix(selVecMPC,:);
summaryTickLabel = summaryTickLabel(selVecMPC);
summaryTickLabel(1) = strrep(summaryTickLabel(2),'8','4');

orderCriteria1 = [4,5,1,2,3]; %Blades, Twower, Power, power std, Ac
orderCriteria2 = [4,5,1,2,3,6];
legCell = {'C_{mean(P)}','C_{std(P)}','C_{Act}','C_{DEL,b}','C_{DEL,t}','max(Blade_{Dev}), °'};
%legCell = {'C_{mean(P)}','C_{std(P)}','C_{Act}','C_{Bl}','C_{Tow}','max(\Delta_{Bl}), °'};

fcnPlotBar(K_matrix(:,orderCriteria1),strK,summaryTickLabel,printTitle,figNo,1.3,legCell(orderCriteria1))
set(gca, 'FontSize', fs);
set(findall(gcf,'-property','FontSize'),'FontSize',12)
% pos0 = get(0,'defaultFigurePosition');
% set(gcf,'Position', [pos0(1:2),pos0(3)*1.2,pos0(4)*0.8])

filenamepng = 'afig03_MPCbar';
filepathpng = fullfile(dirFig,filenamepng);
print(gcf,filepathpng, '-dpng');
print(gcf,filepathpng, '-depsc');

%% Plot only IPC criteria
K_matrixIPC0 = [1-aTableBest.bladeIPC/100,1-aTableBest.twrIPC/100, ...
    aTableBest.meanGenPwr,aTableBest.stdGenPwr, aTableBest.actEne];

summaryTickLabelIPC0 = strrep(summaryTickLabel,'MPC','IPC'); % [{'IPC, limit 4 °/s'}, ...
    %{'IPC, limit 8 °/s'}, {'IPC, limit 13 °/s'}];

fcnPlotBar(K_matrixIPC0,strK,summaryTickLabelIPC0,0,figNo+1000,1.3,legCell(orderCriteria1));
set(gca, 'FontSize', fs);
set(findall(gcf,'-property','FontSize'),'FontSize',12)

filenamepng = 'afig01_IPCbar';
filepathpng = fullfile(dirFig,filenamepng);
print(gcf,filepathpng, '-dpng');
print(gcf,filepathpng, '-depsc');


%% Plot IPC and MPC criteria
K_matrixIPCMPC = [K_matrixIPC0; K_matrix(:,orderCriteria1)];

summaryTickLabelIPCMPC = [summaryTickLabelIPC0, summaryTickLabel ]; % [{'IPC, limit 4 °/s'}, ...
    %{'IPC, limit 8 °/s'}, {'IPC, limit 13 °/s'}];

fcnPlotBar(K_matrixIPCMPC,strK,summaryTickLabelIPCMPC,0,figNo+1001,1.7,legCell(orderCriteria1));
set(gca, 'FontSize', 14);
set(findall(gcf,'-property','FontSize'),'FontSize',13)

filenamepng = 'afig01_IPCMPCbar';
filepathpng = fullfile(dirFig,filenamepng);
print(gcf,filepathpng, '-dpng');
print(gcf,filepathpng, '-depsc');


if onlyDissPics == 1
    return;
end


%% Plot MPC and IPC criteria
dirFig = dirFigDiss;

K_matrixIPC = [K_matrix(:,orderCriteria1);...
    K_matrixIPC0];

summaryTickLabelIPC = [summaryTickLabel,{'IPC, limit 4 °/s'}, ...
    {'IPC, limit 8 °/s'}, {'IPC, limit 13 °/s'}];

fcnPlotBar(K_matrixIPC,strK,summaryTickLabelIPC,0,figNo+1,2,legCell(orderCriteria1));
set(gca, 'FontSize', fs);

filenamepng = 'fig04_MPCIPCbar';
filepathpng = fullfile(dirFig,filenamepng);
print(gcf,filepathpng, '-dpng');
print(gcf,filepathpng, '-depsc');


%% Plot IPC criteria
K_matrixIPC1 = [K_matrixIPC0, aTableBest.bladeDev/10];

K_matrixIPC = [[K_matrix(1,orderCriteria1),0];...
    K_matrixIPC1];

summaryTickLabelIPC1 = [summaryTickLabel(1),{'IPC, limit 4°/s'}, ...
    {'IPC, limit 8°/s'}, {'IPC, limit 13 °/s'}];

fcnPlotBar(K_matrixIPC,strK,summaryTickLabelIPC1,0,figNo+3,1.3,legCell(orderCriteria2));

filenamepng = 'fig05_IPCbar';
filepathpng = fullfile(dirFig,filenamepng);
print(gcf,filepathpng, '-dpng');
print(gcf,filepathpng, '-depsc');

%% Plot IPC criteria
fcnPlotBar(K_matrixIPC1,strK,summaryTickLabelIPC1(2:end),0,figNo+4,1,legCell(orderCriteria2));
set(gca, 'FontSize', fs);
pos0 = get(0,'defaultFigurePosition');
set(gcf,'Position', [pos0(1:2),pos0(3)*1.2,pos0(4)*0.9])

filenamepng = 'fig06_IPCbar';
filepathpng = fullfile(dirFig,filenamepng);
print(gcf,filepathpng, '-dpng');
print(gcf,filepathpng, '-depsc');

%% Plot IPC and MPC criteria with blade deviation

K_matrixIPCMPC = [[K_matrix(1,orderCriteria1),0];...
    K_matrixIPC1; [K_matrix(2:end,orderCriteria1),zeros(size(K_matrix,1)-1,1)]];

summaryTickLabelIPCMPC =[{'PID CPC, limit 8°/s'}, ...
    {'PID IPC, limit 4°/s'}, ...
    {'PID IPC, limit 8°/s'}, {'PID IPC, limit 13°/s'}, {'MPC CPC, limit 4°/s'},...
    {'MPC CPC, limit 13°/s'}];
fcnPlotBar(K_matrixIPCMPC,strK,summaryTickLabelIPCMPC,0,figNo+5,2.6,legCell(orderCriteria2));
set(gca, 'FontSize', fs);

% Define x-coordinate for the vertical line and % Add vertical line
ylimVec = get(gca,'YLim');
xPos = [1.5; 4.5];
textContr = {'PID IPC','MPC CPC'};
for idx = 1: 2
    xline(xPos(idx), '--k', 'LineWidth', 1.5,'HandleVisibility', 'off');  % Dashed black vertical line

    % Add text near the vertical line

    text(xPos(idx) + 0.1, ylimVec(2)*0.9, textContr{idx}, 'FontSize', fs+4, 'Color', 'k');
end

pos1 = get(gcf,'Position');
pos0 = get(0,'defaultFigurePosition');
set(gcf,'Position', [pos1(1:2),pos1(3),pos0(4)*0.825])

filenamepng = 'fig07_IPCbar';
filepathpng = fullfile(dirFig,filenamepng);
print(gcf,filepathpng, '-dpng');
print(gcf,filepathpng, '-depsc');


%% Plot IPC and MPC criteria without blade deviation for dissertation
K_matrixIPCMPCdiss = K_matrixIPCMPC(:,1:end-1);
legCellOrd = legCell(orderCriteria2);
%legCellDiss = legCellOrd(1:end-1);

legCellDiss = {'C_{DEL,b}','C_{DEL,t}','C_{Pbar}','C_{\sigma(P)}','C_{Delta \beta}'};
summaryTickLabelIPCMPC =[{'PI CPC, limit 8°/s'}, ...
    {'PID IPC, limit 4°/s'}, ...
    {'PID IPC, limit 8°/s'}, {'PID IPC, limit 13°/s'}, {'MPC CPC, limit 4°/s'},...
    {'MPC CPC, limit 13°/s'}];

fcnPlotBar(K_matrixIPCMPCdiss,strK,summaryTickLabelIPCMPC,0,figNo+6,2.7,legCellDiss);
set(findall(gcf,'-property','FontSize'),'FontSize',12.5)

% Define x-coordinate for the vertical line and % Add vertical line
ylimVec = get(gca,'YLim');
xPos = [1.5; 4.5];
textContr = {'PID IPC','MPC CPC'};
for idx = 1: 2
    xline(xPos(idx), '--k', 'LineWidth', 1.5,'HandleVisibility', 'off');  % Dashed black vertical line

    % Add text near the vertical line

    text(xPos(idx) + 0.1, ylimVec(2)*0.9, textContr{idx}, 'FontSize', fs+4, 'Color', 'k');
end

pos1 = get(gcf,'Position');
pos0 = get(0,'defaultFigurePosition');
set(gcf,'Position', [pos1(1:2),pos1(3),pos0(4)*0.825])

filenamepng = 'fig08_IPCbar';
filepathpng = fullfile(dirFig,filenamepng);
print(gcf,filepathpng, '-dpng');
print(gcf,filepathpng, '-depsc');


%%


%
% [~, idxmin] =  min((1- aTableBest.twrIPC/100) .^2 + (1- aTableBest.bladeIPC/100)  .^2);
%
% rateLimStr = sprintf('MaxRate%d',aTableBest{idxmin,'RateLim'});
% KpStr = sprintf('Kp%s',aTableBest{idxmin,'KpName'}{:});
%
% DELratio_twrIPC = aTableBest{idxmin,'twrIPC'};
% DELratio_bladeIPC = aTableBest{idxmin,'bladeIPC'};
%
% idxRateKpFiles = contains(allMatFiles,rateLimStr) & contains(allMatFiles,KpStr);
%
% anOutTableIPCName = allMatFiles{idxIPCFiles & idxRateKpFiles};
% anOutTableCPCName = allMatFiles{~idxIPCFiles & idxRateKpFiles};
%
% load(anOutTableIPCName,'OutTable0');
% OutTableIPC = OutTable0(idxTime,:); % assign table as IPC table
% Mblade_all_IPC = OutTableIPC(:,idxBlM);
% Mblade_pitch_IPC = OutTableIPC(:,idxBlMZb);
% Mtower_IPC = OutTableIPC(:,idxT);
% Blade_pitch_IPC = OutTableIPC(:,idxPitch);
%
% load(anOutTableCPCName,'OutTable0');
% OutTableCPC = OutTable0(idxTime,:);
% Mblade_all_CPC = OutTableCPC(:,idxBlM);
% Mblade_pitch_CPC = OutTableCPC(:,idxBlMZb);
% Mtower_CPC = OutTableCPC(:,idxT);
% Blade_pitch_CPC = OutTableCPC(:,idxPitch);
%
% filenamepng = 'fig02_barPlots';
% filepathpng = fullfile(dirFig,filenamepng);
% print(gcf,filepathpng, '-dpng');
% print(gcf,filepathpng, '-depsc');

%% First visualizations

plotFASTresults(OutTableCPC,figNo,'CPC, NTM mean wind 16 m/s'); figNo = figNo + 1;
plotFASTresults(OutTableIPC,figNo,'IPC limit 13 °/s, NTM mean wind 16 m/s'); figNo = figNo + 1;

tend = 800;

OutTableIPC_t0 = OutTableIPC;
OutTableIPC_t0.Time = OutTableIPC_t0.Time - OutTableIPC_t0.Time(1);
OutTableIPC_t0 = OutTableIPC_t0(OutTableIPC_t0.Time <= tend,:);

OutTableCPC_t0 = OutTableCPC;
OutTableCPC_t0.Time = OutTableCPC_t0.Time - OutTableCPC_t0.Time(1);
OutTableCPC_t0 = OutTableIPC_t0(OutTableCPC_t0.Time <= tend,:);

plotFASTresults(OutTableIPC,figNo,'IPC limit 13°/s, NTM mean wind 16 m/s',OutTableCPC); figNo = figNo + 1;


%% Plot comparison in one plot
figure(figNo); figNo = figNo + 1;
plot(OutTableIPC.Time, OutTableCPC.TwrBsMyt,OutTableIPC.Time, OutTableIPC.TwrBsMyt,'--')
axis tight; grid on;
title(sprintf('DEL reduction tower %2.2f%%', DELratio_twrIPC));

figure(figNo); figNo = figNo + 1;
for idx = 1:3
    subplot(3,1,idx)
    rootStr = sprintf('Root%s%d',cStrAll,idx);
    plot(OutTableCPC.Time, OutTableCPC.(rootStr), OutTableIPC.Time, OutTableIPC.(rootStr));
    axis tight; grid on; ylabel([rootStr,' (Nm)']);

    if idx == 1
        title(sprintf('DEL reduction blades %2.2f%%', DELratio_bladeIPC));
    end
end
xlabel('time (s)')

figure(figNo); figNo = figNo + 1;
plot(OutTableIPC.Time, abs(OutTableIPC.(blStr1) - OutTableIPC.(blStr2)));
hold on;  plot(OutTableIPC.Time, abs(OutTableIPC.(blStr1) - OutTableIPC.(blStr3)));
hold on;  plot(OutTableIPC.Time, abs(OutTableIPC.(blStr2)- OutTableIPC.(blStr3)));
xlabel('time (s)'); ylabel('Deviation Angles (°)')
axis tight;


% %Get mean generator power
% meanGenPwr = nan( 1,lenComparison);
% for idx = 1: lenComparison
%     meanGenPwr(idx) = mean(cellOfTable{idx}.GenPwr(idxTime)/1000);
% end
%
% %Get std generator power
% stdGenPwr = nan( 1,lenComparison);
% for idx = 1: lenComparison
%     stdGenPwr(idx) = std(cellOfTable{idx}.GenPwr(idxTime));
% end


%figure; plot(OutTable.Time(1:end-1),diff(OutTable.(blStr1))/0.008);
% axis tight;
% hold on;  plot(OutTable.Time, abs(OutTable.(blStr1) - OutTable.(blStr3)));
% hold on;  plot(OutTable.Time, abs(OutTable.(blStr2)- OutTable.(blStr3)));
% xlabel('time (s)'); ylabel('Deviation Angles (°)')
% axis tight;


% figure; plot(OutTable0a.RootMyb2);hold on; plot(OutTable.RootMyb2)
% figure; plot(OutTable0a.RootMyb3);hold on; plot(OutTable.RootMyb3)

%% Unused code
% theta = OutTable.BldPitch1;

% Convert the azimuth angle to radians
% theta_rad = °2rad(theta);

% Define the rotation matrix
% R = [cos(theta_rad) sin(theta_rad) 0;
%     - sin(theta_rad)  cos(theta_rad) 0;
%     0               0              1];
% R = [cos(theta_rad) sin(theta_rad) ;
%    - sin(theta_rad)  cos(theta_rad)];

% M_b = Mt_b.(bStr){:,1:2};
% n = 2000; % height(OutTable);
%
% M_r = nan(2,n);
% t1 = Mt_c.(cStr);
% varnames = t1.Properties.VariableNames;
%
% if plotMbMcTrans == 1
%
%     for idx = 1: n
%         anM_b = M_b(idx,:)';
%
%         aTheta_rad = theta_rad(idx);
%
%         R = [cos(aTheta_rad) sin(aTheta_rad) ;
%             - sin(aTheta_rad)  cos(aTheta_rad)];
%
%         M_r(:,idx) = R*anM_b;%transformBladeToRotor(anM_b, theta);
%     end
%
%     figure;
%     subplot(2,1,1);
%     plot(M_r(1,1:n)); hold on; plot(Mt_c.(cStr){1:n,2},'--');
%     legend('R*M_b','M_c')
%     title(varnames{2})
%     subplot(2,1,2);
%     plot(M_r(2,1:n)); hold on; plot(Mt_c.(cStr){1:n,1},'--');
%     title(varnames{1})
%     legend('R*M_b','M_c')
% end

