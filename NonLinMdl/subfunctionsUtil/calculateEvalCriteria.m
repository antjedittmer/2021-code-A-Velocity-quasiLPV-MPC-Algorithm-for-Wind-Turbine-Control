function [K_matrix, strK] = calculateEvalCriteria(cellOfTable,summaryTickLabel,printTitle)
% calculate evaluation criteria for controllers 

if nargin < 3
    printTitle = 1;
end

% epected variablenames: GenPwr,GenSpeed,TwrBsM,BlPitchC
outTable1 = cellOfTable{1};
lenComparison = length(cellOfTable);
% Get quantitative criteria for controller evaluation
idxTime = outTable1.Time > 60; % 10 : height(outTable1); %;  %
varNames = outTable1.Properties.VariableNames;

isFAST = nan(lenComparison,1);
for idx = 1: lenComparison
    varNames = cellOfTable{idx}.Properties.VariableNames;
isFAST(idx) = any(contains(varNames,'GenSpeed'));
end

if nargin <2
    summaryTickLabel = {'CPC', 'IPC_{fminsearch}','IPC_{SimAnn}'};
    %summaryTickLabel = {'NREL_{FAST}', 'NREL_{OC3}','DLR','TUDelft'};
end

%% Get evaluation criteria

%Get mean generator power
meanGenPwr = nan( 1,lenComparison);
for idx = 1: lenComparison
    meanGenPwr(idx) = mean(cellOfTable{idx}.GenPwr(idxTime)/1000);
end

%Get std generator power
stdGenPwr = nan( 1,lenComparison);
for idx = 1: lenComparison
    stdGenPwr(idx) = std(cellOfTable{idx}.GenPwr(idxTime));
end

%Get std generator speed power

stdGenSpeed = nan( 1,lenComparison);
% if isFAST
% for idx = 1: lenComparison
%     stdGenSpeed(idx) = std(cellOfTable{idx}.GenSpeed(idxTime));
% end
% end

%C_damageBlades calculated via Rainflow count w flap. blade RootMyb signals
% Antje ToDo Using flapwise moments only: Also edgewise(x) and pitch (z)?
% Antje ToDo: Transform from rotating into fixed reference system?
if isFAST
    idxRootMyb = contains(varNames,'RootMyb') | contains(varNames,'RootMFlp');
    if any(idxRootMyb)
        m = 10;
        damageBlades = nan( 1,lenComparison);
        for idx = 1: lenComparison
            damageBlades(idx) = calculateDamageBlades(cellOfTable{idx}(idxTime,:),m,idxRootMyb);
        end
    end
    
else
    damageBlades = nan( 1,lenComparison);
    for idx = 1: lenComparison
        damageBlades(idx) = std(cellOfTable{idx}.zeta_dotdot(idxTime,:));
    end
end

% C_damageTower calculated via Rainflow count with tower base moments
if isFAST
    idxTowerM = contains(varNames,'TwrBsM');
    if any(idxTowerM) % Antje ToDo What to use for offshore/floating wind turbines?
        m = 4;
        damageTower = nan( 1,lenComparison);
        for idx = 1: lenComparison
            damageTower(idx) = calculateDamageBlades(cellOfTable{idx}(idxTime,:),m,idxTowerM);
        end
    else
        damageTower = nan(size(damageBlades));
    end
else
    damageTower = nan( 1,lenComparison);
    for idx = 1: lenComparison
        damageTower(idx) = std(cellOfTable{idx}.NcIMUTAxs(idxTime,:));
    end
end

% C_actuator is the magnitude of the power consumed by the actuators
idxBlPitchC = contains(varNames,'BlPitch');
 %flapwise moment (i.e., the moment caused by flapwise forces) at the blade root
% AD ToDo: Should this be in blade/rotating frame or fixed? Definition P_{act,ref} = P_{act,0} 
meanActPwr = nan( 1,lenComparison);
if isFAST
    idxRootMyz = contains(varNames,'RootMz'); %'RootMzc1') | contains(varNames,'RootMzb1'); RootFzb1
    for idx = 1: lenComparison
        if any(idxRootMyz)
            meanActPwr(idx) = ...
                mean(mean(abs(diff(cellOfTable{idx}{idxTime,idxBlPitchC}).*cellOfTable{idx}{idxTime,idxRootMyz}(1:end-1,:))));
        else
            meanActPwr(idx) = ...
                mean(abs(diff(cellOfTable{idx}{idxTime,idxBlPitchC})));
        end
    end
else
    for idx = 1: lenComparison
        meanActPwr(idx) = ...
            mean(abs(diff(cellOfTable{idx}{idxTime,idxBlPitchC})));
    end
end

% Opt Criteria
refPw = abs(meanGenPwr(1)/5 - 1);
refStdGenSpeed = stdGenSpeed(1); % max(stdGenSpeed);
refStdGenPwr = stdGenPwr(1); % max(stdGenPwr);
refMeanActPwr = meanActPwr(1); % max(meanActPwr); %Definition P_{act,ref} = P_{act,0} ?
refDamageBlades = damageBlades(1); %max(damageBlades);
refDamageTower = damageTower(1); %max(damageTower);

K_GenPwr  = abs(meanGenPwr/refPw - 1);
K_stdGenSpeed = stdGenSpeed/refStdGenSpeed;
K_stdGenPwr = stdGenPwr/refStdGenPwr;
K_ActPwr = meanActPwr/refMeanActPwr;
K_damageBlades = damageBlades/refDamageBlades;
K_damageTower = damageTower/refDamageTower;


strMeanGenPwr = sprintf('%2.3f, ', meanGenPwr);
strK.powerabs = sprintf('%2.3f, ', K_GenPwr);
%strStdGenSpeed = sprintf('%2.3f, ', stdGenSpeed);

K_matrix = [K_GenPwr/K_GenPwr(1); K_stdGenPwr; K_ActPwr; K_damageBlades; K_damageTower]';

strK.power = sprintf('%2.3f, ', K_matrix(:,1));
%strK.StdGenSpeed = sprintf('%2.3f, ',K_stdGenSpeed);
strK.StdPowSpeed = sprintf('%2.3f, ',K_matrix(:,2)); %K_stdGenPwr);
strK.ActPwr = sprintf('%2.3f, ',K_matrix(:,3)); %K_ActPwr);
strK.DamageBlades = sprintf('%2.3f, ',K_matrix(:,4)); %K_damageBlades);
strK.DamageTower = sprintf('%2.3f, ',K_matrix(:,5)); %K_damageTower);

% figure;
% bar([K_GenPwr/max(K_GenPwr);K_stdGenSpeed]');
% axis tight; grid on;
% titleCell = {'Comparison PI Torque/CPC Controller',...
%     ['Mean GenPwr [MW]: ',strMeanGenPwr(1:end-2)],...
%     ['C_{P}: ', strK.powerabs(1:end-2)],...
%     ['Power tracking C_{P}', '{\color[rgb]{',num2str(cl(1,:)),'}',strK.power(1:end-2),'}'],...
%     ['C_{StdGenSpeed} = |STD_{ctrl}/STD_{max}|: ',  '{\color[rgb]{',num2str(cl(2,:)),'}',strK.StdGenSpeed(1:end-2),'}']};
% title(titleCell)
% xlabel('Controller');
% ylabel('Minimization criteria [-]');
% legend('C_{Power}','C_{StdGenSpd}','Location','SouthEast')
% ax = gca;
% ax.XTickLabel = summaryTickLabel;


subPlotBar(K_matrix,strK,summaryTickLabel,printTitle);

function subPlotBar(K_matrix,strK,summaryTickLabel,printTitle)
figure;
cl1 = colormap('lines');
cl = [0.5,0.5,0.5;
    cl1(4,:);
    [1,0,0];
    cl1(1,:);
    cl1(5,:)];

b = bar(K_matrix);
legend('C_{P}','C_{\sigma(P)}','C_{Actuator}','C_{D,Blades}','C_{D,Tower}','Location','NorthWest');

for k = 1:5
    b(k).FaceColor = cl(k,:);
end



axis tight; grid on;
titleCellNew = {['Power tracking C_{P}: ', '{\color[rgb]{',num2str(cl(1,:)),'}',strK.power(1:end-2),'}'],...
    ['Power variance C_{\sigma(P)}: ',  '{\color[rgb]{',num2str(cl(2,:)),'}',strK.StdPowSpeed(1:end-2),'}'], ...
    ['Actuator power C_{Actuator}: ',  '{\color[rgb]{',num2str(cl(3,:)),'}',strK.ActPwr(1:end-2),'}'], ...
    ['Blade damage C_{D,Blades}: ',  '{\color[rgb]{',num2str(cl(4,:)),'}',strK.DamageBlades(1:end-2),'}'],...
    ['Tower damage C_{D,Tower}: ',  '{\color[rgb]{',num2str(cl(5,:)),'}',strK.DamageTower(1:end-2),'}']};

if printTitle
title(titleCellNew)
else
%     for idx = 1: length(b)
%         xtips1 = b(idx).XEndPoints;
%         ytips1 = b(idx).YEndPoints;
%         labels1 = string(round(b(idx).YData,1));
%         text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%             'VerticalAlignment','bottom')
%     end
end

ax = gca;
ax.XLim = ax.XLim + [-0.1 0.1];
ax.YLim = ax.YLim + [0 0.2];
ax.XTickLabelRotation = 22.5;

ax.XTickLabel = summaryTickLabel;