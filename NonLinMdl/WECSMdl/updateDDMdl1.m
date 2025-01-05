% script to update data dictionary
% clear;
function updateDDMdl1(multRb)

if ~nargin
    multRb = 0.75;
end

% Update DD
nameCell = {'beta';'Ce';'Cq';'Ct';'K5'; 'lambda'; 'M';'Q'; 'rho';'wecs'};

loadEqPoints = 0; % load equilibrium/reference points for wind speeds from file
if  loadEqPoints == 1    
    load('NREL5MW_linearised_4to25.mat','Lin');
    Lin_points = Lin;
    nameCell{end+1} = 'Lin_points';
end

plotOn = 0;
Rotor_Lamda = ''; Rotor_Pitch = ''; Rotor_cQ = ''; Rotor_cT = '';
figDir = ''; titleOn = 0;

[wecs, M, Ce, K, Q, L, rho, tau, kappa, lambda, beta,Cq,Ct,Q3] = ....
    initModel5MWNREL(plotOn, Rotor_Lamda, Rotor_Pitch, Rotor_cQ, Rotor_cT, figDir,titleOn,multRb); %#ok<*ASGLU> 
K5 = [K,-K(:,4)];%#ok<*ASGLU> 
DDNameCell = {'DD_test.sldd'; 'DD_Mdl1.sldd'};


for idxDD = 1: length(DDNameCell)
    DDName =  DDNameCell{idxDD}; %'DD_test.sldd'; %'DD_Mdl1.sldd'; %
    Q = Q3;
    DDDataNames = nameCell;

    myDictionaryObj = Simulink.data.dictionary.open(DDName);
    dDataSectObj = getSection(myDictionaryObj,'Design Data');

    for idx = 1: length(DDDataNames)
        try
            tempObj = getEntry(dDataSectObj,DDDataNames{idx});

            eval(['setValue(tempObj, ',DDDataNames{idx},')']);
        catch
        end
    end

    saveChanges(myDictionaryObj)
end

% listEntry(myDictionaryObj)
%   Design Data   beta              DD_Mdl1.sldd   2021-04-02 13:18   ditt_aj          double   
%   Design Data   Ce                DD_Mdl1.sldd   2021-04-02 13:18   ditt_aj          double   
%   Design Data   Cq                DD_Mdl1.sldd   2021-04-02 13:18   ditt_aj          double   
%   Design Data   Ct                DD_Mdl1.sldd   2021-04-02 13:18   ditt_aj          double   
%   Design Data   K5                DD_Mdl1.sldd   2021-04-02 13:18   ditt_aj          double   
%   Design Data   lambda            DD_Mdl1.sldd   2021-04-02 13:18   ditt_aj          double   
%   Design Data   M                 DD_Mdl1.sldd   2021-04-02 13:18   ditt_aj          double   
%   Design Data   Q                 DD_Mdl1.sldd   2021-04-02 13:18   ditt_aj          double   
%   Design Data   rho               DD_Mdl1.sldd   2021-04-02 13:18   ditt_aj          double   
%   Design Data   wecs              DD_Mdl1.sldd   2021-04-02 13:18   ditt_aj          struct   

% tmp = load('NREL5MW_linearised.mat');
% load('OutDataStep.mat');
% idxT = OutTable.Time >= 100;
% OutTableT = OutTable(idxT,:);
% 
% OutTableT.Time =  OutTableT.Time - OutTableT.Time(1);
% time = OutTableT.Time;
% 
% tmp.Lin.V = 4:25;
% 
% for idx = 1 : length(tmp.Lin.V)
%     timeIdx = time <= idx*100-10 & time >= idx*100-30;
%     Lin_points1.V(idx) = mean(OutTableT.Wind1VelX(timeIdx));
% 
%     Lin_points1.RSpeed(idx) = mean(OutTableT.GenSpeed(timeIdx)/97/(60/(2*pi)));
%     Lin_points1.Pitch(idx) = mean(OutTableT.BlPitch1(timeIdx)/(180/(pi)));
%     Lin_points1.Torque(idx) = mean(OutTableT.GenTq(timeIdx));
% 
% end
% 
% Lin.RSpeed = Lin_points1.RSpeed;
% Lin.Pitch = Lin_points1.Pitch;
% Lin.Torque = Lin_points1.Torque;
% Lin.V = Lin_points1.V;
% 
% Lin_points = Lin;
