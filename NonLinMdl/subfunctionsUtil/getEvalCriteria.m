function [W,K,ref,weightVec] = getEvalCriteria(aTable,ref,weightVec)
% calculate evaluation criteria for controllers 
% inputs: * Table from simulation, has to contain variable names:
%           GenPwr,GenSpeed,TwrBsM,BlPitchC,RootMyb/Flp
%         * reference value for optimization criteria (opt.)
%         * weightVec weights optimization criteria
%         
% outputs: * W: weighted sum of optimization criteria
%          * K: structure with optimization criteria
%        

%% Prepare inputs

% Handle different input numbers
if nargin <2 || isempty(ref)
        ref = [];
end

if nargin <3 || isempty(weightVec)
        weightVec = ones(5,1); % for test: weightVec = [1.25,1,0.5,1,1.25];
end

% Ignore artifacts from initialization and get variable names 
idxTime = aTable.Time > 15; % 10 : height(outTable1); %;  %
varNames = aTable.Properties.VariableNames;


%% Calculate evaluation criteria
%Get mean generator power
meanGenPwr = mean(aTable.GenPwr(idxTime)/1000);

%Get std generator power
stdGenPwr = std(aTable.GenPwr(idxTime));

%Get std generator speed power (not used for now)
%stdGenSpeed = std(aTable.GenSpeed(idxTime));

%C_damageBlades calculated via Rainflow count w flap. blade RootMyb signals
% Antje ToDo Using flapping moments only: Also edgewise(x) and pitch (z)?
idxRootMyb = contains(varNames,'RootMyb') | contains(varNames,'RootMFlp');
m = 10;
damageBlades = calculateDamageBlades(aTable{idxTime,idxRootMyb},m);

%C_damageTower calculated via Rainflow count with tower base moments 
idxTowerM = contains(varNames,'TwrBsM');
mT = 4;
damageTower = calculateDamageBlades(aTable{idxTime,idxTowerM},mT);

% C_actuator is the magnitude of the power consumed by the actuators
idxBlPitchC = contains(varNames,'BlPitch');
%flapwise moment (i.e., the moment caused by flapwise forces) at the blade root
% AD ToDo: Should this be in blade/rotating frame or fixed? Definition P_{act,ref} = P_{act,0} 
meanActPwr = mean(mean(abs(diff(aTable{idxTime,idxBlPitchC})))); %.*aTable{idxTime,idxRootMyb}(1:end-1,:))));

%% Set reference
refOut.Pw = abs(meanGenPwr/5 - 1); % this only is ok for above rated test signals
% refOut.StdGenSpeed = stdGenSpeed;% max(stdGenSpeed);
refOut.StdGenPwr = stdGenPwr;% max(stdGenPwr);
refOut.MeanActPwr = meanActPwr;% max(meanActPwr); %Definition P_{act,ref} = P_{act,0} ?
refOut.DamageBlades = damageBlades;%max(damageBlades);
refOut.DamageTower = damageTower;%max(damageTower);

%% Opt Criteria
if isempty(ref)
    ref = refOut;
end

K.GenPwr  = weightVec(1)*abs(meanGenPwr/5 - 1)/ref.Pw;
%K.stdGenSpeed = stdGenSpeed/ref.StdGenSpeed;
K.stdGenPwr = weightVec(2)*stdGenPwr/ref.StdGenPwr;
K.ActPwr = weightVec(3)*meanActPwr/ref.MeanActPwr;
K.damageBlades = weightVec(4)*damageBlades/ref.DamageBlades;
K.damageTower = weightVec(5)*damageTower/ref.DamageTower;

%% Calculate the weighted sum
W = sum(cell2mat(struct2cell(K)))/sum(weightVec);

 
