function damageBlades = calculateDamageBlades(My,mCarbonFibre,idxRootMy)

if nargin == 1
    mCarbonFibre = 10;
end

if nargin == 3 % table contains more than three rows -> index selected
    My = My{:,idxRootMy};
end


damageBlades = calculateDamage(My(:,1),mCarbonFibre) + ... 
    calculateDamage(My(:,2),mCarbonFibre) + ...
    calculateDamage(My(:,3),mCarbonFibre);

% theta = outTableNRELorig.Azimuth * pi/180;
% Mout= [cos(theta) cos((theta+2*pi/3)) cos((theta+4*pi/3));...
%       sin(theta) sin((theta+2*pi/3)) sin((theta+4*pi/3))] .* [My1,My2,My3];
% M_yaw  = 2/3*Mout(1);
% M_tilt = 2/3*Mout(2);

function damage = calculateDamage(sigt,m)
% S-N curve parameters:
% Inputs 
% sigt: signal
% m: slope of the curve blades: 10 (CFK), tower 4 (steel)

% tp = sig2ext(sigt);    % turning points
% rf = rainflow2(tp);    % rainflow
% CycleRate = rf(3,:);   % number of cycles
% siga = rf(1,:);        % cycle amplitudes
% 
% % calculation of the damage
% damage = sum((CycleRate).*((siga).^m));

% expected time to failure in seconds
% T=To/damage;

%figure, rfhist(rf,30,'ampl')
%figure, rfhist(rf,30,'mean')
%figure, rfmatrix(rf,30,30)

% Postprocessing Parameter
WoehlerExponent                 = m;            % [-]   for steel: 4; for cfk: 10
N_REF                           = 2e6/(20*8760);% [-]   fraction of 2e6 in 20 years for 1h
c                       = rainflow(sigt);
Count                   = c(:,1);
Range                   = c(:,2);
damage         = (sum(Range.^WoehlerExponent.*Count)/N_REF).^(1/WoehlerExponent);


% expected time to failure in days
%disp(['Calculated fatigue life in days: ' num2str(T/3600/24)])






