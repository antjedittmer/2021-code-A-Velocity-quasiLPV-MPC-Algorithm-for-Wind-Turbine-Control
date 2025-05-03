function [wecs, M, Ce, K, Q, L, rho, tau, kappa, lambda, pitch, Cq, Ct, Q3, Cp] = initModel5MWNREL(plotOn, Rotor_Lamda, Rotor_Pitch, Rotor_cQ, Rotor_cT, figDir,titleOn,multRb)
% initModel5MWNREL initializes parameters of NREL 5 MW turbine.
% Aerodynamic force and thrust LUT can be plotted.
% All inputs are optional.
% - plotOn: Plots aerodynamic torque and force dependent on pitch and tip
%   speed ratio (Default: 0)
% - Rotor_Lamda: tip speed ratio vector (from mat file NREL5MW_CPdata)
% - Rotor_Pitch: blade pitch vector (from mat file NREL5MW_CPdata)
% - Rotor_cQ: aerodynamic torque look-up table
% - Rotor_cT: aerodynamic force look-up table
% - figDir: output figure directory (Default: 1)
%
% Outputs:
% - wecs: Structure with NREL 5 WM information
% - M: mass matrix Lagrange's equation
% - Ce: damping matrix Lagrange's equation
% - K: stiffness matrix Lagrange's equation
% - Q: term for input matrix based on Lagrange's equation
% - L: term for system matrix based on Lagrange's equation
% - rho: air density
% - tau: time constant pitch actuator
% - kappa: time constant torque actuator
% - Cq: aerodynamic torque
% - Ct: aerodynamic force

%% Set path for input data and output figure directory

debugOn = 0; % this is only used to check that cQ = cP/lambda
scaleLoopUp = 1; % this gives the possibility to scale cQ

% Set path to input data directory
workDir = fileparts(mfilename('fullpath'));
mainDir = fileparts(workDir);

% Switch for plots
if ~nargin
    plotOn = 1;
end

% Load data for cT, cQ LUT if not passed as input
if nargin < 5 || isempty(Rotor_Lamda) || isempty(Rotor_Pitch) || ....
        isempty(Rotor_cQ) || isempty(Rotor_Pitch)
    dataInDir = fullfile(mainDir,'dataIn');
    load(fullfile(dataInDir,'FASTToolCPnew.mat'),... % FASTToolCPnew','NREL5MW_CPdata.mat', 'NREL5MW_CPdata_BS.mat'
        'Rotor_Lamda','Rotor_Pitch','Rotor_cQ','Rotor_cT','Rotor_cP');
else
    dataInDir = fullfile(mainDir,'dataIn');
    load(fullfile(dataInDir,'NREL5MW_CPdata.mat'),...
        'Rotor_cP');
end


% Set path to figure directory
if nargin < 6 || isempty(figDir)
    figDir = fullfile(mainDir,'figDir');
    if ~isfolder(figDir)
        mkdir(figDir)
    end
end

% Put a plot title
if nargin < 7  || isempty(titleOn)
    titleOn = 0;
end

% Set a multiplication factor for the aerodynamic center on the blade
if nargin < 8  || isempty(multRb)
    multRb = 0.75;
end

%% Define parameters
% 5MW FAST WIND TURBINE (NREL/TP-500-38060) Definition of a 5-MW Reference
% Wind Tudeltaine for Offshore System Development
% Cut-In, Rated, Cut-Out Wind Speed 3 m/s, 11.4 m/s, 25 m/s
% Cut-In, Rated Rotor Speed 6.9 rpm, 12.1 rpm
% Tower equivalent mass, MT 438,000 kg
% Tower equivalent damping, CT 6421 wecs.N s/m
% Tower equivalent stiffness, KT 1,846,000 wecs.N/m

% Constants: Air density and actuator time constants
rho = 1.225; % Air density (kg/m^3)
tau = 0.1; % time constant pitch actuator
kappa = 0.01; % time constant torque actuator
% Bg = 0.9; % Tg= Bg(wg - wz). Unused because we use Tg as input
%x = [0.3540, 6.9975e+03]; % [0.608380, 7128.534769];
x = [1,1,1];
%x = [3.2950    3.9394    0.1974];
%x =   [4.1814    0.7651];
% x = [1.2745    0.6876];
%x =  [2.1336    0.8631   100];


% Turbine constants
% Turbine constants
wecs.N = 3;% Number of blades, Table 1-1
wecs.Ng = 97; % Gearbox ratio
wecs.mh = 56780; % kg;  Hub mas, Table 4-1
wecs.mbl = 17740; %kg;  Mass of each blade Table 2-2
wecs.mb = x(1)*wecs.mbl*0.25; %kg;  Modal mass of each blade
wecs.mn = 240000; % kg;  Nacelle mass Table 1-1
wecs.mtower = 347460; % kg; Tower mass Table 1-1
wecs.mr = wecs.mh + wecs.N*wecs.mbl; % kg; Rotor Mass: 110000 Table 1-1
wecs.mt = 0.25*wecs.mtower + wecs.mn + wecs.mr;
wecs.mtb = wecs.mt + wecs.N * wecs.mb; %tower modal mass + modal mass of blades

wecs.H =  90; % 87.6;

wecs.Jg = 534.116; %  kg*m^2; Inertia of the generator
wecs.Jr = x(2) * 3.8759e+07; % 3.8759e+07; %
wecs.JrL = 115926 + 3 * 11.776e6; % kg*m^2; Inertia of the rotor (Hub inertia + 3 blades)
wecs.Js = wecs.Jr + wecs.Ng^2*wecs.Jg;
f0 = 0.324;  % Hz, First natural tower fore-aft frequency
f0sw = 0.3120;% ; % First natural tower sidewards frequency
wecs.wnb = 0.6993 * 2*pi; % rad/s First natural blade frequency
wecs.wnt = f0 * 2*pi; % wecs.wnb; 0.3240
wecs.wntsw = f0sw * 2*pi; % wecs.wnb; 0.3240
wecs.zetat = 1/100; % damping ratio of tower (Table 6.2)
wecs.zetab = x(3) * 0.477465/100;% damping ratio of blade


wecs.Kt = wecs.wnt^2 * wecs.mt; % Stiffness of the tower s^2 + B/Ms + K/m
wecs.Bt = 2 * wecs.zetat * wecs.wnt *wecs.mt; % tower 2*6421;

wecs.Ktsw = wecs.wntsw^2 * wecs.mt; % Stiffness of the tower s^2 + B/Ms + K/m
wecs.Btsw = 2 * wecs.zetat * wecs.wntsw *wecs.mt; % tower 2*6421;

wecs.Kb = wecs.wnb^2 * wecs.mb; %Stiffness of each blade
wecs.Bb = 2 * wecs.zetab* wecs.wnb *wecs.mb; %Damping of the blade
wecs.Ks = 867637000; %Nm/rad Stiffness of the transmission
% 2*zeta*wn = B
wecs.Bs = 6215000;  %Nm/rad/sec %Damping of the transmission
wecs.Rr = 63; % m length rotor radius
wecs.rb = wecs.Rr*multRb; %*0.75; % m aerodynamic center on blade radius
wecs.etag = 0.944; %Drivetrain.Generator.Efficiency: 0.944;

%% Lagrange's Model matrices
% Force input w = [Ft_fa,Ft_sw,Tr,Tg]; Ft_fa = Ft, Ft_sw = 3/2*Tg
% States q: xdot_fa, zeta, xdot_sw, omega_r, omega_gr

M =[wecs.mtb wecs.N*wecs.mb*wecs.rb 0 0 0; %
    wecs.N*wecs.mb*wecs.rb  wecs.N*wecs.mb*wecs.rb^2 0 0 0;
    0 0 wecs.mt 0 0;
    0 0 0 wecs.JrL 0;
    0 0 0 0 wecs.Jg*wecs.Ng^2];

Ce = [wecs.Bt 0 0 0 0;
    0 wecs.N*wecs.Bb*wecs.rb^2 0 0 0 ; %
    0 0 wecs.Btsw 0 0
    0 0 0  wecs.Bs -wecs.Bs;
    0 0 0 -wecs.Bs wecs.Bs];

K = [wecs.Kt 0 0 0;
    0 wecs.N*wecs.Kb*wecs.rb^2 0 0; %
    0 0 wecs.Ktsw 0
    0 0 0 wecs.Ks;
    0 0 0 -wecs.Ks];

Q = [1 0 0 0;
    wecs.rb 0 0 0;
    0 1 0 -wecs.Ng;
    0 0 1 0;
    0 0 0  -wecs.Ng];

Q3 = [1 0 0;
    wecs.rb 0 0;
    0 2/(3*wecs.H) -2*wecs.Ng/(3*wecs.H);
    0 1 0;
    0 0 -wecs.Ng];

L = [eye(4), [0;0;0;-1]];

%% Ct/Cq for LUT for Model
idxPitch = Rotor_Pitch >= -2;
idxTSR = Rotor_Lamda < 20;

betaDeg = Rotor_Pitch(idxPitch);
pitch = betaDeg * pi/180;
lambda = Rotor_Lamda(idxTSR);

Cq = Rotor_cQ(idxPitch,idxTSR);
Ct = Rotor_cT(idxPitch,idxTSR);
Cp = Rotor_cP(idxPitch,idxTSR);

if debugOn == 0
    Cq1 = Cp ./ lambda';
end

if scaleLoopUp == 1

    % Ct: Find the minimum and maximum values
    lookupTable = Ct;
    minValue = min(Ct(:)); % Find the minimum and maximum values
    maxValue = max(Ct(:));
    newMaxValue = maxValue; % Reduce the highest point if desired

    % Ct: Apply linear scaling to adjust all values smoothly
    scaledTable = (lookupTable - minValue) / (maxValue - minValue); % Normalize to [0, 1]
    Ct1 = scaledTable * (newMaxValue - minValue) + minValue; % Scale back to new range
    Ct = Ct1;

    % Cq: Find the minimum and maximum values
    lookupTable = Cq;
    minValue = min(lookupTable(:));
    maxValue = max(lookupTable(:));
    newMaxValue = maxValue * 0.9475; % Reduce the highest point by 0.95%  * 0.9475

    % Cq: Apply linear scaling to adjust all values smoothly
    scaledTable = (lookupTable - minValue) / (maxValue - minValue); % Normalize to [0, 1]
    Cq1 = scaledTable * (newMaxValue - minValue) + minValue; % Scale back to new range
    Cq = Cq1;

end


%% Plot aerodynamic torque and force dependent on pitch and tip speed ratio

% Plot not generated by default
if plotOn
    [Xq,Yq] = meshgrid(lambda,betaDeg);
   
    % Plot aerodynamic force C_T
    zlabelStr = 'Thrust coefficient c_T [-]';
    figStr = 'Ct_NRRLFAST5MW'; figN = 1; % Plot aerodynamic torque C_T
    plotCoefficients(Xq,Yq,Ct,zlabelStr,figN,titleOn,figDir,figStr);

    zlabelStr = 'Power coefficient c_P [-]'; % Plot aerodynamic torque C_P
    figStr = 'Cp_NRRLFAST5MW'; figN = 100;
    plotCoefficients(Xq,Yq,Cp,zlabelStr,figN,titleOn,figDir,figStr);
    
    zlabelStr = 'Torque coefficient c_Q [-]'; % Plot aerodynamic torque C_Q
    figStr = 'Cq_NRRLFAST5MW'; figN = 2;
    plotCoefficients(Xq,Yq,Cq,zlabelStr,figN,titleOn,figDir,figStr);
   
    if debugOn == 1
        figStr = 'Cq1_NRRLFAST5MW'; figN = 200;
        plotCoefficients(Xq,Yq,Cq1,zlabelStr,figN,titleOn,figDir,figStr);
    end
end
end

function plotCoefficients(Xq,Yq,C,zlabelStr,figN,titleOn,figDir,figStr)

xlabelStr = 'Tip speed ratio \lambda [-]';
ylabelStr = 'Pitch angle \beta_0 [deg]';

idnan = isnan(C);
C(idnan) = 0;
idx0 = C <= eps;
Cnan = C;
Cnan(idx0) = nan;

fs = 12; % front size
sw = 3; % step width for look-up grid

titleStrData1Interp = ['Aerodynamique ',zlabelStr,'  from NREL FAST 5 MW'];
if titleOn == 1, faceAlphaVal = 1; else, faceAlphaVal = 0.9; end

figure(figN ); %For debugging
filtX1 = 1:1:size(Xq,1);
filtX2 = 1:1:size(Xq,2);
surf(Xq(filtX1,filtX2),Yq(filtX1,filtX2),Cnan(filtX1,filtX2),...
    'FaceColor','interp','EdgeColor','none','FaceAlpha', faceAlphaVal); hold on
axis tight;
xlabel(xlabelStr); ylabel(ylabelStr ); zlabel(zlabelStr);
set(gca, 'YDir', 'reverse');

if titleOn == 1
    title(titleStrData1Interp)
    view(60,30)
else
    hold on;
    filtX1 = 1:sw:size(Xq,1);
    filtX2 = 1:sw:size(Xq,2);
    surf(Xq(filtX1,filtX2),Yq(filtX1,filtX2),C(filtX1,filtX2),'FaceColor','none');
    set(gca,'Fontsize',fs)
    hold off
    % view(-160,60)
    colorbar
    set(gca,'Fontsize',fs)
end

print(gcf,fullfile(figDir,figStr), '-dpng');
print(gcf,fullfile(figDir,figStr), '-depsc');
end
