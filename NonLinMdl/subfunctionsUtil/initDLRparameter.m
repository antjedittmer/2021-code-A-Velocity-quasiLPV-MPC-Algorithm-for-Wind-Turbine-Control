function initDLRparameter(onIPC,switchIPC)

if nargin < 1
    onIPC = 0;
end

if nargin < 2
    switchIPC = 3;
end

Cut_In_Speed          = 70.16224; %#ok<*NASGU> % in rad/s, ursprünglich VS_CtInSp. Ab dieser Windgeschwindigkeit wird Energie gewonnen -> Transitional generator speed (HSS side)
Speed_Start_Region1_5 = 91.21091; % in rad/s, ursprünglich VS_Rgn2Sp. Ab hier Regeln auf optimaler Leistungskennlinie (Region 2) Transitional generator speed (HSS side)
Speed_Region2_5       = 110.6186; % in rad/s, ursprünglich VS_SySp. Differenz aus dieser und istgeschwindigkeit ist eingangsgröße für slope 25 ->Synchronous speed of region 2 1/2 induction
Speed_Start_Region2   = 119.0138; % in rad/s, ursprünglich VS_TrGnSp Transitional generator speed (HSS side)
Speed_Start_Region3   = 121.6805; % in rad/s(entspricht 1161,9 rpm), ursprünglich VS_RtGnSp. Generatornenngeschwindigkeit -> Rated generator speed (HSS side)
v_Gen_nenn            = 122.9096; % in rad/s(entspricht 1173   rpm), ursprünglich PC_RefSpd -> Desired (reference) HSS speed for pitch controller
slope15    = 921.83015253225492464369482727412; % Torque/speed slope of region 1 1/2 cut-in t
slope25    = 3.9350e+03;  % Torque/speed slope of region 2 1/2 inductio
VS_Rgn2K   = 2.332287;    % Generator torque constant in Region 2 (HSS

%%%%%%%%%%%%%%%%%%%%% Pitchregelung %%%%%%%%%%%%%%%%%%%%%
MaxPit = 1.570796;    % PC_MaxPit -> Zur Berechnung der oberen Pitchwinkelbegrenzung -> Maximum pitch setting in pitch controller
PC_KI  = 0.008068634; % PC_KI
MinPit = 0;           %PC_MinPit -> Untere Pitchwinkelbegrenzung -> Minimum pitch setting in pitch controller
Pitch_change_max = 0.1396263; % Begrenzung der Veränderung des Pitchwinkels (in absolute value)
PC_KK = 0.1099965;   % Faktor zur Berechnung des Gain Scheduling (Pitch angle were the the derivative of the ..)
PC_KI = 0.008068634; % Gain des Integralanteils im Pitch-Regler im Nennbereich
PC_KP = 0.01882681;  % Gain des Proportionalanteils im Pitch-Regler im Nennbereich
VS_Rgn3MP = 0.01745329;  % Minimum pitch angle at which the torque is

%%%%%%%%%%%%%%%%%%%%%% Generatormomentregelung %%%%%%%%%%%%%%%%%%%%%%
M_Gen_max = 47402.91;       % Obere  Begrenzung des Generatormoments in Nm (VS_MaxTq)
M_Gen_min = 0;              % Untere Begrenzung des Generatormoments in Nm
M_Gen_change_max = 15000.0; % Begrenzung der Veränderung des Generatormoments in Nm/s (VS_MaxRat)
P_nenn = 5296610;           % Nennleistung in W (VS_RtPwr)

% Tiefpass für Generatorgeschwindigkeit:
a1 = 1;
a0 = 1;
b1 = 1.1;
b0 = 1;
b1Torque = b1;
b0Torque = b0;

%%%%%%%%%%%%%%%%%%%%%% Sonstige %%%%%%%%%%%%%%%%%%%%%
i_Getriebe = 97;        % Getriebeübersetzung, dimensionslos
GK_Start   = 1;         % Startparameter Gain Scheduling
eta        = 0.944;     % Wirkungsgrad
offest_1p  = 25*pi/180; % Offset für Theta IPC 1P und 2P
IPC_K_1P   = 10^-6;     % IPC Gain
T_Turmfilter = 0.1304*pi; % Zeitkonstante Turmfilter
c1           = 2*pi/(4.95);
c2           = 0;
d1           = 1;
d3           = 4*pi/(4.95)*pi/(4.95);
e0=100;
e1=1;
e2=20;
e3=100;

% Regler 2
% IPC_KI_1p = 0.00096466;
% IPC_KP_1p = 0.00022946;
% d2        = 7.9345;
% kp_Turm   = 6.63e-05;
% IPC_K_1P  = 4.0006e-06;

%%%%%% Regler 1 %%%%%
IPC_K_1P  = 10^(-06);
k_Pitch   = 1;

IPC_An = onIPC;
if switchIPC == 1
    %25:
    IPC_KI_1p = 1.0495e-05;
    IPC_KP_1p = 0.00042667;
    d2        = 1.4983;
    kp_Turm   = 1e-05;
elseif switchIPC == 2
    %12
    IPC_KI_1p = 0.00028;
    IPC_KP_1p = 0.02786;
    d2        = 0.3120;
    kp_Turm   = 0.0086982;
    %
    %15:
else
    IPC_KI_1p = 4.45e-05;
    IPC_KP_1p = -0.00010109;
    d2        = 1.925;
    kp_Turm   = 1.9362e-05;
    
end

IPC_KI_1p_tilt = IPC_KI_1p;
IPC_KP_1p_tilt = IPC_KP_1p;

%% Assign this to workspace
% ToDo should this be moved to DD?
hEntry = who;
for idx = 1: length(hEntry)
    assignin('base',hEntry{idx},eval(hEntry{idx}))
end
