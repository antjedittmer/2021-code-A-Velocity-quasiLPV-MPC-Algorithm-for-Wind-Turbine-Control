function [sysOut,gapCell,sysCmp] = compareLinearModels(speedVec,figFolder,useActuatorStates,figNoAdd,createBodePlots,plotVisible,noOut,xWeights,allPlots)
% compareLinearModels compares two sets of linearized turbine models with
% linearized FASTTool models at different wind speeds.
% All inputs are optional.
%
% Mdl1 is a 5th order model (Tu Delft) with state vector :
% x = [omega_r yt xt yt_dot xt_dot]'
% omega_r is the rotor speed and yt and xt are the axial tower fore-aft and
% sidewards displacement at nacelle height.
% The input and output vectors are:
% u = [Tg_ref beta_ref V ]'
% y = [yt_dot delta_dot omega_r omega_gr]'
%
% Mdl2 is a 9th order model with state vector:
% x = [yt delta xt theta_s yt_dot delta_dot xt_dot omega_r omega_gr]'
% delta the angular blade out-of-plane displacement and theta_s is the slip
% angle. omega_gr is the generator speed in the low speed shaft reference
% frame. Input u and output y are defined as above
%
% Inputs:
% - speedVec: index vector for 22 windspeeds, 4 to 25 m/s (Default: [1, 22])
% - figFolder: figure folder (Default: figDir in mainDir)
% - useActuatorStates: include (1)/exclude(0) actuator states (Default: 0)
% - figNoAdd: Add integer to figure number (100 + index) (Default: 0)
% - createBodePlots: Bode plots for analysis (Default: 1)
% - plotVisible: Figure visible (Default: 'on')
% - noOut: Rotor speed and Tower Foreaft Acc (2), additional Side-Side (3, default)
% - xWeights: Weights for optimization: Modal blade mass, generator
%   inertia, and damping ratio of blade (Default: [1,1,1])
% - allPlots: Generate all plots, including norm plots (Default: 0)
% - plotOnCpCt: Plot the Cp and Ct
%
% Outputs:
% - sysOut: Mdl2 with actuator states and outputs beta, Tg (last wind speed)
% - gapCell: gap metrics (allPlots = 1)
% - sysCmp: cell (index of the wind speed) with struct V, Mdl1, Mdl2, FAST:
%   linear models with the inputs [Tg (kNm), beta (rad), V (m/s)] and the
%   outputs [omega_r (rad/s), TwrAcc_fa (m/s^2), TwrAcc_sw (m/s^2)]

% Bindu Sharan, Antje Dittmer, ICS TUHH
% TU Delft model (ll. 175- 222): Atindriyo K. Pamososuryo and Jan-Willem
% van Wingerden (with small additions: Antje Dittmer)

%% Handle optional inputs
% The default inputs are provided here.

% Get relative paths to main directory (for data input dir)
workDir = fileparts(mfilename('fullpath'));
mainDir = fileparts(workDir);

% speedVec: index vector for 22 windspeeds, 4 to 25 m/s (Default: [1, 22])
if nargin == 0 || isempty(speedVec)
    speedVec = 8; %[1,8,9,22]; % variations of wind speed
end

if nargin < 2 || isempty(figFolder),figFolder = fullfile(mainDir,'figDir');
end

% Create output folder
if ~isfolder(figFolder),mkdir(figFolder); end

if nargin < 3 || isempty(useActuatorStates), useActuatorStates = 0; end

if nargin < 4 || isempty(useActuatorStates), figNoAdd = 0; end

if nargin < 5 || isempty(createBodePlots), createBodePlots = 1; close all;
end

if nargin < 6 || isempty(plotVisible),  plotVisible = 'on'; end

if nargin < 7 || isempty(noOut), noOut = 3; end

if nargin < 8  || isempty(xWeights), xWeights = [1,1,1]; end

if nargin < 9 || isempty(allPlots), allPlots = 0; end

if nargin < 10, plotOnCpCt = 0; end



%% Load and initialize name of models and load into workspace
addpath(fullfile(mainDir,'dataIn'));
load('NREL5MW_CPdata','Rotor_Lamda', 'Rotor_Pitch', 'Rotor_cQ', 'Rotor_cT','Rotor_cP');

% Initialize name of models and load into workspace
idxModel = 1; %in case several linearized models are provided
modelNames = {'NREL5MW_linearised_4to25'}; %{'NREL5MW_linearised_AD'}; %
modelNames = modelNames(idxModel);

modelsL = {'FAST 30 states'};
modelsL = modelsL(idxModel);

lenModelNames = length(modelNames);
sysmCell = cell(lenModelNames,1);

load(modelNames{1},'sysm','Lin'); % 1st model: state space and input OP
sysmCell{1} = sysm;
for idxModels = 2 : lenModelNames
    load(modelNames{idxModels},'sysm');
    sysmCell{idxModels} = sysm;
end

%% Get parameters
% 5MW FAST WIND TURBINE (NREL/TP-500-38060) Definition of a 5-MW Reference
% Wind Turbine for Offshore System Development
% Cut-In, Rated, Cut-Out Wind Speed 3 m/s, 11.4 m/s, 25 m/s
% Cut-In, Rated Rotor Speed 6.9 rpm, 12.1 rpm

% Wind energy conversion system (WECS) parameters
[wecs, M, Ce, K, Q, L, rho, tau, kappa, lambda, pitch, Cq, Ct ] = ...
    initModel5MWNREL(plotOnCpCt, Rotor_Lamda, Rotor_Pitch, Rotor_cQ, Rotor_cT, Rotor_cP, figFolder,0.75,xWeights); %  
wecs.Js = wecs.Jr + wecs.Ng^2 * wecs.Jg;

% Blade model of Mdl2 (getBladeModelVersion): 'original' (thesis, default)
% or variant B ('B1', 'B2': parameters in wecs.blade*, initModel5MWNREL);
% with bladeModelCompare the original Mdl2 is shown additionally (Mdl2
% variant B: yellow, dash-dotted)
[bladeVersion, bladeCompare] = getBladeModelVersion;
if strcmp(bladeVersion, 'original')
    mdl2Versions = {'original'};
else
    if bladeCompare
        mdl2Versions = {'original', bladeVersion};
    else
        mdl2Versions = {bladeVersion};
    end
end

%% Compute aerodynamic force and torque gradients
ksw = 1/(2/3* wecs.H);
dLambda = mean(diff(lambda));
dBeta = mean(diff(pitch));
[Cqdlambda,Cqdbeta] = gradient(Cq,dLambda,dBeta);
[Ctdlambda,Ctdbeta] = gradient(Ct,dLambda,dBeta);

% Constant matrices
% Set airspeed independent system and input matrix
Qxt = zeros(size(Q,1), size(Q,2)-1); % Qxt: FT, Tr, Tg
Qxt(3,2) = Q(3,2)*ksw; % Q3: FT, Ftx, Tr, Tg
Qxt(3,3) = Q(3,4)*ksw;

Q3 = Qxt; % Q3: FT, Tr, Tg
Q3(1:2,1) = Q(1:2,1);
Q3(4:5,:) = Q(4:5,2:end);
B1_3 = [zeros(4,3); M\Q3];

% Options for Bode plots
opts = bodeoptions('cstprefs');
opts.FreqUnits = 'Hz'; opts.MagUnits ='Abs'; opts.MagScale ='log'; opts.XLim ={[1e-2 2]};
if isappdata(groot, 'ShowThesisTitle') && ~getappdata(groot, 'ShowThesisTitle') % title switch of runGenerateAllThesisFigures
    opts.Title.String = ''; % no default title 'Bode Diagram'
end

% Constant matrices
% Set airspeed independent system and input matrix
A1 = [zeros(4), L;...
    -M\K, -M\Ce];

% Initialize natural frequencies, damping, pole vectros
lenLinV = length(Lin.V);
sysCell = cell(lenLinV,1);
u_bar_cell = cell(lenLinV,1);
x_bar_cell = cell(lenLinV,1); % X(s) = (sI - A)^1 *B *U(s), xdot = 0 = Ax +Bu
y_bar_cell = cell(lenLinV,1);

% Get and check FAST input and output names
idxFASTOutputAll = [38 12 13];
idxFASTOutput = idxFASTOutputAll(1:noOut);
idxFASTInput = [8 9 1];

% Set plot outputs
sysOutputname =[{'RotSpd \omega_r (rad/s)'};{'TwrAcc_{fa} (m/s^2)'}; {'TwrAcc_{sw} (m/s^2)'}];
sysOutputname = sysOutputname(1:noOut);
sysInputname =[{'GenTq T_g (kNm)'}; {'BlPitch \beta_0 (rad)'};{'Wind V_{\infty} (m/s)'}];

% Check FAst names
sysOutNamesOrig = sysm{1}.OutputName(idxFASTOutput);
disp(['TU Delft sysOutputnames: ', sprintf('%s, ', sysOutputname{:})])
disp('FAST sysOutputnames:')
fprintf('%s\n', sysOutNamesOrig{:});
fprintf('\n')
disp(['TU Delft sysInputnames: ', sprintf('%s, ', sysInputname{:})])
sysInNamesOrig = sysm{1}.InputName(idxFASTInput);
disp('FAST sysInputnames: ');
fprintf('%s\n', sysInNamesOrig{:})

% Set output gains for FAST accordingly
multFASTOutputAll = [pi/30,1,1]; % ED RotSpeed, (rpm)
multFASTOutput = multFASTOutputAll(1:noOut);
multFASTInput = [1000,1,1];
multFASTInput0 = multFASTInput;

% Initialize gap metrics for quantitative evaluation
noVelMdl = length(speedVec);
% noOutMdl = length(sysOutputname);
% noInMdl = length(sysOutputname);

gapCell.gap9DoF = cell(noVelMdl,1);
gapCell.nugap9DoF = cell(noVelMdl,1);
gapCell.gap5DoF = cell(noVelMdl,1);
gapCell.nugap5DoF = cell(noVelMdl,1);
gapCell.norm9DoF = cell(noVelMdl,1);
gapCell.norm5DoF = cell(noVelMdl,1);
sysCmp = cell(lenLinV,1);


%% Loop over all air speeds in LinV
for index =  speedVec

    %% Variables dependent on wind speed V
    Vbar = Lin.V(index); % hub height wind speed
    % Rotor-effective wind speed for the aerodynamics (rotor tilt, precone,
    % rotor averaging), as in the Simulink models: Ve = kWindEff*V
    kWind = 1;
    if isfield(wecs, 'kWindEff'), kWind = wecs.kWindEff; end
    Vebar = kWind*Vbar;
    Lambda_bar = Lin.RSpeed(index)*wecs.Rr/Vebar; % (-)
    omegabar = Lin.RSpeed(index); % 1/rad
    beta_bar = Lin.Pitch(index); % rad

    %% Grid points
    Cqb = interp2(lambda,pitch,Cq,Lambda_bar,beta_bar,'cubic');
    Ctb = interp2(lambda,pitch,Ct,Lambda_bar,beta_bar,'cubic');
    dCqdlambdab = interp2(lambda,pitch,Cqdlambda,Lambda_bar,beta_bar,'cubic');
    dCtdlambdab = interp2(lambda,pitch,Ctdlambda,Lambda_bar,beta_bar,'cubic');
    dCqdbetab = interp2(lambda,pitch,Cqdbeta,Lambda_bar,beta_bar,'cubic');
    dCtdbetab = interp2(lambda,pitch,Ctdbeta,Lambda_bar,beta_bar,'cubic');

    %% Tip speed ratio lambda = (Rr * omega)/Ve =  (Rr * omega)/(V - ydotfa)
    dlambdadomega = wecs.Rr/Vebar;
    dlambdadV = - omegabar*wecs.Rr/Vebar^2; % derivative w.r.t. Ve = kWind*V - ydot
    dlambdadydotfa =  - dlambdadV; %s/m   ydotfa m/s
    % dlambdaddelta =  - dlambdadV; %s/m   ydotfa m/s

    %% Ft = 0.5*rho*pi*Rr^2*Ct*v^2
    kCT = 0.5*rho*pi*wecs.Rr^2; % constant for CT: 0.5*rho*A
    dFtdomega = kCT *dCtdlambdab*dlambdadomega*Vebar^2;
    dFtdV = kWind * kCT * Vebar * (dCtdlambdab*dlambdadV*Vebar + 2*Ctb); % dVe/dV = kWind;
    dFtdydotfa = kCT * Vebar *(dCtdlambdab*dlambdadydotfa*Vebar - 2*Ctb); % dVe/dydot = -1;
    dFtdbeta = kCT*dCtdbetab*Vebar^2;

    %% Tr = 0.5*rho*pi*Rr^3*Cq*v^2
    kCQ = kCT*wecs.Rr;
    dTrdomega = kCQ *dCqdlambdab*dlambdadomega*Vebar^2;
    dTrdV = kWind * kCQ *Vebar *(dCqdlambdab*dlambdadV*Vebar + 2*Cqb); % dVe/dV = kWind
    dTrdydotfa = kCQ * Vebar*(dCqdlambdab*dlambdadydotfa*Vebar - 2*Cqb);
    dTrdbeta = kCQ *dCqdbetab*Vebar^2;

    %% Linear Model Mdl1: System, input matrices and state space model
    % States omega ydotfa xdotsw xfa xss

    % dotomega = (Tr - Ng * Tg)/J
    A11 = 1/wecs.Js * dTrdomega;
    A12 = 1/wecs.Js * dTrdydotfa;

    % ydotdotfa = 1/Mt( -Bt * ydotfa - Kt* yfa +Ft)
    % Mdl1 has no blade state: the tower-top mass wecs.mt already contains the
    % full rotor mass (as in Mdl_TuDelftOL). The blade modal mass in
    % wecs.mtb is only used in Mdl2, where the blade is a separate state.
    A21 = 1/wecs.mt * dFtdomega;
    A22 = 1/wecs.mt * dFtdydotfa;

    % xdotdotsw = 1/Mt( -Btsw * xdotsw - Ktsw* xsw + F_Ty) with
    % F_Ty = -3/(2*H)*M_x and the torsional nacelle moment
    % M_x = (1 - cMx)*Tr + cMx*Ng*Tg, cMx = (Jr + Ng*Jg)/(Jr + Ng^2*Jg)
    cMx = wecs.cMx;
    A31 = -(1 - cMx)*ksw/wecs.mt * dTrdomega;
    A32 = -(1 - cMx)*ksw/wecs.mt * dTrdydotfa;

    A = [A11  A12                  0                 0                0;
        A21 (A22 - wecs.Bt/wecs.mt) 0                -wecs.Kt/wecs.mt  0;
        A31  A32                  -wecs.Btsw/wecs.mt  0               -wecs.Ktsw/wecs.mt;
        [0 1 0 0 0; 0 0 1 0 0]];

    % Inputs: Torque [Nm], pitch angle beta [rad], wind speed V [m/s]
    B11 = -1/wecs.Js * wecs.Ng;
    B12 = 1/wecs.Js * dTrdbeta;
    B13 = 1/wecs.Js * dTrdV;

    B22 = 1/wecs.mt * dFtdbeta;
    B23 = 1/wecs.mt * dFtdV;

    B31 = -cMx*ksw/wecs.mt * wecs.Ng;
    B32 = -(1 - cMx)*ksw/wecs.mt * dTrdbeta;
    B33 = -(1 - cMx)*ksw/wecs.mt * dTrdV;

    B = [B11 B12 B13;
        0 B22 B23;
        B31 B32 B33;
        zeros(2,3);]*diag(multFASTInput0);

    % Outputs: omega [rad/s], fore-aft dx_dot_f, fore-aft dx_dot_fa [m/s^2]
    nO = min(noOut,3);
    C(2:nO,:) = A(2:nO,:);
    C(1,1) = 1;
    D(2:nO,:)= B(2:nO,:);

    sys5DoF = ss(A,B,C,D);

    % For debugging
    % testMatrix = ones(3); testMatrix(3,1) = 100;
    % tfSYS =  testMatrix .* tf(sys5DoF); sys5DoF = ss(tfSYS);

    sys5DoF.OutputName = sysOutputname;
    sys5DoF.InputName = sysInputname;

    %% Linear Model Mdl2: System and input matrices
    % Set up system and input matrices dependent on airspeed, pitch and rotor speed  x_FA

    sysLagrangeAll = cell(size(mdl2Versions));
    for idxBl = 1:numel(mdl2Versions)
    % Blade model matrices (original: M, K, Ce, B1_3, wecs.rb; variant B:
    % at the operating point, see bladeMatricesMdl2)
    Ftbar = kCT*Ctb*Vebar^2; % thrust at the operating point (N)
    [M2, K2, Ce2, B1_3m, rb2, dG2dsched] = bladeMatricesMdl2(mdl2Versions{idxBl}, ...
        wecs, M, K, Ce, Q3, B1_3, Lambda_bar, beta_bar, Ftbar);

    C5aero =  [B1_3m(5,1) * [-dFtdV, -dFtdV*rb2, 0, dFtdomega 0]; ... % yFAdot
        B1_3m(6,1) * [-dFtdV, -dFtdV*rb2, 0, dFtdomega 0]; ... % deltaDot
        B1_3m(7,2) * [-dTrdV, -dTrdV*rb2, 0, dTrdomega,0];... % xSWdot
        B1_3m(8,2) * [-dTrdV, -dTrdV*rb2, 0, dTrdomega,0];... % omega_r
        zeros(1,5)];  % omega_gr

    % B1_3 has the order: F_T, T_r, T_g
    QTildaTg =  [0;0; B1_3m(7,3); 0; B1_3m(9,3)];
    QTildaBeta = [B1_3m(5,1) *dFtdbeta; B1_3m(6,1)*dFtdbeta; B1_3m(7,2) * dTrdbeta;...
        B1_3m(8,2)*dTrdbeta;0];
    QTildaV = [B1_3m(5,1) * dFtdV; B1_3m(6,1) * dFtdV; B1_3m(7,2) * dTrdV;...
        B1_3m(8,2)*dTrdV;0];

    A9DoF = zeros(size(A1) + 2); % 9 states + 2 ctrl input
    A9DoF(1:4,:) = [A1(1:4,:) ,zeros(4,2)]; % states 1 - 4: y_FA,delta,x_sw, Theta_s
    A9DoF(5:9,1:4) = -M2\K2;
    A9DoF(5:9,5:9) = -M2\Ce2 + C5aero; % states 5 - 9: y_FA,delta,x_sw, w_r, w_gr
    A9DoF(5:9,10:11) = [QTildaTg, QTildaBeta];
    %A9DoF(10:11,:) = zeros(2,11); % states 10 - 11: T_g, beta
    A9DoF(10,10) =  -1/kappa;
    A9DoF(11,11) =  -1/tau;

    B9DoF = zeros(size(A1,1) + 2,3); % 3 inputs -> 9 states + 2 ctrl input
    B9DoF(5:9,3) = QTildaV;

    % Variant B: generalized blade force r_b(s)*F_T - N*Kb(s)*r_ref^2*theta
    % depends on the scheduling variables s = (lambda0, beta), lambda0 =
    % Rr*omega_r/(Ve - ydot_fa) (zero for the original blade model)
    sBl = M2\[0; 1; 0; 0; 0];
    A9DoF(5:9,8) = A9DoF(5:9,8) + sBl*dG2dsched(1)*dlambdadomega; % omega_r
    A9DoF(5:9,5) = A9DoF(5:9,5) + sBl*dG2dsched(1)*dlambdadydotfa; % ydot_fa
    A9DoF(5:9,11) = A9DoF(5:9,11) + sBl*dG2dsched(2); % beta
    B9DoF(5:9,3) = B9DoF(5:9,3) + sBl*dG2dsched(1)*dlambdadV*kWind; % V
    B9DoF(10,1) = 1/kappa;
    B9DoF(11,2) = 1/tau;
    A7DoF = A9DoF(1:9,1:9);

    B7DoF = B9DoF(1:9,:);
    B7DoF(:,1) = A9DoF(1:9,8+2); % Torque Tg as input
    B7DoF(:,2) = A9DoF(1:9,9+2); % Pitch beta as input

    B7DoF(:,1) = A9DoF(1:9,8+2); % Torque Tg as input

    if noOut == 2
        idxOut = 5; % ouput: ydotdot_fa
    else
        idxOut = [5,7];% ouput: ydotdot_fa, xdotdot_sw
    end
    C9DoF = [zeros(1,5+2),     1     zeros(1,3) % state 8: omega_r
        A9DoF(idxOut,:)]; % state 5: xdot_FA: output xdotdot_FA
    C7DoF = C9DoF(:,1:9);

    D9DoF = [zeros(1,3); B9DoF(idxOut,:)];
    D7DoF = [zeros(1,3); B7DoF(idxOut,:)];

    %% Create state-space model for Mdl2
    if useActuatorStates == 1
        sysLagrange = ss(A9DoF,B9DoF*diag(multFASTInput0),C9DoF,D9DoF*diag(multFASTInput0)); % to give input Tg in kN/m
        Mdl2Str = ' Rot,Twr,Bld+Act ';
    else
        sysLagrange = ss(A7DoF,B7DoF*diag(multFASTInput0),C7DoF,D7DoF*diag(multFASTInput0)); % to give input Tg in kN/m
        Mdl2Str = ' Rot,Twr,Gen+Bld ';
    end

    sysLagrange.Outputname = sysOutputname;
    sysLagrange.Inputname = sysInputname;
    sysLagrangeAll{idxBl} = sysLagrange;
    end % blade model versions: the last one is Mdl2 (outputs, metrics)

    Cout = [C9DoF;
        [zeros(1,7+2) 1 zeros(1,1)];  % beta
        [zeros(1,8+2) 1]];% torque
    Dout = [D9DoF; zeros(2,3)];
    sysOut = ss(A9DoF,B9DoF*diag(multFASTInput0),Cout,Dout); % to give input Tg in kN/m

    sysOut.Outputname = [sysOutputname;{'beta'};{'Tg'}];
    sysOut.Inputname = sysInputname;

    sysCell{index} = sysOut;
    % Linear models for further comparisons (third output): Mdl1, Mdl2 and FAST
    modelFAST = diag(multFASTOutput)*sysmCell{1}{index}(idxFASTOutput,idxFASTInput)*diag(multFASTInput);
    modelFAST.OutputName = sysOutputname;
    modelFAST.InputName = sysInputname;
    sysCmp{index} = struct('V', Vbar, 'Mdl1', sys5DoF, 'Mdl2', sysLagrange, 'FAST', modelFAST);
    Torque_g_bar = 0.5*rho*pi*wecs.Rr^3*Cqb*Vebar^2/wecs.Ng; % Nm, Tr = 0.5*rho*pi*Rr^3*Cq*V^2 (as kCQ)
    u_bar_cell{index} = [Torque_g_bar, beta_bar,Vbar]'; % beta_bar is already in rad
    x_bar_cell{index}  = - sysLagrange.A\sysLagrange.B *u_bar_cell{index}; % X(s) = (sI - A)1 *B *U(s) (xdot = 0 = Ax +Bu)
    y_bar_cell{index}  = sysLagrange.C * x_bar_cell{index} + sysLagrange.D * u_bar_cell{index};

    %% Make a Bode plot
    if createBodePlots == 1
        ff = figure(index*100+3 + figNoAdd);
        ff.Position = [520   267   560   531];
        ff.Color = 'white'; ff.Visible = plotVisible;
        %wVec = logspace(opts.XLim{1}(1),opts.XLim{1}(2),100);

        bodemag(sysLagrangeAll{1},opts); hold on; grid on; % original Mdl2 if bladeCompare
        freq9DoF = linspace(opts.XLim{:}(1),opts.XLim{:}(2),100);
        mag9DoF = bode(sysLagrange,freq9DoF);
        bodemag(sys5DoF,'k--',opts); hold on; grid on;
        mag5DoF = bode(sys5DoF,freq9DoF);

        for indexModels=1:lenModelNames
            sysm =  sysmCell{indexModels};
            model = diag(multFASTOutput)*sysm{index}(idxFASTOutput,idxFASTInput)* diag(multFASTInput);
            model.OutputName = sysOutputname;
            model.InputName = sysInputname;
            bodemag(model,opts); hold on
            magFAST = bode(model,freq9DoF);
            grid on
        end

        % Review of the blade model: Mdl2 variant B in addition (yellow, dash-dotted)
        if numel(sysLagrangeAll) > 1
            linesBefore = findall(ff, 'Type', 'line');
            bodemag(sysLagrangeAll{end}, '-.', opts);
            clB = lines(3);
            set(setdiff(findall(ff, 'Type', 'line'), linesBefore), 'Color', clB(3,:));
        end

        vecStates = [1:5,16:18];

        Atest = model.A(vecStates,vecStates); %#ok<NASGU>
        sstest = ss(model.A(vecStates,vecStates),model.B(vecStates,:),model.C(:,vecStates),model.D); %#ok<NASGU>

        cl = lines;
        legCell = [ {'\color{black}Mdl1: Rot+Twr '};
            {['\color[rgb]{',num2str(cl(1,:)),'}Mdl2:',Mdl2Str]};
            modelsL];

        strExtr = '';
        for ilegC = 1: length(modelsL)
            strExtr = [strExtr,...
                '\color[rgb]{',num2str(cl(ilegC+1,:)),'}',' ',modelsL{ilegC}]; %#ok<AGROW> This fast
        end

        % Gap and nugap metric results
        tempNormMatrix5 = nan(size(sysLagrange));
        tempNormMatrix9 = nan(size(sysLagrange));
        if allPlots == 1
            for idxOut = 1:size(sysLagrange,1)
                for idxIn = 1:size(sysLagrange,2)
                    [gapCell.gap9DoF{index}(idxOut,idxIn),gapCell.nugap9DoF{index}(idxOut,idxIn)] = gapmetric(sysLagrange(idxOut,idxIn),model(idxOut,idxIn));
                    [gapCell.gap5DoF{index}(idxOut,idxIn),gapCell.nugap5DoF{index}(idxOut,idxIn)] = gapmetric(sys5DoF(idxOut,idxIn),model(idxOut,idxIn));

                    %% For debugging
                    % [gapCell.gap59DoF{index}(idxOut,idxIn),gapCell.nugap59DoF{index}(idxOut,idxIn)] =
                    % ...
                    % gapmetric(sysLagrange(idxOut,idxIn),sys5DoF(idxOut,idxIn));

                    mag_9DOF = 20*log10(squeeze(mag9DoF(idxOut,idxIn,:)));    % Squeeze to get a 1D array
                    mag_5DOF = 20*log10(squeeze(mag5DoF(idxOut,idxIn,:)));
                    mag_FAST = 20*log10(squeeze(magFAST(idxOut,idxIn,:)));

                    % figure; semilogx(freq9DoF,mag_9DOF,freq9DoF,mag_5DOF,'k--',freq9DoF,mag_FAST)
                    tempNormMatrix9(idxOut,idxIn) = norm(mag_9DOF - mag_FAST)/norm(mag_FAST);
                    tempNormMatrix5(idxOut,idxIn) = norm(mag_5DOF - mag_FAST)/norm(mag_FAST);

                end
            end

            gapCell.norm9DoF{index} = tempNormMatrix9;
            gapCell.norm5DoF{index} = tempNormMatrix5;
        end


        % Title string
        titleStr = sprintf('Bode mag. plot for V = %d, pitch = %2.2f, and TSR = %2.4f',Vbar,beta_bar,Lambda_bar);
        titleCell = {titleStr, [legCell{1},legCell{2},strExtr]};
        if ~isappdata(groot, 'ShowThesisTitle') || getappdata(groot, 'ShowThesisTitle') % title switch of runGenerateAllThesisFigures
            title(['V_{\infty}',sprintf(' = %d: %s', Vbar, titleCell{2})]); %title(titleCell)
        end

        % Figure information
        % file names with the suffix _Blade<version> for the blade model variant B
        strBlade = '';
        if ~strcmp(bladeVersion, 'original'), strBlade = ['_Blade', bladeVersion]; end
        figStr = [strrep(titleStr,'Bode mag. plot for V =','Figure_Bode V ='),strBlade,'.png'];
        figFolderStr = fullfile(figFolder,figStr);

        figFolderStrEps = fullfile(figFolder,[sprintf('BodeV%d',Vbar),strBlade]);
        posaxes =   get(0,'defaultFigurePosition');
        set(gcf,'Position',[posaxes(1:3),posaxes(4)*1.1]);
        setFigureSize(gcf, [1 1.1]); % export size independent of the screen

        % Line width 1 and frequency ticks 0.01, 0.3, 0.7 (first blade
        % frequency), 1 and 2 Hz (labels only on the axes that show x tick
        % labels, i.e. the bottom row)
        set(findall(ff, 'Type', 'line'), 'LineWidth', 1);
        fTicks = [0.01 0.3 0.7 1 2];
        axBode = findall(ff, 'Type', 'axes', 'Visible', 'on');
        for idxAx = 1:numel(axBode)
            hasLabels = ~isempty(get(axBode(idxAx), 'XTickLabel'));
            set(axBode(idxAx), 'XTick', fTicks);
            if hasLabels
                % Labels rotated by 45 deg (horizontally, the 2 and the 0.01
                % of neighbouring axes overlap) and closer to the axis
                set(axBode(idxAx), 'XTickLabel', {'0.01', '0.3', '0.7', '1', '2'}, ...
                    'XTickLabelRotation', 45);
                axBode(idxAx).XAxis.TickLabelGapOffset = -2; % default 2 points
            else
                set(axBode(idxAx), 'XTickLabel', []);
            end
        end

        printFigure(gcf, figFolderStr, {'png'});
        printFigure(gcf, figFolderStrEps, {'eps'});


    end

end


end % compareLinearModels

function [M2, K2, Ce2, B1_3m, rb2, dG2dsched] = bladeMatricesMdl2(bladeVersion, wecs, M, K, Ce, Q3, B1_3, lam0, beta0, Ftbar)
% bladeMatricesMdl2 returns the Lagrange matrices of Mdl2 for the blade
% model bladeVersion at the operating point (lam0, beta0) with the thrust
% Ftbar (N):
% - 'original': M, K, Ce, B1_3 = [0; M\Q3] and r_b = wecs.rb unchanged
% - 'B1', 'B2' (variant B, bladeParamsMdl2 with the parameters wecs.blade*,
%   as the Simulink models Mdl_BianchiOL_BladeB*): modal blade mass at
%   r_ref, thrust load point r_b(lambda, beta), stiffness Kb(lambda, beta)
%   (B2: mass mb = Kb/wnb^2). dG2dsched: derivative of the generalized
%   blade force r_b*F_T - N*Kb*r_ref^2*theta with respect to (lambda0,
%   beta) at the equilibrium theta = F_T*r_b/(N*Kb*r_ref^2) (zero for
%   'original'). Mass and damping changes do not contribute (zero
%   accelerations and velocities at the equilibrium).
if strcmp(bladeVersion, 'original')
    M2 = M; K2 = K; Ce2 = Ce; B1_3m = B1_3; rb2 = wecs.rb; dG2dsched = [0, 0];
    return
end
N = wecs.N;
[rr, rb2, Kb, mb, Bb, dG2dsched] = bladeParamsMdl2(wecs, lam0, beta0, Ftbar);
M2 = M;
M2(1,1) = wecs.mt + N*mb;
M2(1,2) = N*mb*rr; M2(2,1) = M2(1,2);
M2(2,2) = N*mb*rr^2;
K2 = K; K2(2,2) = N*Kb*rr^2;
Ce2 = Ce; Ce2(2,2) = N*Bb*rr^2;
Q3m = Q3; Q3m(2,1) = rb2;
B1_3m = [zeros(4,3); M2\Q3m];
end
