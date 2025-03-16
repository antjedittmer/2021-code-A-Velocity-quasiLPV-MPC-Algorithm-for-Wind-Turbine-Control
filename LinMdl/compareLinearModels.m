function [sysOut,gapCell] = compareLinearModels(speedVec,figFolder, useActuatorStates,figNoAdd,createBodePlots,plotVisible,noOut)
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
%
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
    speedVec = [1,8,9,22]; % variations of wind speed
end

if nargin < 2 || isempty(figFolder)
    figFolder = fullfile(mainDir,'figDir');
end
% Create output folder
if ~isfolder(figFolder)
    mkdir(figFolder);
end

if nargin < 3, useActuatorStates = 0; end

if nargin < 4, figNoAdd = 0; close all; end

if nargin < 5, createBodePlots = 1; close all;  end

if nargin < 6,  plotVisible = 'on';  end

if nargin < 7, noOut = 3;  end


%% Load and initialize name of models and load into workspace
addpath(fullfile(mainDir,'dataIn'));
load('NREL5MW_CPdata','Rotor_Lamda', 'Rotor_Pitch', 'Rotor_cQ', 'Rotor_cT','Rotor_cP');

% Initialize name of models and load into workspace
idxModel = 1; %in case several linearized models are provided
modelNames = {'NREL5MW_linearised_4to25'};
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
    initModel5MWNREL(0, Rotor_Lamda, Rotor_Pitch, Rotor_cQ, Rotor_cT, Rotor_cP, figFolder);
wecs.Js = wecs.Jr + wecs.Ng^2 * wecs.Jg;

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


%% Loop over all air speeds in LinV
for index =  speedVec

    %% Variables dependent on wind speed V
    Vbar = Lin.V(index);
    Lambda_bar = Lin.RSpeed(index)*wecs.Rr/Vbar; % (-)
    omegabar = Lin.RSpeed(index); % 1/rad
    beta_bar = Lin.Pitch(index); % rad

    %% Grid points
    Cqb = interp2(lambda,pitch,Cq,Lambda_bar,beta_bar,'cubic');
    Ctb = interp2(lambda,pitch,Ct,Lambda_bar,beta_bar,'cubic');
    dCqdlambdab = interp2(lambda,pitch,Cqdlambda,Lambda_bar,beta_bar,'cubic');
    dCtdlambdab = interp2(lambda,pitch,Ctdlambda,Lambda_bar,beta_bar,'cubic');
    dCqdbetab = interp2(lambda,pitch,Cqdbeta,Lambda_bar,beta_bar,'cubic');
    dCtdbetab = interp2(lambda,pitch,Ctdbeta,Lambda_bar,beta_bar,'cubic');

    %% Tip speed ratio lambda = (rb * omega)/Ve =  (rb * omega)/(V - ydotfa)
    dlambdadomega = wecs.Rr/Vbar;
    dlambdadV = - omegabar*wecs.Rr/Vbar^2; %Vbar = Ve = V - ydot
    dlambdadydotfa =  - dlambdadV; %s/m   ydotfa m/s
    % dlambdaddelta =  - dlambdadV; %s/m   ydotfa m/s

    %% Ft = 0.5*rho*pi*rb^2*Ct*v^2
    kCT = 0.5*rho*pi*wecs.Rr^2; % constant for CT: 0.5*rho*A
    dFtdomega = kCT *dCtdlambdab*dlambdadomega*Vbar^2;
    dFtdV = kCT * Vbar * (dCtdlambdab*dlambdadV*Vbar + 2*Ctb); % dVe/dV = 1;
    dFtdydotfa = kCT * Vbar *(dCtdlambdab*dlambdadydotfa*Vbar - 2*Ctb); % dVe/dydot = -1;
    dFtdbeta = kCT*dCtdbetab*Vbar^2;

    %% Tr = 0.5*rho*pi*rb^3*Cq*v^2
    kCQ = kCT*wecs.Rr;
    dTrdomega = kCQ *dCqdlambdab*dlambdadomega*Vbar^2;
    dTrdV = kCQ *Vbar *(dCqdlambdab*dlambdadV*Vbar + 2*Cqb);
    dTrdydotfa = kCQ * Vbar*(dCqdlambdab*dlambdadydotfa*Vbar - 2*Cqb);
    dTrdbeta = kCQ *dCqdbetab*Vbar^2;

    %% Linear Model Mdl1: System, input matrices and state space model
    % States omega ydotfa xdotsw xfa xss

    % dotomega = (Tr - Ng * Tg)/J
    A11 = 1/wecs.Js * dTrdomega;
    A12 = 1/wecs.Js * dTrdydotfa;

    % ydotdotfa = 1/Mt( -Bt * ydotfa - Kt* yfa +Ft)
    wecs.mt = wecs.mtb;
    A21 = 1/wecs.mt * dFtdomega;
    A22 = 1/wecs.mt * dFtdydotfa;

    % xdotdotsw = 1/Mt( -Bt * xdotsw - Kt* xsw + 3/(2*H)Tg)
    A31 = ksw/wecs.mt * dTrdomega;
    A32 = ksw/wecs.mt * dTrdydotfa;

    A = [A11  A12                  0                 0                0;
        A21 (A22 - wecs.Bt/wecs.mt) 0                -wecs.Kt/wecs.mt  0;
        A31  A32                  -wecs.Bt/wecs.mt  0               -wecs.Kt/wecs.mt;
        [0 1 0 0 0; 0 0 1 0 0]];

    % Inputs: Torque [Nm], pitch angle beta [rad], wind speed V [m/s]
    B11 = -1/wecs.Js * wecs.Ng;
    B12 = 1/wecs.Js * dTrdbeta;
    B13 = 1/wecs.Js * dTrdV;

    B22 = 1/wecs.mtb * dFtdbeta;
    B23 = 1/wecs.mtb * dFtdV;

    B31 = -ksw/wecs.mt * wecs.Ng;
    B32 = ksw/wecs.mt * dTrdbeta;
    B33 = ksw/wecs.mt * dTrdV;

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

    C5aero =  [B1_3(5,1) * [-dFtdV, -dFtdV*wecs.rb, 0, dFtdomega 0]; ... % yFAdot
        B1_3(6,1) * [-dFtdV, -dFtdV*wecs.rb, 0, dFtdomega 0]; ... % deltaDot
        B1_3(7,2) *[-dTrdV, -dTrdV*wecs.rb,0,dTrdomega,0];... % xSWdot
        B1_3(8,2) *[-dTrdV, -dTrdV*wecs.rb,0,dTrdomega,0];... % omega_r
        zeros(1,5)];  % omega_gr
    QTildaTg = [0;0; B1_3(7,3);B11; -B11 ]; %B1_3(9,3);
    QTildaBeta = [B1_3(5,1) *dFtdbeta; B1_3(6,1)*dFtdbeta; B1_3(7,2) * dTrdbeta;...
        B1_3(8,2)*dTrdbeta;0];
    QTildaV = [B1_3(5,1) * dFtdV; B1_3(6,1) * dFtdV; B1_3(7,2) * dTrdV;...
        B1_3(8,2)*dTrdV;0];


    A9DoF = zeros(size(A1) + 2); % 9 states + 2 ctrl input
    A9DoF(1:4,:) = [A1(1:4,:) ,zeros(4,2)]; % states 1 - 4: y_FA,delta,x_sw, Theta_s
    A9DoF(5:9,1:4) = -M\K;
    A9DoF(5:9,5:9) = -M\Ce + C5aero; % states 5 - 9: y_FA,delta,x_sw, w_r, w_gr
    A9DoF(5:9,10:11) = [QTildaTg, QTildaBeta];
    %A9DoF(10:11,:) = zeros(2,11); % states 10 - 11: T_g, beta
    A9DoF(10,10) =  -1/kappa;
    A9DoF(11,11) =  -1/tau;


    B9DoF = zeros(size(A1,1) + 2,3); % 3 inputs -> 9 states + 2 ctrl input
    B9DoF(5:9,3) = QTildaV;
    B9DoF(10,1) = 1/kappa;
    B9DoF(11,2) = 1/tau;

    %     A9DoForig = [A1(1:4,:) ,zeros(4,2); % states 1 - 4: y_FA,delta,x_sw, Theta_s
    %         ...
    %         [A1(5,1:4), ... % state 5: xdot_FA (from  states 1 - 4: y_FA,delta,x_sw, Theta_s)
    %         A1(5,5) - B1(5,1)*dFtdV,... % state 5-> 5: xdot_FA <- xdot_FA
    %         A1(5,6) - B1(5,1)*dFtdV*wecs.rb, ...% state 6-> 5: xdot_FA <- delta
    %         A1(5,7), ... % state 7 > 5: xdot_FA <- xdot_sw
    %         A1(5,8) + B1(5,1)*dFtdomega, ... % state 8 > 5: xdot_FA <- omega_r
    %         A1(5,9), ...  % state 9 > 5: xdot_FA <- xdot_sw % state 10, 11 > 5: xdot_FA <- Tg, beta
    %         0 , B1(5,1)*dFtdbeta];
    %         ...
    %         [A1(6,1:4), ...% state 6: delta (from  states 1 - 4: y_FA,delta,x_sw, Theta_s)
    %         A1(6,5) - B1(6,1)*dFtdV,... % state 5> 6: delta<- xdot_FA
    %         A1(6,6) - B1(6,1)*dFtdV*wecs.rb,... % state 5> 6: delta<- xdot_FA
    %         A1(6,7), ... % delta<- xdot_FA
    %         A1(6,8) + B1(6,1)*dFtdomega, ... % state 5> 6: delta<- xdot_FA
    %         A1(6,9), ... % state 5> 6: delta<- xdot_FA
    %         0 ,B1(6,1)*dFtdbeta];
    %         ...
    %         [A1(7,1:4), ... % state 7: x_sw (from  states 1 - 4: y_FA,delta,x_sw, Theta_s)
    %         A1(7,5) - ksw * B1(7,2)*dTrdV,...
    %         A1(7,6) - ksw * B1(7,2)*dTrdV*wecs.rb,...
    %         A1(7,7), ...
    %         A1(7,8) + ksw * B1(7,2)*dTrdomega, ...
    %         A1(7,9),...
    %         B31, ksw * B1(7,2)*dTrdbeta]; %  B31 = wecs.Ng/wecs.mtb *3/(2*wecs.H); Tg -> x_sw
    %         ...
    %         [A1(8,1:4), ... % state 8: omega_r (from  states 1 - 4: y_FA,delta,x_sw, Theta_s)
    %         A1(8,5) - B1(8,3)*dTrdV,...
    %         A1(8,6) - B1(8,3)*dTrdV*wecs.rb, ....
    %         A1(8,7), ...
    %         A1(8,8) + B1(8,3)*dTrdomega, ...
    %         A1(8,9), ...
    %         B11, B1(8,3)*dTrdbeta]; %  B11 = - wecs.Ng/wecs.Jr; Tg -> omega_r (instead of B1(8,1))
    %         ...
    %         [A1(9,1:6), A1(9,7) , A1(9,8:9) , -B11 ,0]; %B1(9,4)
    %         zeros(1,9),  -1/kappa, 0;...
    %         zeros(1,10),  -1/tau];
    %
    A7DoF = A9DoF(1:9,1:9);
    %
    %     B9DoF = [zeros(4,3);... % states 1 - 4: y_FA,delta,x_sw, Theta_s
    %         zeros(1,2),B1(5,1)*dFtdV;... % state 5: ydot_FA
    %         zeros(1,2),B1(6,1)*dFtdV;...% state 6: xdot_FA
    %         B31,0,B1(7,2)*ksw *dTrdV;...% state 7: xdot_SW
    %         B11,0,B1(8,3)*dTrdV;... state 8: omega_r
    %         zeros(1,3);... state 9: omega_gr
    %         1/kappa, 0,0;...
    %         0,1/tau,0];

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

    Cout = [C9DoF;
        [zeros(1,7+2) 1 zeros(1,1)];  % beta
        [zeros(1,8+2) 1]];% torque
    Dout = [D9DoF; zeros(2,3)];
    sysOut = ss(A9DoF,B9DoF*diag(multFASTInput0),Cout,Dout); % to give input Tg in kN/m

    sysOut.Outputname = [sysOutputname;{'beta'};{'Tg'}];
    sysOut.Inputname = sysInputname;

    sysCell{index} = sysOut;
    Torque_g_bar = 0.5*rho*pi*wecs.rb^3*Cqb*Vbar^2/wecs.Ng; % Nm
    u_bar_cell{index} = [Torque_g_bar, beta_bar * pi/180,Vbar]';
    x_bar_cell{index}  = - sysLagrange.A\sysLagrange.B *u_bar_cell{index}; % X(s) = (sI - A)1 *B *U(s) (xdot = 0 = Ax +Bu)
    y_bar_cell{index}  = sysLagrange.C * x_bar_cell{index} + sysLagrange.D * u_bar_cell{index};

    %% Make a Bode plot
    if createBodePlots == 1
        ff = figure(index*100+3 + figNoAdd);
        ff.Position = [520   267   560   531];
        ff.Color = 'white'; ff.Visible = plotVisible;
        %wVec = logspace(opts.XLim{1}(1),opts.XLim{1}(2),100);

        bodemag(sysLagrange,opts); hold on; grid on;
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

        
        % Title string
        titleStr = sprintf('Bode mag. plot for V = %d, pitch = %2.2f, and TSR = %2.4f',Vbar,beta_bar,Lambda_bar);
        titleCell = {titleStr, [legCell{1},legCell{2},strExtr]};
        title(['V_{\infty}',sprintf(' = %d: %s', Vbar, titleCell{2})]); %title(titleCell)

        % Figure information
        figStr = [strrep(titleStr,'Bode mag. plot for V =','Figure_Bode V ='),'.png'];
        figFolderStr = fullfile(figFolder,figStr);

        figFolderStrEps = fullfile(figFolder,sprintf('BodeV%d',Vbar));
        posaxes =   get(0,'defaultFigurePosition');
        set(gcf,'Position',[posaxes(1:3),posaxes(4)*1.1]);

        print(figFolderStr, '-dpng');
        print(figFolderStrEps, '-depsc');

        
    end

end

