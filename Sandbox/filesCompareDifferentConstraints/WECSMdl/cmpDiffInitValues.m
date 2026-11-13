close all;

% Set path to input data directory
workDir = fileparts(mfilename('fullpath'));
mainDir = fileparts(workDir);

% Data vor cQ and cT LUT
dataInDir = fullfile(mainDir,'dataIn');
load(fullfile(dataInDir,'NREL5MW_CPdata.mat'),...
    'Rotor_Lamda','Rotor_Pitch','Rotor_cQ','Rotor_cT');


plotOn = 0;
Rotor_Lamda = ''; Rotor_Pitch = ''; Rotor_cQ = ''; Rotor_cT = '';
figDir = ''; titleOn = 0;

multRb = 0.93;
[wecs, M, Ce, K, Q, L, rho, tau, kappa, lambda, beta,Cq,Ct,Q3] = ....
    initModel5MWNREL(plotOn, Rotor_Lamda, Rotor_Pitch, Rotor_cQ, Rotor_cT, figDir,titleOn,multRb);

A1 = [zeros(4), L;...
    -M\K, -M\Ce];
eigA1 = eig(A1);

multRb = 1;
[wecs, M, Ce, K, Q, L, rho, tau, kappa, lambda, beta,Cq,Ct,Q3] = ....
    initModel5MWNREL(plotOn, Rotor_Lamda, Rotor_Pitch, Rotor_cQ, Rotor_cT, figDir,titleOn,multRb);


A2 = [zeros(4), L;...
    -M\K, -M\Ce];
eigA2 = eig(A2);

figure; plot(real(eigA1), imag(eigA1),'k*',real(eigA2), imag(eigA2),'ro' )





