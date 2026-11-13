
function qWeight = getInfluenceBladeParams(x,freqDomain)

if ~nargin
    x = [1,1];
end

if nargin <2
    freqDomain = 1;
end


if freqDomain

    %x0 = [1,1];


    % fun = @getInfluenceBladeParams;
    % [x,fval,exitflag,output] = fminsearch(fun,x0,options);
    %
    % save('x','x');

    %x =  [0.608380, 7128.534769];

    %x =   [4.1814    0.7651];

    % Get relative paths to main directory (for data input dir)
    workDir = mfilename('fullpath');
    mainDir = workDir;

    speedVec = [1,8,9,22]; % [1,22]; % variations of wind speed 8,9,

    figFolder = fullfile(mainDir,'figDir');

    useActuatorStates = 0;

    figNoAdd = 0;
    createBodePlots = 1;
    plotVisible = 'off';
    noOut = 3;

    createNugapPlot = 0;

    [~,gapCell] = compareLinearModels(speedVec,figFolder,useActuatorStates,figNoAdd,createBodePlots,plotVisible,noOut,x);
    aMatrix = plotNormBodePlots(gapCell,speedVec,figFolder,createNugapPlot);

    qWeight = norm(aMatrix(5:8,1:6));


else
    % [wecs, M, Ce, K, Q, L, rho, tau, kappa, lambda, pitch, Cq, Ct, Q3, Cp] = initModel5MWNRELTune(x);
    %updateDDMdl1(0.75);

    multRb = 0.75;
    useTurn = 1;
    updateDDMdl1(multRb,useTurn, x)

    DT = 0.008; %#ok<NASGU>

    try
        normStruct = runCompareModel2;
        normArray = struct2array(normStruct);
        qWeight1 = normArray([1,5,6]) * [100,1,1]';
    catch
        qWeight1 = 10;
    end

    try
        normStruct2 = runCompareModel2('Sweep');
        normArray2 = struct2array(normStruct2);
        qWeight2 = normArray2([1,5,6]) * [100,1,1]';
    catch
        qWeight2 = 10;
    end


    qWeight = norm([qWeight1, qWeight2]);
end


% 0.3540  6.9975e+03

