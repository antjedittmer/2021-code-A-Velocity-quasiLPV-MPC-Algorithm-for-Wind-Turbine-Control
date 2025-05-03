function aMatrix = plotNormBodePlots(gapCell,speedVec,figFolder,createNugapPlot)
% plotNormBodePlots plots the norm of the difference of the signals in the
% bode plots.
% Inputs gapCell,speedVec,figFolder

%% Handle inputs

if nargin < 2 || isempty(speedVec)
    speedVec = [1,8,9,22]; % variations of wind speed
end

if nargin < 3 || isempty(figFolder)
    % Get relative paths to main directory (for data input dir)
    workDir = fileparts(mfilename('fullpath'));
    mainDir = fileparts(workDir);
    figFolder = fullfile(mainDir,'figDir');
end

if ~isfolder(figFolder), mkdir(figFolder); end

if nargin < 4 || isempty(figFolder)
    createNugapPlot = 0;
end

%% Create matrix with values for bar plots
fieldsGapAll = fieldnames(gapCell); % get the field names
fieldsGap = fieldsGapAll(contains(fieldsGapAll,'n')); % only use Vinnicombe/nu-gap and norm for now
sortFieldsGap = sort(fieldsGap); % Mdl1 (5 states) should be listed before Mdl2 (9 states)

cnt = 1; % this is the counter for constructing
aMatrix = nan(4*length(speedVec), 9); %the matrix with the values

for idxGapN = 1:length(fieldsGap)
    aFieldName = sortFieldsGap{idxGapN};
    aCellField = gapCell.(aFieldName);
    fprintf('\n %s\n---\n',aFieldName)

    for idx = speedVec
        str1 = sprintf('%2.2f & ',aCellField{idx});
        fprintf('%d m/s & %s \n',idx+3, str1); % simply write the tables out
        tmp = aCellField{idx}'; %Transpose: Matlab is column first
        aMatrix(cnt,:) = tmp(:); % but the subplots are row first
        cnt = cnt + 1; % next row of the matrix
    end
end


%% Set informations for both plots

lc = lines; %line colors
lc1 = [0*[0.5, 0.5,0.5];lc];
len = length(speedVec);
vec = 1:len;
idxVec = 1:9; % this is the number of input/output combination

tmpStr = sprintf('%d m/s,', speedVec + 3); %wind speed infor
windTickLabel = regexp(tmpStr(1:end-1), ',','split');

% sysOutputname =[{'RotSpd \omega_r (rad/s)'};{'TwrAcc_{fa} (m/s^2)'}; {'TwrAcc_{sw} (m/s^2)'}];
% sysInputname =[{'GenTq T_g (kNm)'}; {'BlPitch \beta_0 (rad)'};{'Wind V_{\infty} (m/s)'}];

sysOutputname =[{'RotSpd \omega_r (-)'};{'TwrAcc_{fa} (-)'}; {'TwrAcc_{sw} (-)'}];
sysInputname =[{'GenTq T_g (-)'}; {'BlPitch \beta_0 (-)'};{'Wind V_{\infty} (-)'}];

Mdl2Str = ' Rot,Twr,Gen+Bld';

%% Set values for first, norm plot
aMat1 = [aMatrix(vec,idxVec),aMatrix(vec+len,idxVec)];
aMax = max(aMat1(:))* 1.2;
legCell = ['                        Norm: \color{black}Mdl1: Rot+Twr ',...
    '\color[rgb]{',num2str(lc(1,:)),'}Mdl2:',Mdl2Str];

scaleF = 100;
figure(1);
pos0 = get(0,'defaultFigurePosition');
set(gcf,'Position',[pos0(1)-pos0(3)*0.7,pos0(2),pos0(3)*1.5,pos0(4)])
t = tiledlayout(3,3);

t.TileSpacing = 'compact';
t.Padding = 'compact';

for idx = 1:9
    %subplot(3,3,idx)
    nexttile;
    b = bar([aMatrix(vec,idx),aMatrix(vec+len,idx)]*scaleF);

    for k = 1:length(b), b(k).FaceColor = lc1(k,:); end
    set(gca,'XTickLabel',windTickLabel);
    set(gca,'YLim',[0,aMax*scaleF]);

    if idx == 1
        title( legCell)
    end
    % if printTitle
    %     title(titleCellNew)
    % else

    % end
for idxB = 1: length(b)
    xtips1 = b(idxB).XEndPoints;
    ytips1 = b(idxB).YEndPoints;
    yData = b(idxB).YData;

    if idxB <= 5
        labels1 = string(round(yData,1));
    else
        labels1 = string(round(yData*10,2));
    end

    % Adjustments
    yOffset = 0.02 * aMax * scaleF;  % Move text up a bit more
    if idxB == 1
        xOffset = 0.; % Left bar, shift slightly left
    else
        xOffset = +0.15; % Right bar, shift slightly right
    end

    % Apply text with rotation and offset
    text(xtips1 + xOffset, ytips1 + yOffset, labels1, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'Rotation', 30, ...
        'FontSize', 9)
end


    if mod(idx-1,3)== 0
        anOutputStr = sysOutputname{(idx-1)/3 + 1};
        ylabel(anOutputStr);
    end

    if idx >= 7
        anInputStr = sysInputname{idx-6};
        xlabel(anInputStr);
    end
end

figStr = 'Norm';
figFolderStr = fullfile(figFolder,figStr);
figFolderStrEps = figFolderStr;
print(figFolderStr, '-dpng');
print(figFolderStrEps, '-depsc');

if ~createNugapPlot
    return;
end


%% Set values for second, nugap plot
aMat1 = [aMatrix(vec+2*len,idxVec),aMatrix(vec+3*len,idxVec)];
aMax = max(aMat1(:))* 1.01;
legCell = ['                          Nugap: \color{black}Mdl1: Rot+Twr ',...
    '\color[rgb]{',num2str(lc(1,:)),'}Mdl2:',Mdl2Str];

figure(2);
for idx = 1:9
    subplot(3,3,idx)
    b = bar([aMatrix(vec+ 2*len,idx),aMatrix(vec+3*len,idx)]);

    for k = 1:length(b), b(k).FaceColor = lc1(k,:); end
    set(gca,'XTickLabel',windTickLabel);
    set(gca,'YLim',[0,aMax]);

    if idx == 1
        title( legCell)
    end

    if mod(idx-1,3)== 0
        anOutputStr = sysOutputname{(idx-1)/3 + 1};
        ylabel(anOutputStr);
    end

    if idx >= 7
        anInputStr = sysInputname{idx-6};
        xlabel(anInputStr);
    end
end

figStr = 'Nugap';
figFolderStr = fullfile(figFolder,figStr);
figFolderStrEps = figFolderStr;
print(figFolderStr, '-dpng');
print(figFolderStrEps, '-depsc');