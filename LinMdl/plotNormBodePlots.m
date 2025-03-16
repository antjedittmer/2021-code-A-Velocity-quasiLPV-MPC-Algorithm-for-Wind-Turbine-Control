function plotNormBodePlots(gapCell,speedVec,figFolder)
% plotNormBodePlots plots the norm of the difference of the signals in the
% bode plots.
% Inputs 

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
aMax = max(aMat1(:))* 1.01;
legCell = ['                        Norm: \color{black}Mdl1: Rot+Twr ',...
    '\color[rgb]{',num2str(lc(1,:)),'}Mdl2:',Mdl2Str];

figure(1);
for idx = 1:9
    subplot(3,3,idx)
    b = bar([aMatrix(vec,idx),aMatrix(vec+len,idx)]);

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

figStr = 'Norm';
figFolderStr = fullfile(figFolder,figStr);
figFolderStrEps = figFolderStr;
print(figFolderStr, '-dpng');
print(figFolderStrEps, '-depsc');


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