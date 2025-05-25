function plotNormTimePlots(normStruct,figNo1,figDirStr)
%plotNormTimePlots plots the norm of the difference of the signals in the
% time plots.

if nargin <2 || isempty(figNo1)
    figNo1 = 100;
end

if nargin < 3
    figDirStr = 'figDir';
end


fieldsGapAll = fieldnames(normStruct); %get the field names
sortFieldsGap = fieldsGapAll; % Mdl1 (5 states) should be listed before Mdl2 (9 states)

sortFieldsSignals =  [{'GenTq'}  {'BlPitch1'}  {'RotSpeed'}    {'GenPwr'}   {'TwrAccFA'}    {'TwrAccSW'}];

% sysOutputname =[{'RotSpd \omega_r (-)'};{'TwrAcc_{fa} (-)'}; {'TwrAcc_{sw} (-)'}];
% sysInputname =[{'GenTq T_g (-)'}; {'BlPitch \beta_0 (-)'};{'Wind V_{\infty} (-)'}];

% NcIMUTAxs
%windTickLabel = sortFieldsGap;


%% Set informations for both plots

lc = lines; %line colors
lc1 = [0*[0.5, 0.5,0.5];lc];
%len = 6;
%vec = 1:len;
%idxVec = 1:len; % this is the number of input/output combination


%% Set values for first, norm plot


Mdl2Str = ' Rot,Twr,Gen+Bld';
legCell = ['                        Norm: \color{black}Mdl1: Rot+Twr ',...
    '\color[rgb]{',num2str(lc(1,:)),'}Mdl2:',Mdl2Str];

figure(figNo1);
for idx = 1:6
    subplot(3,2,idx)
    aFieldname = sortFieldsSignals{idx};
    aFieldn = aFieldname;
    if strcmp(aFieldn,'TwrAccFA')
        aFieldn = 'NcIMUTAxs';
    elseif strcmp(aFieldn,'TwrAccSW')
         aFieldn = 'NcIMUTAys';
    end

    b = bar([normStruct.Sweep.(aFieldn); normStruct.NTM18.(aFieldn)]);

    for k = 1:length(b), b(k).FaceColor = lc1(k,:); end
    set(gca,'XTickLabel',sortFieldsGap);
    %set(gca,'YLim',[0,aMax]);

    if idx == 1
        title( legCell)
    end
% if printTitle
%     title(titleCellNew)
% else
%     for idx = 1: length(b)
%         xtips1 = b(idx).XEndPoints;
%         ytips1 = b(idx).YEndPoints;
%         if idx <= 5
%         labels1 = string(round(b(idx).YData,2));
%         else
%             labels1 = string(round(b(idx).YData*10,2));
%         end
%         text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
%             'VerticalAlignment','bottom')
%     end
% end

    ylabel([aFieldname,' (-)']);

end

figFolder = figDirStr;
figStr = 'NormTime';
figFolderStr = fullfile(figFolder,figStr);
figFolderStrEps = figFolderStr;
print(figFolderStr, '-dpng');
print(figFolderStrEps, '-depsc');