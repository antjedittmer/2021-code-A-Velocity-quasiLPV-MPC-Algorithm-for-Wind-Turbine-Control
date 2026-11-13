function fcnPlotBar(K_matrix,strK,summaryTickLabel,printTitle,figNo)
% fcnPlotBar plots the bar plot for the values given in the K_matrix.
% fcnPlotBar(K_matrix,strK,summaryTickLabel,printTitle,figNo)


if nargin < 4
    printTitle = 0;
end

if nargin < 5
    figNo = 1;
end


figure(figNo);
cl1 = colormap('lines');
cl = [0.5,0.5,0.5;
    cl1(4,:);
    [1,0,0];
    cl1(1,:);
    cl1(5,:)];

b = bar(K_matrix);
legend('C_{P}','C_{\sigma(P)}','C_{Actuator}','C_{D,Blades}','C_{D,Tower}','Location','NorthWest');

for k = 1:5
    b(k).FaceColor = cl(k,:);
end

axis tight; grid on;
titleCellNew = {['Power tracking C_{P}: ', '{\color[rgb]{',num2str(cl(1,:)),'}',strK.power(1:end-2),'}'],...
    ['Power variance C_{\sigma(P)}: ',  '{\color[rgb]{',num2str(cl(2,:)),'}',strK.StdPowSpeed(1:end-2),'}'], ...
    ['Actuator power C_{Actuator}: ',  '{\color[rgb]{',num2str(cl(3,:)),'}',strK.ActPwr(1:end-2),'}'], ...
    ['Blade damage C_{D,Blades}: ',  '{\color[rgb]{',num2str(cl(4,:)),'}',strK.DamageBlades(1:end-2),'}'],...
    ['Tower damage C_{D,Tower}: ',  '{\color[rgb]{',num2str(cl(5,:)),'}',strK.DamageTower(1:end-2),'}']};

if printTitle
    title(titleCellNew)
else
    for idx = 1: length(b)
        xtips1 = b(idx).XEndPoints;
        ytips1 = b(idx).YEndPoints;
        labels1 = string(round(b(idx).YData,1));
        text(xtips1,ytips1,labels1,'HorizontalAlignment','center',...
            'VerticalAlignment','bottom')
    end
end

ax = gca;
ax.XLim = ax.XLim + [-0.1 0.1];
ax.YLim = ax.YLim + [0 0.2];
ax.XTickLabelRotation = 22.5;

ax.XTickLabel = summaryTickLabel;