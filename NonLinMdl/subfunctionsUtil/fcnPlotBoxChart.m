function fcnPlotBoxChart(tableObj,nameVec,vecRateLim,fs,ylabelVec,figNo,tickVec,cl,noText)

tableObj0 = tableObj;
vecBl = (100 - tableObj.bladeIPC)/100;
vecTwr = (100 - tableObj.twrIPC)/100;


if nargin <= 5
    figNo = 100;
end

fsl = fs -2;
if nargin <= 6
    legCell = regexp(sprintf('Rate lim. %d°/s,', vecRateLim),',','split');
    legCell0 = legCell(1:end-1);
    tickVec{1} = {'','',''};
    tickVec{2} = legCell0;
    fsl = fs;

end

if nargin <=7 || isempty(cl)
     cl= [0*ones(1,3); lines]; %cl = lines;
end

if nargin <= 8
    noText = 0;
end

figure(figNo);
tiledlayout(length(nameVec),1,"TileSpacing","compact")

for idx = 1: length(nameVec)
    data = tableObj0.(nameVec{idx})';

    yLimVec = max(data(:));

    ax1(idx) = nexttile; %subplot(2,1,idx); %#ok<SAGROW>
    hold on;
    for i = 1:size(data, 2)

        aDataAll = data(:, i);
        aData = aDataAll(~isnan(aDataAll));
        boxchart(ones(size(aData, 1), 1) * i, aData, ...
            'BoxFaceColor', cl(i, :), ...
            'MarkerStyle', '+',...
            'MarkerColor', cl(i, :));  % Add this line

        % Compute statisticsaData
        medianValue = median(aData);
        meanValue = mean(aData);
        maxValue = max(aData);
        maxDisp = max(maxValue,3.25);


        dataRange = max(data(:, i)) - min(data(:, i));
        offsetMax = 0.05 * dataRange;   % above max
        offsetMed = 0.025 * dataRange;  % above median, but well below max


        % Annotate maximum value near the max point
        if maxValue < 10
            strMax = sprintf('Max: %.1f', maxValue);
            strMed = sprintf('Med: %.1f', medianValue);
    
        else
            strMax = sprintf('Max: %d', round(maxValue));
            strMed = sprintf('Med: %d', round(medianValue));

        end
        
       if noText == 0

        text(i + 0.075, maxDisp, strMax, ...
            'HorizontalAlignment', 'left', 'FontSize', fs-1, 'Color', 'r');
        % Annotate maximum value near the max point

        if maxValue < 1000
        text(i + 0.075, maxDisp - yLimVec/(4*size(data, 2)), strMed, ...
            'HorizontalAlignment', 'left', 'FontSize', fs-1, 'Color', 'k');

        elseif contains(nameVec,'power')
              text(i + 0.075, maxDisp - yLimVec/(8*(size(data, 2))), strMed, ...
            'HorizontalAlignment', 'left', 'FontSize', fs-1, 'Color', 'k');
        else
             text(i + 0.075, maxDisp - yLimVec/(size(data, 2)), strMed, ...
            'HorizontalAlignment', 'left', 'FontSize', fs-1, 'Color', 'k');
        end

       end

        % Annotate max value above whisker
        % text(i + 0.075, maxValue + offsetMax, strMax, ...
        %     'HorizontalAlignment', 'left', 'FontSize', fs-1, 'Color', 'r');
        %
        % % Annotate median value below the box
        % text(i + 0.075, medianValue + offsetMed, strMed, ...
        %     'HorizontalAlignment', 'left', 'FontSize', fs-1, 'Color', 'k');

    end
    hold off;
    axis tight; grid on;
    yLimVec = get(ax1(idx),"YLim");
    xLimVec = get(ax1(idx),"XLim");
    if length(nameVec) > 1
        set(ax1(idx),"YLim",[yLimVec(1) yLimVec(2)*1.1]);
    else
        set(ax1(idx),"YLim",[yLimVec(1)*0.999 yLimVec(2)*1.001]);
    end
    ylabel(ylabelVec{idx});
    set(ax1(idx),"XLim",[xLimVec(1) xLimVec(2)*1.05]);
    set(gca,'FontSize',fs);
    set(gca,'XTick',1:size(data, 2),'XTickLabel',tickVec{idx})
    set(gca,'XTickLabelRotation',0)
end


%legend(legCell0,'Location','northwest')

pos0 = get(0,'defaultFigurePosition');
set(gcf,'Position', [pos0(1:3),pos0(4)*1])
