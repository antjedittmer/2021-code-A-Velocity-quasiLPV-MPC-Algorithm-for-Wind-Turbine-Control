function setTimeLabelFontSize(fig)
% setTimeLabelFontSize sets the font size of the time axis labels of the
% figure fig (x labels of axes and tiled layouts that start with 'Time') to
% the root application data TimeLabelFontSize if set
% (runGenerateAllThesisFigures, timeLabelFontSize), otherwise 14. The other
% texts keep their size (getTimeSeriesFontSize). Call it after the font
% sizes of the figure are set.
if nargin < 1 || isempty(fig)
    fig = gcf;
end
fs = 14;
if isappdata(groot, 'TimeLabelFontSize') && ~isempty(getappdata(groot, 'TimeLabelFontSize'))
    fs = getappdata(groot, 'TimeLabelFontSize');
end
h = [findall(fig, 'Type', 'axes'); findall(fig, 'Type', 'tiledlayout')];
for k = 1:numel(h)
    lbl = h(k).XLabel;
    if startsWith(strtrim(strjoin(cellstr(lbl.String), ' ')), 'Time')
        lbl.FontSize = fs;
    end
end
end
