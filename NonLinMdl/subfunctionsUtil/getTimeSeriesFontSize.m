function fs = getTimeSeriesFontSize(fsDefault)
% getTimeSeriesFontSize returns the font size of the time series figures:
% the root application data TimeSeriesFontSize if set (runGenerateAllThesisFigures,
% timeSeriesFontSize), otherwise fsDefault (default 13).
if nargin < 1 || isempty(fsDefault)
    fsDefault = 13;
end
fs = fsDefault;
if isappdata(groot, 'TimeSeriesFontSize') && ~isempty(getappdata(groot, 'TimeSeriesFontSize'))
    fs = getappdata(groot, 'TimeSeriesFontSize');
end
end
