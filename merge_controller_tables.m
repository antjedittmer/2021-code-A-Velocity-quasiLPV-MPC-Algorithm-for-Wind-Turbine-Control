function merge_controller_tables(file_sweep, file_ntw, file_out, file_ntw_filt, useFilteredNtw)
% merge_controller_tables merges the table rows written by runCompareCtrl
% for the stepped wind sweep and the turbulent wind into one LaTeX table.
%
% merge_controller_tables('results_tableSweep.tex', ...
%                         'results_tableNTW18.tex', ...
%                         'results_tableMerged.tex', ...
%                         'results_tableNTW18_windFilt.tex')
%
% file_ntw_filt (optional): rows of the turbulent wind with the qLMPC also
% on the filtered plant and scheduled with the filtered wind (runCompareCtrl
% with useFiltered = 2), added as a third column.
% useFilteredNtw (optional, default 1): useFiltered of the turbulent wind
% column (runCompareCtrl), for the caption: 1 PI on the filtered plant,
% 0 both controllers on the same plant without filter.
%
% Expected row format of the input files (written by runCompareCtrl):
%   Tower fore-aft acc. (m/s$^2$) & $0.0514 / 0.0319 = 1.61$ \\
% The cells are copied unchanged, so the number formats set in
% runCompareCtrl (four digits for the tower acceleration, two otherwise)
% are kept. Rows are matched by their label, not by their position.

if nargin < 1 || isempty(file_sweep), file_sweep = 'results_tableSweep.tex';  end
if nargin < 2 || isempty(file_ntw),   file_ntw   = 'results_tableNTW18.tex';  end
if nargin < 3 || isempty(file_out),   file_out   = 'results_tableMerged.tex'; end
if nargin < 4, file_ntw_filt = ''; end
if nargin < 5 || isempty(useFilteredNtw), useFilteredNtw = 1; end
hasFilt = ~isempty(file_ntw_filt);

[labels_sw,  cells_sw]  = parseRows(file_sweep);
[labels_ntw, cells_ntw] = parseRows(file_ntw);
if hasFilt
    [labels_filt, cells_filt] = parseRows(file_ntw_filt);
end

fid = fopen(file_out, 'w');
if fid < 0
    error('Could not open %s for writing.', file_out);
end
cleanupFid = onCleanup(@() fclose(fid));

fprintf(fid, '\\begin{table}[t]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption[Standard deviation comparison: Baseline PI vs.\\ qLMPC for stepped wind sweep and turbulent wind.]%%\n');
fprintf(fid, '{Standard deviation comparison: Baseline PI vs.\\ qLMPC for stepped wind sweep and\n');
fprintf(fid, 'turbulent wind (18\\,m/s mean, IEC class B).\n');
fprintf(fid, 'The standard deviations obtained with the PI controller are divided by those obtained with the qLMPC;\n');
if hasFilt
    fprintf(fid, 'quotients above~1 indicate improvement. Ratios are computed from unrounded values.\n');
    if useFilteredNtw == 0
        fprintf(fid, 'Stepped wind sweep and turbulent wind: both controllers on the same plant without the\n');
        fprintf(fid, 'rotor-effective wind filter, qLMPC scheduled with the unfiltered wind; filtered: both\n');
        fprintf(fid, 'controllers on the plant with the filter, qLMPC scheduled with the filtered wind.}\n');
    else
        fprintf(fid, 'Turbulent wind: PI on the plant with the rotor-effective wind filter, qLMPC on the plant\n');
        fprintf(fid, 'without filter and scheduled with the unfiltered wind; filtered: qLMPC also on the plant\n');
        fprintf(fid, 'with the filter and scheduled with the filtered wind.}\n');
    end
    fprintf(fid, '\\label{tab:ControllerComparisonMerged}\n');
    fprintf(fid, '\\begin{tabular}{l rrr}\n');
    fprintf(fid, '\\toprule\n');
    fprintf(fid, '{Standard Deviation of} & {Stepped wind sweep} & {Turbulent wind} & {Turbulent wind, filtered} \\\\\n');
    fprintf(fid, '{Signal (unit)} & \\multicolumn{3}{c}{std PI / std qLMPC = ratio} \\\\\n');
else
    fprintf(fid, 'quotients above~1 indicate improvement. Ratios are computed from unrounded values.}\n');
    fprintf(fid, '\\label{tab:ControllerComparisonMerged}\n');
    fprintf(fid, '\\begin{tabular}{l rr}\n');
    fprintf(fid, '\\toprule\n');
    fprintf(fid, '{Standard Deviation of} & {Stepped wind sweep} & {Turbulent wind} \\\\\n');
    fprintf(fid, '{Signal (unit)} & \\multicolumn{2}{c}{std PI / std qLMPC = ratio} \\\\\n');
end
fprintf(fid, '\\midrule\n');

for i = 1:numel(labels_sw)
    idx = find(strcmp(labels_ntw, labels_sw{i}), 1);
    if isempty(idx)
        error('Row "%s" not found in %s.', labels_sw{i}, file_ntw);
    end
    if hasFilt
        idxF = find(strcmp(labels_filt, labels_sw{i}), 1);
        if isempty(idxF)
            error('Row "%s" not found in %s.', labels_sw{i}, file_ntw_filt);
        end
        fprintf(fid, '%s & %s & %s & %s \\\\\n', labels_sw{i}, cells_sw{i}, cells_ntw{idx}, cells_filt{idxF});
    else
        fprintf(fid, '%s & %s & %s \\\\\n', labels_sw{i}, cells_sw{i}, cells_ntw{idx});
    end
end

fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');

fprintf('Merged table written to %s\n', file_out);
end

% =========================================================================
% Local function
% =========================================================================
function [labels, cells] = parseRows(filename)
% parseRows returns the labels and cells of the rows 'label & cell \\'.
if ~isfile(filename)
    error('File %s not found.', filename);
end
fileLines = splitlines(fileread(filename));
labels = {};
cells  = {};
for i = 1:numel(fileLines)
    lineStr = strtrim(fileLines{i});
    if isempty(lineStr) || startsWith(lineStr, '%') || ~contains(lineStr, '&')
        continue % skip empty lines, LaTeX comments and non-data lines
    end
    parts = strtrim(split(lineStr, '&'));
    if numel(parts) ~= 2
        error('Unexpected row format in %s: %s', filename, lineStr);
    end
    labels{end+1} = parts{1};                                   %#ok<AGROW>
    cells{end+1}  = strtrim(regexprep(parts{2}, '\\\\\s*$', '')); %#ok<AGROW> remove trailing \\
end
if isempty(labels)
    error('No table rows found in %s.', filename);
end
end