function merge_controller_tables(file_sweep, file_ntw, file_out)
% merge_controller_tables merges the table rows written by runCompareCtrl
% for the stepped wind sweep and the turbulent wind into one LaTeX table.
%
% merge_controller_tables('results_tableSweep.tex', ...
%                         'results_tableNTW18.tex', ...
%                         'results_tableMerged.tex')
%
% Expected row format of the input files (written by runCompareCtrl):
%   Tower fore-aft acc. (m/s$^2$) & $0.0514 / 0.0319 = 1.61$ \\
% The cells are copied unchanged, so the number formats set in
% runCompareCtrl (four digits for the tower acceleration, two otherwise)
% are kept. Rows are matched by their label, not by their position.

if nargin < 1 || isempty(file_sweep), file_sweep = 'results_tableSweep.tex';  end
if nargin < 2 || isempty(file_ntw),   file_ntw   = 'results_tableNTW18.tex';  end
if nargin < 3 || isempty(file_out),   file_out   = 'results_tableMerged.tex'; end

[labels_sw,  cells_sw]  = parseRows(file_sweep);
[labels_ntw, cells_ntw] = parseRows(file_ntw);

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
fprintf(fid, 'quotients above~1 indicate improvement. Ratios are computed from unrounded values.}\n');
fprintf(fid, '\\label{tab:ControllerComparisonMerged}\n');
fprintf(fid, '\\begin{tabular}{l rr}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, '{Standard Deviation of} & {Stepped wind sweep} & {Turbulent wind} \\\\\n');
fprintf(fid, '{Signal (unit)} & \\multicolumn{2}{c}{std PI / std qLMPC = ratio} \\\\\n');
fprintf(fid, '\\midrule\n');

for i = 1:numel(labels_sw)
    idx = find(strcmp(labels_ntw, labels_sw{i}), 1);
    if isempty(idx)
        error('Row "%s" not found in %s.', labels_sw{i}, file_ntw);
    end
    fprintf(fid, '%s & %s & %s \\\\\n', labels_sw{i}, cells_sw{i}, cells_ntw{idx});
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