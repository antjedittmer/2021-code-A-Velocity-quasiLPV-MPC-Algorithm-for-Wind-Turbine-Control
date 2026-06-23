function merge_controller_tables(file_sweep, file_ntw, file_out)
% merge_controller_tables('results_tableSweep.tex',
%                         'results_tableNTW18.tex',
%                         'results_tableMerged.tex')
if ~nargin
    file_sweep = 'results_tableSweep.tex';
    file_ntw   = 'results_tableNTW18.tex';
    file_out   = 'results_tableMerged.tex';
end

% --- Parse a table file, returns cell array of rows {metric, unit, stdPI, stdqLPV, ratio}
function rows = parseTable(filename)
    fid  = fopen(filename, 'r');
    rows = {};
    while ~feof(fid)
        line = strtrim(fgetl(fid));
        % Data rows contain & and \\ but not { (header/toprule lines do)
        if contains(line, '&') && contains(line, '\\') && ~contains(line, '{$')
            parts = strsplit(line, '&');
            if numel(parts) == 5
                metric = strtrim(parts{1});
                unit   = strtrim(parts{2});
                stdPI  = str2double(strtrim(parts{3}));
                stdqLP = str2double(strtrim(parts{4}));
                ratio  = str2double(strtrim(strrep(parts{5}, '\\', '')));
                rows{end+1} = {metric, unit, stdPI, stdqLP, ratio}; %#ok<AGROW>
            end
        end
    end
    fclose(fid);
end

rows_sw  = parseTable(file_sweep);
rows_ntw = parseTable(file_ntw);

% --- torque row index (zero values in NTW18 = rated torque constant)
torque_idx = 3;

% --- Helper: format one cell as "a / b = r" or "0 / 0 = ---"
function s = fmtCell(pi_val, mpc_val, ratio_val)
    if pi_val == 0 && mpc_val == 0
        s = '$0 / 0 = {-}$';
    else
        s = sprintf('$%.2f / %.2f = %.2f$', pi_val, mpc_val, ratio_val);
    end
end

fid = fopen(file_out, 'w');
fprintf(fid, '\\begin{table}[h]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption[Controller comparison: PI baseline vs.\\ qLMPC for both test cases.]%%\n');
fprintf(fid, '{Controller comparison: PI baseline vs.\\ qLMPC for stepped wind sweep and\n');
fprintf(fid, 'turbulent wind (18\\,m/s mean, IEC class B).\n');
fprintf(fid, 'Each cell shows $\\sigma_\\text{PI} / \\sigma_\\text{qLMPC} = \\text{ratio}$;\n');
fprintf(fid, 'values above~1 indicate improvement.}\n');
fprintf(fid, '\\label{tab:ControllerComparisonMerged}\n');
fprintf(fid, '\\begin{tabular}{ll cc}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, '{Signals} & {Unit} & {Stepped wind sweep} & {Turbulent wind} \\\\\n');
fprintf(fid, ' & & \\multicolumn{2}{c}{$\\sigma_\\text{PI} / \\sigma_\\text{qLMPC} = \\text{ratio}$} \\\\\n');
fprintf(fid, '\\midrule\n');

for i = 1:numel(rows_sw)
    sw  = rows_sw{i};
    ntw = rows_ntw{i};

    cell_sw = fmtCell(sw{3},  sw{4},  sw{5});

    if i == torque_idx
        cell_ntw = '$0 / 0 = {-}$';
    else
        cell_ntw = fmtCell(ntw{3}, ntw{4}, ntw{5});
    end

    fprintf(fid, '%s & %s & %s & %s \\\\\n', sw{1}, sw{2}, cell_sw, cell_ntw);
end

fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
fprintf('Merged table written to %s\n', file_out);
end