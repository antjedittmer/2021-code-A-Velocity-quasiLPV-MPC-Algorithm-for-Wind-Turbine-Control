function merge_stats_tables(file_mpc, file_pi, file_out)
% merge_stats_tables('statsTabletestConstr.tex',
%                    'statsTabletestConstrPI.tex',
%                    'statsTableMerged.tex')
if nargin < 1
    fd ='figDir';
    file_mpc = fullfile(fd,'statsTabletestConstr.tex');
    file_pi  = fullfile(fd,'statsTabletestConstrPI.tex');
    file_out = fullfile(fd,'statsTableMerged.tex');
end

% --- Parse std table from a file (second table = stddev) ---
function [cols, data] = parseStdTable(filename)
    fid  = fopen(filename, 'r');
    text = fread(fid, '*char')';
    fclose(fid);

    % Find second \begin{tabular} block (stddev table)
    starts = strfind(text, '\begin{tabular}');
    ends   = strfind(text, '\end{tabular}');
    block  = text(starts(2):ends(2));

    lines = strsplit(block, '\n');
    cols  = [];
    data  = {};

    for i = 1:numel(lines)
        line = strtrim(lines{i});

        % Parse header to get column order
        if contains(line, 'low') && contains(line, 'high') && contains(line, '&')
            parts = strsplit(line, '&');
            cols  = strtrim(parts(2:end));
            cols  = cellfun(@(s) strrep(s,'\\',''), cols, 'UniformOutput', false);
            cols  = strtrim(cols);

        % Parse data rows
        elseif contains(line, '&') && contains(line, '\\') && ~contains(line, 'tabular')
            parts = strsplit(line, '&');
            if numel(parts) == 4
                vals = cellfun(@(s) str2double(strtrim(strrep(s,'\\',''))), parts(2:end));
                data{end+1} = vals; %#ok<AGROW>
            end
        end
    end
end

[cols_mpc, data_mpc] = parseStdTable(file_mpc);
[cols_pi,  data_pi ] = parseStdTable(file_pi);

% --- Reorder columns to: low (4), standard (8), high (13) ---
desired = {'low (4)', 'standard (8)', 'high (13)'};

function idx = colOrder(cols, desired)
    idx = zeros(1, numel(desired));
    for k = 1:numel(desired)
        match = find(strcmpi(strtrim(cols), desired{k}));
        if isempty(match)
            error('Column "%s" not found in: %s', desired{k}, strjoin(cols, ', '));
        end
        idx(k) = match;
    end
end

idx_mpc = colOrder(cols_mpc, desired);
idx_pi  = colOrder(cols_pi,  desired);

% Reorder each row
for i = 1:numel(data_mpc)
    data_mpc{i} = data_mpc{i}(idx_mpc);
end
for i = 1:numel(data_pi)
    data_pi{i} = data_pi{i}(idx_pi);
end

% --- Row display names (in desired output order) ---
rowNames = {
    'Tower fore-aft acc.\ (m/s$^2$)', ...
    'Generated power (kW)', ...
    'Rotor speed (rpm)', ...
    'Blade pitch rate (\textdegree/s)'
};

% --- Baseline: PI at standard rate limit (8 deg/s = column 2) ---
iBase = 2;
baseline = cellfun(@(v) v(iBase), data_pi);

% --- Write merged table ---
fid = fopen(file_out, 'w');
fprintf(fid, '\\begin{table}[htbp]\n');
fprintf(fid, '    \\centering\n');
fprintf(fid, '    \\caption[Standard deviation comparison: Baseline PI vs.\\ qLMPC for turbulent wind with different pitch rate limits.]%%\n');
fprintf(fid, '    {Standard deviation comparison: Baseline PI vs.\\ qLMPC\n');
fprintf(fid, '        for turbulent wind with pitch rate limits. The PI controller is largely insensitive to the rate limit\n');
fprintf(fid, '        whereas the qLMPC uses higher rates for less rotor speed, tower movement and\n');
fprintf(fid, '        power variation. Values in brackets show the ratio to the standard\n');
fprintf(fid, '        (8\\textdegree/s) baseline PI value.}\n');
fprintf(fid, '    \\label{tab:stddev_unified}\n');
fprintf(fid, '    \\begin{tabular}{ll rrr}\n');
fprintf(fid, '        \\toprule\n');
fprintf(fid, '        & & \\multicolumn{3}{c}{Rate limit (\\textdegree/s) (ratio to PI at 8\\textdegree/s)} \\\\\n');
fprintf(fid, '        \\cmidrule(lr){3-5}\n');
fprintf(fid, '        Controller & Signal & low: 4\\textdegree/s & standard: 8\\textdegree/s & high: 13\\textdegree/s \\\\\n');
fprintf(fid, '        \\midrule\n');

writeBlock(fid, 'PI',    rowNames, data_pi,  baseline);
fprintf(fid, '        \\midrule\n');
writeBlock(fid, 'qLMPC', rowNames, data_mpc, baseline);

fprintf(fid, '        \\bottomrule\n');
fprintf(fid, '    \\end{tabular}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
fprintf('Merged table written to %s\n', file_out);

% --- Write one controller block: value (baseline/value) ---
function writeBlock(fid, name, rowNames, data, baseline)
    fprintf(fid, '        \\multirow{%d}{*}{%s}\n', numel(rowNames), name);
    for r = 1:numel(rowNames)
        v = data{r};
        q = baseline(r) ./ v;
        fprintf(fid, '        & %s & %.2f (%.2f) & %.2f (%.2f) & %.2f (%.2f) \\\\\n', ...
            rowNames{r}, v(1), q(1), v(2), q(2), v(3), q(3));
    end
end
end