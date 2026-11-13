function merge_stats_tables_norm(file_mpc, file_pi, file_out)
if nargin < 1
    fd ='figDir';
    file_mpc = fullfile(fd,'statsTabletestConstr.tex');
    file_pi  = fullfile(fd,'statsTabletestConstrPI.tex');
    file_out = fullfile(fd,'statsTableMerged.tex');
end

function [cols, data] = parseStdTable(filename)
    fid  = fopen(filename, 'r');
    text = fread(fid, '*char')';
    fclose(fid);
    starts = strfind(text, '\begin{tabular}');
    ends   = strfind(text, '\end{tabular}');
    block  = text(starts(2):ends(2));
    lines  = strsplit(block, '\n');
    cols   = [];
    data   = {};
    for i = 1:numel(lines)
        line = strtrim(lines{i});
        if contains(line, 'low') && contains(line, 'high') && contains(line, '&')
            parts = strsplit(line, '&');
            cols  = strtrim(parts(2:end));
            cols  = cellfun(@(s) strrep(s,'\\',''), cols, 'UniformOutput', false);
            cols  = strtrim(cols);
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

for i = 1:numel(data_mpc)
    data_mpc{i} = data_mpc{i}(idx_mpc);
end
for i = 1:numel(data_pi)
    data_pi{i} = data_pi{i}(idx_pi);
end

% Normalisation factors: standard (8) PI values (column index 2)
norm_factors = cellfun(@(r) r(2), data_pi);

rowNames = {
    'Tower fore-aft acc.\ (m/s$^2$)', ...
    'Generated power (kW)', ...
    'Rotor speed (rpm)', ...
    'Blade pitch rate (\textdegree/s)'
};

% Rows to normalise (exclude pitch rate)
normalise_rows = [1, 2, 3,4];

% Helper: format one cell as "val (norm)" or just "val" if not normalised
function s = fmtCell(val, nf, doNorm)
    if doNorm
        s = sprintf('%.2f (%.2f)', val, nf/val);  % flipped: PI/val
    else
        s = sprintf('%.2f', val);
    end
end

fid = fopen(file_out, 'w');
fprintf(fid, '\\begin{table}[htbp]\n');
fprintf(fid, '\\centering\n');
fprintf(fid, '\\caption[Standard deviation comparison of different pitch rate limits: PI vs.\\ qLMPC for turbulent wind.]%%\n');
fprintf(fid, '{Standard deviation comparison of different pitch rate limits: PI vs.\\ qLMPC\n');
fprintf(fid, 'for turbulent wind. The PI controller is largely insensitive to the rate limit\n');
fprintf(fid, 'whereas the qLMPC uses higher rates for less rotor speed, tower movement and\n');
fprintf(fid, 'power variation. Values in brackets show the ratio to the standard\n');
fprintf(fid, '(8\\,\\textdegree/s) PI baseline value.}\n');
fprintf(fid, '\\label{tab:stddev_unified}\n');
fprintf(fid, '\\begin{tabular}{ll rrr}\n');
fprintf(fid, '\\toprule\n');
fprintf(fid, ' & & \\multicolumn{3}{c}{Rate limit (\\textdegree/s) (ratio to PI at 8\\,\\textdegree/s)} \\\\\n');
fprintf(fid, '\\cmidrule(lr){3-5}\n');
fprintf(fid, 'Controller & Signal & low (4) & standard (8) & high (13) \\\\\n');
fprintf(fid, '\\midrule\n');

% Write rows for one controller
function writeRows(fid, data, norm_factors, rowNames, normalise_rows)
    for r = 1:numel(rowNames)
        doNorm = ismember(r, normalise_rows);
        nf     = norm_factors(r);
        c1 = fmtCell(data{r}(1), nf, doNorm);
        c2 = fmtCell(data{r}(2), nf, doNorm);
        c3 = fmtCell(data{r}(3), nf, doNorm);
        fprintf(fid, '    & %s & %s & %s & %s \\\\\n', rowNames{r}, c1, c2, c3);
    end
end

fprintf(fid, '\\multirow{4}{*}{PI}\n');
writeRows(fid, data_pi, norm_factors, rowNames, normalise_rows);

fprintf(fid, '\\midrule\n');

fprintf(fid, '\\multirow{4}{*}{qLMPC}\n');
writeRows(fid, data_mpc, norm_factors, rowNames, normalise_rows);

fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
fprintf('Merged table written to %s\n', file_out);
end