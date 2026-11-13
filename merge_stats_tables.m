function merge_stats_tables(file_mpc, file_pi, file_out, captionCase)
% merge_stats_tables('statsTabletestConstr.tex',
%                    'statsTabletestConstrPI.tex',
%                    'statsTableMerged.tex')
% captionCase (optional): 'piNoFilt' for the tables without wind filter
% (runGeneratePicsPaper: statsTabletestConstr_Filt8.tex and
% statsTabletestConstrPI_NoFilt.tex). Both controllers then have the
% filtered case as last column; the qLMPC value of this column is
% normalized with the filtered PI value, the PI value with the unfiltered
% PI at 8 deg/s.
if nargin < 4 || isempty(captionCase)
    captionCase = '';
end
isPINoFilt = strcmpi(captionCase, 'piNoFilt');
if nargin < 1
    fd ='figDir';
    file_mpc = fullfile(fd,'statsTabletestConstr.tex');
    file_pi  = fullfile(fd,'statsTabletestConstrPI.tex');
    file_out = fullfile(fd,'statsTableMerged.tex');
end

% --- Parse std table from a file (second table = stddev) ---
% The unrounded values of the mat file of plotOutTable (same name as the
% tex file) are used if available, so that the ratios are computed from
% unrounded values; otherwise the rounded values of the tex file
function [cols, data] = parseStdTable(filename)
    [fPath, fName] = fileparts(filename);
    matFile = fullfile(fPath, [fName, '.mat']);
    if isfile(matFile)
        S = load(matFile, 'stdDevTable');
        cols = strtrim(S.stdDevTable.Properties.VariableNames(2:end));
        data = arrayfun(@(r) S.stdDevTable{r, 2:end}, 1:height(S.stdDevTable), 'UniformOutput', false);
        return
    end
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
            if numel(parts) >= 4 % three rate limits, optionally the fourth case
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

% Optional fourth case (runGeneratePicsPaper, addFilteredCase), e.g.
% 'filt. (13)' or 'filt. (8)': controller on the filtered plant
iFiltMPC = find(startsWith(strtrim(cols_mpc), 'filt.'), 1);
hasFilt = ~isempty(iFiltMPC);
hasFiltPI = false;
rateFilt = '';
if hasFilt
    nameFilt = strtrim(cols_mpc{iFiltMPC});
    rateFilt = regexp(nameFilt, '\d+', 'match', 'once');
    idx_mpc = [idx_mpc, iFiltMPC];
    % The PI table can have the same fourth case
    iFiltPI = find(strcmpi(strtrim(cols_pi), nameFilt), 1);
    hasFiltPI = ~isempty(iFiltPI);
    if hasFiltPI
        idx_pi = [idx_pi, iFiltPI];
    end
end

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
    'Power tracking error (kW)', ...
    'Rotor speed (rpm)', ...
    'Blade pitch rate (\textdegree/s)'
};

% --- Baseline: PI at standard rate limit (8 deg/s = column 2) ---
nCol = 3 + hasFilt;
iBase = 2;
baseline = cellfun(@(v) v(iBase), data_pi);
baselinePI  = repmat(baseline(:), 1, nCol);
baselineMPC = baselinePI;
if isPINoFilt && hasFiltPI
    % Filtered column of the qLMPC: normalized with the filtered PI (same plant)
    baselineMPC(:, 4) = cellfun(@(v) v(4), data_pi);
end

% --- Write merged table ---
fid = fopen(file_out, 'w');
fprintf(fid, '\\begin{table}[htbp]\n');
fprintf(fid, '    \\centering\n');
if isPINoFilt
    fprintf(fid, '    \\caption[Standard deviation comparison: Baseline PI vs.\\ qLMPC without wind filter for turbulent wind with different pitch rate limits.]%%\n');
    fprintf(fid, '    {Standard deviation comparison: Baseline PI vs.\\ qLMPC\n');
    fprintf(fid, '        for turbulent wind with pitch rate limits. Columns 4, 8 and 13\\textdegree/s: both controllers on the\n');
    fprintf(fid, '        same plant without the rotor-effective wind filter. Values in brackets show the ratio to the\n');
    fprintf(fid, '        standard (8\\textdegree/s) baseline PI value without filter. Ratios are computed from unrounded values.\n');
    if hasFilt
        fprintf(fid, '        Last column: both controllers at %s\\textdegree/s on the plant with the filter, the qLMPC scheduled\n', rateFilt);
        fprintf(fid, '        with the filtered wind. The qLMPC value is normalized with the filtered PI value, the PI value with\n');
        fprintf(fid, '        the PI at 8\\textdegree/s without filter (the only ratio of two different plants).}\n');
    else
        fprintf(fid, '        }\n');
    end
    fprintf(fid, '    \\label{tab:stddev_unified_NoFilt}\n');
else
    fprintf(fid, '    \\caption[Standard deviation comparison: Baseline PI vs.\\ qLMPC for turbulent wind with different pitch rate limits.]%%\n');
    fprintf(fid, '    {Standard deviation comparison: Baseline PI vs.\\ qLMPC\n');
    fprintf(fid, '        for turbulent wind with pitch rate limits. The PI controller is largely insensitive to the rate limit\n');
    fprintf(fid, '        whereas the qLMPC uses higher rates for less rotor speed, tower movement and\n');
    fprintf(fid, '        power variation. Values in brackets show the ratio to the standard\n');
    if hasFilt
        fprintf(fid, '        (8\\textdegree/s) baseline PI value. Last column: qLMPC at %s\\textdegree/s on the plant with the\n', rateFilt);
        fprintf(fid, '        rotor-effective wind filter and scheduled with the filtered wind.}\n');
    else
        fprintf(fid, '        (8\\textdegree/s) baseline PI value.}\n');
    end
    fprintf(fid, '    \\label{tab:stddev_unified}\n');
end
fprintf(fid, '    \\begin{tabular}{ll %s}\n', repmat('r', 1, nCol));
fprintf(fid, '        \\toprule\n');
fprintf(fid, '        & & \\multicolumn{%d}{c}{Rate limit (\\textdegree/s) (ratio to PI at 8\\textdegree/s)} \\\\\n', nCol);
fprintf(fid, '        \\cmidrule(lr){3-%d}\n', 2 + nCol);
if hasFilt
    fprintf(fid, '        Controller & Signal & low: 4\\textdegree/s & standard: 8\\textdegree/s & high: 13\\textdegree/s & filtered: %s\\textdegree/s \\\\\n', rateFilt);
else
    fprintf(fid, '        Controller & Signal & low: 4\\textdegree/s & standard: 8\\textdegree/s & high: 13\\textdegree/s \\\\\n');
end
fprintf(fid, '        \\midrule\n');

writeBlock(fid, 'PI',    rowNames, data_pi,  baselinePI,  nCol);
fprintf(fid, '        \\midrule\n');
writeBlock(fid, 'qLMPC', rowNames, data_mpc, baselineMPC, nCol);

fprintf(fid, '        \\bottomrule\n');
fprintf(fid, '    \\end{tabular}\n');
fprintf(fid, '\\end{table}\n');
fclose(fid);
fprintf('Merged table written to %s\n', file_out);

% --- Write one controller block: value (baseline/value) ---
% (nCol columns, baseline: rows x nCol; a block with fewer values, e.g. the
% PI without the fourth case, gets '--' in the remaining columns)
function writeBlock(fid, name, rowNames, data, baseline, nCol)
    % Number formats of the rows: tower acceleration and rotor speed with
    % four digits, power and pitch rate with two (as plotOutTable)
    rowFormats = {'%.4f', '%.2f', '%.4f', '%.2f'};
    fprintf(fid, '        \\multirow{%d}{*}{%s}\n', numel(rowNames), name);
    for r = 1:numel(rowNames)
        v = data{r};
        fprintf(fid, '        & %s', rowNames{r});
        for c = 1:nCol
            if c <= numel(v)
                fprintf(fid, [' & ', rowFormats{r}, ' (%.2f)'], v(c), baseline(r,c) / v(c));
            else
                fprintf(fid, ' & --');
            end
        end
        fprintf(fid, ' \\\\\n');
    end
end
end
