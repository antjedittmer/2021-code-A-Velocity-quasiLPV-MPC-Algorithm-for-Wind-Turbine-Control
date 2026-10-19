% =========================================================================
% Dynamic Sentence Generator from LaTeX Table (`results_tableMerged.tex`)
% =========================================================================
% Expected row format (label, stepped sweep, turbulent wind):
% Tower fore-aft acc. (m/s$^2$) & $0.0514 / 0.0319 = 1.61$ & $0.1300 / 0.0765 = 1.70$ \\
% Ratios and variance factors are calculated from the standard deviations
% in the table (not from the rounded ratios) and checked against the table.

filename = 'results_tableMerged.tex';

% Read the contents of the LaTeX file
if ~isfile(filename)
    error('File %s not found in the current directory.', filename);
end
fileText = fileread(filename);
tableLines = splitlines(fileText); % do not use 'lines', it shadows a MATLAB function

% Rows to extract: field name and label in the table
rowDefs = {'tower', 'Tower fore-aft acc.';
           'power', 'Power tracking error'};

% Initialize variables
stepped   = struct();
turbulent = struct();

% Loop through lines to find the relevant rows
for i = 1:numel(tableLines)
    lineStr = strtrim(tableLines{i});
    if isempty(lineStr) || startsWith(lineStr, '%')
        continue % skip empty lines and LaTeX comments
    end
    for r = 1:size(rowDefs, 1)
        if contains(lineStr, rowDefs{r,2})
            parts = split(lineStr, '&');
            if numel(parts) < 3
                error('Row "%s" has fewer than two data columns.', rowDefs{r,2});
            end
            % Last two columns: stepped sweep and turbulent wind
            stepped   = parseCell(stepped,   rowDefs{r,1}, parts{end-1});
            turbulent = parseCell(turbulent, rowDefs{r,1}, parts{end});
        end
    end
end

% Check that all rows were found
for r = 1:size(rowDefs, 1)
    if ~isfield(stepped, [rowDefs{r,1}, '_ratio'])
        error('Row "%s" not found in %s.', rowDefs{r,2}, filename);
    end
end

% =========================================================================
% Generate Sentences
% =========================================================================
sentence1 = sprintf([...
    'Overall, the standard deviation ratios for the tower fore-aft ' ...
    'acceleration and power tracking error are %.2f and %.2f ' ...
    '(corresponding to variance reduction factors of %.2f and %.2f), ' ...
    'as shown in Table~\\ref{tab:ControllerComparisonMerged}.'], ...
    stepped.tower_ratio, stepped.power_ratio, ...
    stepped.tower_var, stepped.power_var);

sentence2 = sprintf([...
    'Figure~\\ref{fig:CtrlCmpMdlSweep}(b) shows results for the turbulent wind. ' ...
    'The improvements over the baseline PI are even more pronounced than ' ...
    'in the stepped sweep case. ' ...
    'Regarding the tower fore-aft acceleration and power tracking error, ' ...
    'the qLMPC yields standard deviation ratios of %.2f and %.2f ' ...
    '(equivalent to variance reduction factors of %.2f and %.2f), ' ...
    'as presented in Table~\\ref{tab:ControllerComparisonMerged}.'], ...
    turbulent.tower_ratio, turbulent.power_ratio, ...
    turbulent.tower_var, turbulent.power_var);

% =========================================================================
% Output Results to Command Window
% =========================================================================
disp('--- Stepped Wind Sweep Sentence ---');
disp(sentence1);
fprintf('\n');
disp('--- Turbulent Wind Sentence ---');
disp(sentence2);

% =========================================================================
% Local function
% =========================================================================
function s = parseCell(s, name, cellStr)
% parseCell reads '$a / b = c$' and stores the ratio a/b and the variance
% factor (a/b)^2 in the struct s.
tok = regexp(cellStr, '([0-9.]+)\s*/\s*([0-9.]+)\s*=\s*([0-9.]+)', 'tokens', 'once');
if isempty(tok)
    error('Could not parse "%s" (expected $a / b = c$).', strtrim(cellStr));
end
stdPI    = str2double(tok{1});
stdMPC   = str2double(tok{2});
ratioTab = str2double(tok{3});
ratio    = stdPI / stdMPC;
if abs(ratio - ratioTab) > 0.01
    warning('%s: ratio from standard deviations (%.3f) differs from table ratio (%.2f).', ...
        name, ratio, ratioTab);
end
s.([name, '_ratio']) = ratio;
s.([name, '_var'])   = ratio^2;
end
