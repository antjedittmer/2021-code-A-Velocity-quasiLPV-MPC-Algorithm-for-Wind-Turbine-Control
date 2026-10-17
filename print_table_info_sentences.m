% =========================================================================
% Dynamic Sentence Generator from LaTeX Table (`results_tableMerged.tex`)
% =========================================================================

filename = 'results_tableMerged.tex';

% Read the contents of the LaTeX file
if ~isfile(filename)
    error('File %s not found in the current directory.', filename);
end
fileText = fileread(filename);
lines = splitlines(fileText);

% Helper function to extract the numerical ratio after the '=' sign
extractRatio = @(str) str2double(char(regexp(str, '=\s*([0-9.]+)', 'tokens', 'once')));

% Initialize variables
stepped = struct();
turbulent = struct();

% Loop through lines to find the relevant rows
for i = 1:numel(lines)
    if contains(lines{i}, 'Tower fore-aft acc.')
        parts = split(lines{i}, '&');
        stepped.tower_ratio = extractRatio(parts{3});
        turbulent.tower_ratio = extractRatio(parts{4});
        
    elseif contains(lines{i}, 'Power tracking error')
        parts = split(lines{i}, '&');
        stepped.power_ratio = extractRatio(parts{3});
        turbulent.power_ratio = extractRatio(parts{4});
    end
end

% =========================================================================
% Calculations (Variance reduction factor = ratio^2)
% =========================================================================
stepped.tower_var = stepped.tower_ratio^2;
stepped.power_var = stepped.power_ratio^2;

turbulent.tower_var = turbulent.tower_ratio^2;
turbulent.power_var = turbulent.power_ratio^2;

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
    'in the stepped sweep case, particularly in above-rated conditions. ' ...
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