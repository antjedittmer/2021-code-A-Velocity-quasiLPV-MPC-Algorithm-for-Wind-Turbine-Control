function cleanupSimulinkForThesis(modelName)
% cleanupSimulinkForThesis - Prepare Simulink model for thesis inclusion
% Usage: cleanupSimulinkForThesis('your_model_name')
%
% This function:
%   - Renames variables to match LaTeX notation
%   - Removes colored backgrounds from subsystems
%   - Increases font size for readability
%   - Removes default Mux annotations
%   - Sets consistent block colors

if nargin < 1
    modelName = bdroot;  % use currently open model
end

load_system(modelName);

blocks = find_system(modelName);

%% 1. Set global font properties for readability
set_param(modelName, 'DefaultBlockFontSize', 14);
set_param(modelName, 'DefaultBlockFontName', 'Arial');
set_param(modelName, 'DefaultLineFontSize', 12);
set_param(modelName, 'DefaultAnnotationFontSize', 12);

%% 2. Rename signal labels to LaTeX-friendly notation
% Map of Simulink names -> thesis notation (using Unicode where possible,
% or LaTeX-style names that you can post-process)
 renameMap = {
    'xdotfa',           'xt_dot';
    'zetadot',          'delta_b_dot';
    'omega_r',          'omega_r';
    'x_t_fa',           'x_t';
    'zeta_b',           'delta_b';
    'y_t_sw',           'y_t';
    'omega_g',          'omega_g';
    'Tau_aero',         'T_r';
    'F_T',              'F_T_r';
    'wind [m/s]',       'V_inf [m/s]';
    'beta [rad]',       'beta_0 [rad]';
    'Tg [Nm]',          'T_g [Nm]';
    'omega_r_init [rad/s]', 'omega_r,0 [rad/s]';
};

% Rename Goto/From tags
gotos = find_system(modelName, 'BlockType', 'Goto');
froms = find_system(modelName, 'BlockType', 'From');
allTagBlocks = [gotos; froms];

for k = 1:size(renameMap, 1)
    oldTag = renameMap{k, 1};
    newTag = renameMap{k, 2};
    for b = 1:length(allTagBlocks)
        if strcmp(get_param(allTagBlocks{b}, 'GotoTag'), oldTag)
            set_param(allTagBlocks{b}, 'GotoTag', newTag);
        end
    end
end

% Rename Inport/Outport labels
ports = [find_system(modelName, 'BlockType', 'Inport'); 
         find_system(modelName, 'BlockType', 'Outport')];
for k = 1:size(renameMap, 1)
    for p = 1:length(ports)
        try
        currentName = get_param(ports{p}, 'Name');
        if strcmp(currentName, renameMap{k, 1})
            set_param(ports{p}, 'Name', renameMap{k, 2});
        end
        catch
        end
    end
end

%% 3.5 Normalize Goto/From appearance (size + font consistency)

gotoBlocks = find_system(modelName, 'BlockType', 'Goto');
fromBlocks = find_system(modelName, 'BlockType', 'From');
allGF = [gotoBlocks; fromBlocks];

% Match Inport font size (fallback to 12 if not found)
portFontSize = 14;
ports = find_system(modelName, 'BlockType', 'Inport');
if ~isempty(ports)
    try
        portFontSize = str2double(get_param(ports{1}, 'FontSize'));
    catch
        portFontSize = 14;
    end
end

% Tags associated with F_Tr or T_r
highlightTags = {'F_Tr', 'T_r'};

for i = 1:length(allGF)
    blk = allGF{i};

    % Get tag (GotoTag for Goto, Name for From display consistency)
    try
        tag = get_param(blk, 'GotoTag');
    catch
        tag = get_param(blk, 'Name');
    end

    % Estimate width based on text length (tune factor if needed)
    baseCharWidth = 6;  % pixels per character (empirical)
    minWidth = 30;
    padding = 20;
    w = max(minWidth, length(tag) * baseCharWidth + padding);
    h = 18; % fixed height for consistency
    pos = get_param(blk, 'Position');

    % Keep center fixed, adjust width/height
    cx = (pos(1) + pos(3)) / 2;
    cy = (pos(2) + pos(4)) / 2;
    newPos = [cx - w/2, cy - h/2, cx + w/2, cy + h/2];
    set_param(blk, 'Position', round(newPos));

    % Set font size to match Inport/Outport
    try
        set_param(blk, 'FontSize', portFontSize);
    catch
        % some versions may ignore this for From/Goto
    end

    % Set colors based on tag
    if ~any(strcmp(tag, highlightTags))
        try
            set_param(blk, 'ForegroundColor', 'black');
            set_param(blk, 'BackgroundColor', 'white');
        catch
            % ignore if color setting not supported
        end
    end
end
%% 3. Remove colored backgrounds from subsystems
 
%% 4. Set neutral colors for all blocks
% allBlocks = find_system(modelName, 'Type', 'Block');
% for b = 1:length(allBlocks)
%     try
%         % Only change blocks that have unusual colors
%         bgColor = get_param(allBlocks{b}, 'BackgroundColor');
%         if ~ismember(bgColor, {'white', 'black', 'gray'})
%             set_param(allBlocks{b}, 'BackgroundColor', 'white');
%         end
%         fgColor = get_param(allBlocks{b}, 'ForegroundColor');
%         if ~strcmp(fgColor, 'black')
%             set_param(allBlocks{b}, 'ForegroundColor', 'black');
%         end
%     catch
%         % skip blocks that don't support these properties
%     end
% end

%% 5. Hide block names for cleanly-labeled blocks (optional)
% Uncomment if you want to hide default block names like "Mux5"
% muxBlocks = find_system(modelName, 'BlockType', 'Mux');
% for m = 1:length(muxBlocks)
%     set_param(muxBlocks{m}, 'ShowName', 'off');
% end

%% 6. Update model display
set_param(modelName, 'ShowPageBoundaries', 'off');

fprintf('Cleanup complete for model: %s\n', modelName);
fprintf('Remember to: \n');
fprintf('  1. Manually reposition any overlapping blocks\n');
fprintf('  2. Export as vector format using exportThesisFigure() below\n');

end


function exportThesisFigure(modelName, outputFile)
% Export Simulink model as vector PDF for thesis
% Usage: exportThesisFigure('your_model', 'figure.pdf')

if nargin < 2
    outputFile = [modelName '.pdf'];
end

% Print to PDF (vector format)
print(['-s' modelName], '-dpdf', '-bestfit', outputFile);

% Alternative: EPS format
% print(['-s' modelName], '-depsc2', [modelName '.eps']);

fprintf('Exported to: %s\n', outputFile);
end