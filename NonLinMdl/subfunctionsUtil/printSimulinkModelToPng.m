% Screenshot Simulink model and all subsystems to SVG and PDF
% Usage: Open your model, then run this script
modelName = gcs;                 % Current system (or replace with model name)
outputDir = fullfile(pwd,'picsSimulink');                 % Output directory (current folder)


% Load the model if not already loaded
load_system(modelName);
% Get the root model name
rootModel = bdroot(modelName);
% Find all subsystems in the model
subsystems = find_system(rootModel, 'BlockType', 'SubSystem');
% Combine root model with all subsystems
allSystems = [{rootModel}; subsystems];
fprintf('Found %d systems to capture:\n', length(allSystems));
% Loop through and capture each system
for i = 1:length(allSystems)
    currentSystem = allSystems{i};
    
    % Create safe filename by replacing '/' with '_'
    safeName = strrep(currentSystem, '/', '_');
    baseName = fullfile(outputDir, safeName);
    
    % Open the system to make it visible
    open_system(currentSystem);
    
    % Brief pause to ensure rendering is complete
    pause(0.5);
    
    try
        % Save as SVG
        print(['-s', currentSystem], '-dsvg', [baseName, '.svg']);
        
        % Save as PDF
        print(['-s', currentSystem], '-dpdf', [baseName, '.pdf']);
        
        fprintf('  [%d/%d] Saved: %s\n', i, length(allSystems), safeName);
    catch ME
        fprintf('  [%d/%d] FAILED: %s - %s\n', i, length(allSystems), safeName, ME.message);
    end
end