function getDDEntries(dictName,outMatDir)
%getDDEntries gets variables from Simulink data dictionaries and assigns
%them to the base workspace.
%
% Inputs:
% dict_name - dictionay name (optional)
% outMatDir - output directory for mat files (optional)
currentDir = fileparts(mfilename('fullpath'));
parentDir = fileparts(currentDir);

if nargin <1 || isempty(dictName)
    dictName = 'NREL5WTUDelft';
end

if nargin <2 || isempty(outMatDir)
    outMatDir = fullfile(parentDir,'outMat');
end

%% Load entries from DD
hDict = Simulink.data.dictionary.open(fullfile(parentDir,[dictName,'.sldd']));
hDesignData = hDict.getSection('Global');
childNamesList = hDesignData.evalin('who');
for n = 1:numel(childNamesList)
    hEntry = hDesignData.getEntry(childNamesList{n});
    assignin('base', hEntry.Name, hEntry.getValue);
end

%% Overwrite sampling ime with entry from OutList
if evalin('base','~exist(''DT'')')
    tmp = load(fullfile(outMatDir,'myOutList.mat'));
    assignin('base','DT',tmp.DT);
else
    evalin('base','Control.DT = DT;');
end