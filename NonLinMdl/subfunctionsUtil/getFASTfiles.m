function [FASTfiles,FASTfilesDesc,OnePlot] = getFASTfiles(dataFile,existsFile,FASTparam, userConfig)
% getFASTfiles sets name string of fast files and legend string
%
% Input:
%  - existsFile: struct with information if bin or text data file exist
%  - userConfig: enables to plot comparison with shipped data files

OnePlot = 3;
if existsFile.Data && existsFile.Comparison && userConfig.runCmp
    FASTfiles = {dataFile.data,dataFile.ComparisonFile};
    FASTfilesDesc = {'SFunc','exe'};
elseif existsFile.Bin && existsFile.Comparison && userConfig.runCmp
    FASTfiles = {dataFile.binFile,dataFile.ComparisonFile};
    FASTfilesDesc = {'SFunc','exe'};
elseif existsFile.Bin && existsFile.Comparison && userConfig.runCmp
    FASTfiles = {dataFile.binFile,dataFile.ComparisonFile};
    FASTfilesDesc = {'SFunc','exe'};
elseif existsFile.Bin && existsFile.BinComparisonFile && userConfig.runCmp
    FASTfiles = {dataFile.binFile,dataFile.binComparisonFile};
    FASTfilesDesc = {['SimulinkCtrl',FASTparam.strTest],['FASTCtrl',FASTparam.strTest]};
elseif existsFile.Bin
    FASTfiles = {dataFile.binFile};
    FASTfilesDesc =  {'SFunc'};
    OnePlot = 2;
elseif existsFile.Data
    FASTfiles = {dataFile.data};
    FASTfilesDesc = {'SFunc'};
    OnePlot = 2;
else
    warning('No new simulation results for plotting. No plots generated')
    return;
end