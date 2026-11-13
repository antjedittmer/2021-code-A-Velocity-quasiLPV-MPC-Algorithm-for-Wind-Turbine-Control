function FASTparam = setFASTparam(FAST_InputFileName,FASTparam,dirInfo)
% FASTparam updates the FAST *.dat files. Content of the dat files is read
% into MATLAB with FAST2Matlab, modified and written back. The changes are:
% * ServoDyn: CONTROL via Simulink (4) or DLL (PITCH and Generator/torque)
% * ServoDyn: Ctrl DLL file used
% * ElastoDyn: Add additional Nacelle accelearations 

% ServoDyn
FASTparam.FastInfo = FAST2Matlab(FAST_InputFileName); %FastInfoTable = struct2table(FASTparam.FastInfo); %AD ToDo does this need to be R2007 compatible
FASTparam.idxTMax = ismember(FASTparam.FastInfo.Label,'TMax');
FASTparam.FastInfo.Val{FASTparam.idxTMax} = FASTparam.TMax;

FASTparam.ServoFile = strrep(FASTparam.FastInfo.Val{ismember(FASTparam.FastInfo.Label,'ServoFile')},'"','');
FASTparam.ServoFileName = fullfile(dirInfo.CertTest_directory, FASTparam.ServoFile);
FASTparam.ServoInfo = FAST2Matlab(FASTparam.ServoFileName); %table(FastServoInfo.Label,FastServoInfo.Val)

% ServoDyn: PITCH CONTROL 
FASTparam.idxPContrl = contains(FASTparam.ServoInfo.Label,'PCMode'); % FASTparam.PitchControlMode = FastServoInfo.Val{contains(FastServoInfo.Label,'PCMode')};
FASTparam.ServoInfo.Val{FASTparam.idxPContrl} = FASTparam.PitchControlMode;

% ServoDyn: GENERATOR AND TORQUE CONTROL
FASTparam.idxVSContrl = contains(FASTparam.ServoInfo.Label,'VSContrl'); % FASTparam.GeneratorTorqueControlMode = FastServoInfo.Val{contains(idxVSContrl)}; % 'VSContrl'
FASTparam.ServoInfo.Val{FASTparam.idxVSContrl} = FASTparam.GeneratorTorqueControlMode;

% ServoDyn: Ctrl DLL file
FASTparam.idxDLL = contains(FASTparam.ServoInfo.Label,'DLL_FileName'); %   '"ServoData/DISCON_x64.dll"' 
FASTparam.origDLLPath =  '"ServoData/DISCON_x64.dll"' ;
FASTparam.newDLLPath = strrep(FASTparam.origDLLPath,'DISCON_x64',FASTparam.DLLServoName);
FASTparam.ServoInfo.Val{FASTparam.idxDLL} = FASTparam.newDLLPath;

% ServoDyn: Ctrl DLL in file
FASTparam.idxDLLin = contains(FASTparam.ServoInfo.Label,'DLL_InFile');
FASTparam.origDLLinPath =  '"DISCON.IN"' ;
FASTparam.newDLLin = strrep(FASTparam.origDLLinPath,'DISCON.IN',FASTparam.DLLin);
FASTparam.ServoInfo.Val{FASTparam.idxDLLin} = FASTparam.newDLLin;

%ServoDynFile parameter to add
listServoParam = {'"BlPitchC1"','"BlPitchC2"','"BlPitchC3"'};
listServoComments = {...
    '- Blade 1 pitch angle command',...
    '- Blade 2 pitch angle command',...
    '- Blade 3 pitch angle command'};

for idx = 1: length(listServoParam)
    if ~ismember(listServoParam{idx},FASTparam.ServoInfo.OutList)
        FASTparam.ServoInfo.OutList{end+1} = listServoParam{idx} ;
        FASTparam.ServoInfo.OutListComments{end+1} = listServoComments{idx};
    end
end

% ElastoDynFile
FASTparam.ElastoDynFile = strrep(FASTparam.FastInfo.Val{ismember(FASTparam.FastInfo.Label,'EDFile')},'"','');
FASTparam.ElastoDynFileName = fullfile(dirInfo.CertTest_directory, FASTparam.ElastoDynFile);
FASTparam.ElastoDynInfo = FAST2Matlab(FASTparam.ElastoDynFileName); %table(FastServoInfo.Label,FastServoInfo.Val)

% ElastoDynFile parameter to add
listEDparam = {'"NcIMUTAxs"','"NcIMUTAys"','"RootMyb1"','"RootMyb2"', '"RootMyb3"'};
listEDcomments = {...
    '- Nacelle inertial measurement unit translational acceleration (absolute)',...
    '- Nacelle inertial measurement unit translational acceleration (absolute)',...
    '- In-plane bending, out-of-plane bending, and pitching moments at the root of blade 1',...
    '- In-plane bending, out-of-plane bending, and pitching moments at the root of blade 2',...
    '- In-plane bending, out-of-plane bending, and pitching moments at the root of blade 3'};

for idx = 1: length(listEDparam)
    if ~ismember(listEDparam{idx},FASTparam.ElastoDynInfo.OutList)
        FASTparam.ElastoDynInfo.OutList{end+1} = listEDparam{idx} ;
        FASTparam.ElastoDynInfo.OutListComments{end+1} = listEDcomments{idx};
    end
end


% FASTparam.OutList = [FASTparam.InflowInfo.OutList; FASTparam.ElastoDynInfo.OutList]



Matlab2FAST(FASTparam.FastInfo,FAST_InputFileName,'tmp.dat',2)
Matlab2FAST(FASTparam.FastInfo,'tmp.dat',FAST_InputFileName,2)
delete('tmp.dat')

% Write FAST Servo.dat file with new parameter
Matlab2FAST(FASTparam.ServoInfo,FASTparam.ServoFileName,'tmp.dat',2)
Matlab2FAST(FASTparam.ServoInfo,'tmp.dat',FASTparam.ServoFileName,2)
delete('tmp.dat')

% Write FAST ElastoDynInfo.dat file with new parameter
Matlab2FAST(FASTparam.ElastoDynInfo,FASTparam.ElastoDynFileName,'tmp.dat',2)
Matlab2FAST(FASTparam.ElastoDynInfo,'tmp.dat',FASTparam.ElastoDynFileName,2)
delete('tmp.dat')

% Write FAST Inflow file with new parameter
% Wind Info
if ~ isempty(FASTparam.InflowCell)
    FASTparam.InflowFile = strrep(FASTparam.FastInfo.Val{ismember(FASTparam.FastInfo.Label,'InflowFile')},'"','');
    FASTparam.InflowFileName = fullfile(dirInfo.CertTest_directory, FASTparam.InflowFile);
    FASTparam.InflowInfo = FAST2Matlab(FASTparam.InflowFileName); %table(FASTparam.InflowInfo.Label,FASTparam.InflowInfo.Val)
    
    FASTparam.idxInflowFilename = contains(FASTparam.InflowInfo.Label,'Filename');
    FASTparam.idxInflowWindType = contains(FASTparam.InflowInfo.Label,'WindType');
    FASTparam.InflowInfo.Val(FASTparam.idxInflowFilename) = FASTparam.InflowCell;
    FASTparam.InflowInfo.Val{FASTparam.idxInflowWindType} = FASTparam.WindType;
    Matlab2FAST(FASTparam.InflowInfo,FASTparam.InflowFileName,'tmp.dat',2)
    Matlab2FAST(FASTparam.InflowInfo,'tmp.dat',FASTparam.InflowFileName,2)
    delete('tmp.dat')
end
end
