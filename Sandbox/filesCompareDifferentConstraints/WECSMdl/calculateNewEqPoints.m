
tmp = load('NREL5MW_linearised.mat');
Lin_points0 = tmp.Lin;


for idx = 1: length(Lin_points0.V)
    aV = Lin_points0.V(idx);
    inTimeWind = [OutTable.Time aV*ones(size(OutTable.Time))];
    sim('testOL_SimulinkMdl2');
    OutDataTable = array2table(OutDataTest, 'VariableNames', {'Wind','RotSpd', 'GenPwr', 'GenTq', 'Pitch', 'Twr_{FA}', 'Twr_{SW}'});
    Lin_points1.RSpeed(idx) = OutDataTable.RotSpd(end)/(60/(2*pi));
    Lin_points1.Pitch(idx) = OutDataTable.Pitch(end)/(180/(pi));
    Lin_points1.Torque(idx) = OutDataTable.GenTq(end)/1000;
end

Lin_points0.RSpeed = Lin_points1.RSpeed;
Lin_points0.Pitch = Lin_points1.Pitch;
Lin_points0.Torque = Lin_points1.Torque;