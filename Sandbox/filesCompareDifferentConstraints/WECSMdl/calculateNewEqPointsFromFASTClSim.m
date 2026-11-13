
load('OutDataStep.mat');
idxT = OutTable.Time >= 100;
OutTableT = OutTable(idxT,:);

OutTableT.Time =  OutTableT.Time - OutTableT.Time(1);
time = OutTableT.Time;

tmp.Lin.V = 4:25;

for idx = 1 : length(tmp.Lin.V)
    timeIdx = time <= idx*100-10 & time >= idx*100-30;
    V(idx) = mean(OutTableT.Wind1VelX(timeIdx));

    RSpeed(idx) = mean(OutTableT.GenSpeed(timeIdx)/97/(60/(2*pi)));
    Pitch(idx) = mean(OutTableT.BlPitch1(timeIdx)/(180/(pi)));
    Torque(idx) = mean(OutTableT.GenTq(timeIdx));

end

Lin.RSpeed = RSpeed;
Lin.Pitch = Pitch;
Lin.Torque = Torque;