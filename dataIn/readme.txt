In order to run the code, matfiles contained in folder dataIn are used. 
They were generated with
- FASTTool (https://github.com/TUDelft-DataDrivenControl/FASTTool)

The mat files are:
- Lin_points_PICL.mat: Lin_points: struct with fields for equilibrium point 
   values for V, Torque, Pitch, GSpeed, RSpeed, x_op, y_op, u_op. Note that 
   the points are not directly taking from FAST, but from the closed loop 
   FASTTool simulation.
- NREL5MW_CPdata.mat, FASTToolCPnew.mat, NREL5MW_CPdata_BS.mat: contain five 
  numerical matrices
   - Rotor_Lamda: tip speed ratio vector 
   - Rotor_Pitch: blade pitch vector
   - Rotor_cQ: aerodynamic torque look-up table matrix 
   - Rotor_cT: aerodynamic force look-up table matrix 
   - Rotor_cP: aerodynamic power look-up table matrix (not used in simulation) 
 Currently FASTToolCPnew.mat obtained with the FASTTool is called in 
 initModel5MWNREL       
- NREL5MW_linearised_4to25.mat: FASTTool output of linearization sweep from
   4 to 25 m/s.
- OutDataSweep.mat: FASTTool simulation with wind sweep from 4 to 25 m/s, 
  step length 20 s
- OutDataWind18NTW.mat: FASTTool simulation with normal turbulent wind with 
   18 m/s mean   
- OutDataStep.mat: FAST v8.16 simulation (S-function and FAST_Library of 
  this repository, AeroDyn v15.03, PI baseline test_FAST_Baseline, same 
  FAST input files as OutDataSweep.mat) with a uniform stepped wind from 3 
  to 25 m/s, step length 100 s (steps at t = k*100 s, 4 m/s from 100 s), 
  2300 s. This is used to obtain references for rotor speed, pitch and 
  torque at these wind speeds (Lin_points, qLPVMPC_RefLUT, updateDDMdl1). 
  Replaces the earlier FASTTool simulation, which started from standstill 
  and had 3-5 % more aerodynamic power below rated than FAST v8.16 of this 
  repository.
