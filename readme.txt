This folder contains the files to reproduce the results and figures presented in

A. Dittmer, B. Sharan and H. Werner, "A Velocity quasiLPV-MPC Algorithm for 
Wind Turbine Control", 2021 European Control Conference (ECC) 

The figures can be reproduced using the file 'runGeneratePicsPaper.m'. 
The data for Figure 1 is generated using 'compareLinearModels.m' in folder 
'LinMdl'. Figure 2 is plotted with 'runCompareModels.m' and for Figures 3 
and 4 using 'runCompareCtrl.m'. Both files are in the folder 'NonLinMdl'.

In order to run the code, matfiles contained in the folder dataIn are used. 
They were generated with
- FASTTool (https://github.com/TUDelft-DataDrivenControl/FASTTool)

The original code in Tag V1.0 has been tested using:
- Matlab R2019b 
The current code (January 2025) has been tested with Matlab 2021b.

It requires MATLAB, Simulink, and Control System Toolbox.

The Simulink project file 'Ecc21.prj' is provide as well as the necessary 
'resources' folder.

The code was updated in January 2025 with the following changes:
- The side-to-side force is related to the torque acting at hub height, 
  the difference between rotor torque and generator torque: 
  F_xt = ksw * (Tr - Ng * Tg), with rotor torque Tr, generator torque Tg, 
  and gear ratio Ng. The lever is ksw = 2/3 Ht, with tower height Ht, the 
  effective point of action for a triangularly distributed load 
- The aerodynamic center rb was set to rb = 0.75R, with rotor radius R
- A scaling factor 0.94 for the maximum torque coefficient Cq look-up-table was 
  introduced. This scaling factor improved the model fit to the FAST model 
  around rated speed.
- The Bode and timeseries plot now also depict the side-to-side motion

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
FITNESS FOR A PARTICULAR PURPOSE.
