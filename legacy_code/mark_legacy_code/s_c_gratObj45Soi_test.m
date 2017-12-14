% s_c_gratObj45Soi_test.m

% fillVector = 0.4:0.01:1.0;
% ratioVector = 0.1:0.01:3.0;
% periodVector = 0.4:0.04:1.0;
% offsetVector = 0.0:0.04:0.8; 

fillVector = 0.5;
ratioVector = 0.25;
periodVector = 0.4;
offsetVector = 0.25; 

optimal_angle = 20.0; 

lambda0 = 1.28; 
d = 0.02;

inputs = {'lambda0',lambda0,'d',d,'fillVector',fillVector,'ratioVector',ratioVector,'periodVector',periodVector,'offsetVector',offsetVector,'optimal_angle',optimal_angle};

gratObj = c_gratObj45Soi(inputs{:});
gratObj = gratObj.setParams;
gratObj = gratObj.runParameterSweep;

save('sim/gratObjFineSweep_2','gratObj')