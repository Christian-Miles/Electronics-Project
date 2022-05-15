VDC = 24; 
V1=15.7656;
V3=0.5912; 
V5 = 2.1877;

 f = @(x) [(4/(pi))* VDC *(sin(x(1))-sin(x(2))+sin(x(3)))-V1; 
           (4/(3*pi))* VDC *(sin(3*x(1))-sin(3*x(2))+sin(3*x(3)))-V3; 
           (4/(5*pi))* VDC *(sin(5*x(1))-sin(5*x(2))+sin(5*x(3)))-V5]; 
alphavalues = fsolve(f,[pi/6 pi/6 pi/6]) * 180/pi
