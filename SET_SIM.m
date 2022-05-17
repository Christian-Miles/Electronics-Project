%% EGH448-2021 SIMULINK M-FILE
clc
clear all
close all

%% SETTINGS (THE FOLLOWING VALUES MUST BE SET TO THE VALUES WHICH ARE GIVEN FOR THE STUDENT GROUPS IN TABLE 4 OF THE INSTRUCTION MANUAL.)
VDC=24;   % DC-LINK VOLTAGE
V1=15.7656;	% THE FIRST HARMONIC OF THE OUTPUT VOLTAGE
V3=0.5912;	% THE THIRD HARMONIC OF THE OUTPUT VOLTAGE
V5=2.1877;		% THE FIFTH HARMONIC OF THE OUTPUT VOLTAGE
RL=33;   % RESISTANCE OF THE SERIES R-L LOAD
LL=1.2e-3; % INDUCTANCE OF THE SERIES R-L LOAD
%% SETTINGS (HERE THE CODES RELATED TO THE NON-LINEAR SET OF EQUATIONS SHOULD BE WRITTEN AND a1, a2, AND a3 SHOULD BE OBTAINED FROM IT USING THE GIVEN V1, V3, V5 AND VDC VALUES FOR EACH GROUP)
a(1)=18.0003; % α1 (THE FIRST ANLGE OF THE THREE-ANLGE CONTROL SCHEME)
a(2)=39.6000; % α2 (THE SECOND ANLGE OF THE THREE-ANLGE CONTROL SCHEME)
a(3)=57.6002; % α3 (THE THIRD ANLGE OF THE THREE-ANLGE CONTROL SCHEME)

%% THESE VALUES IS GIVEN IN THE PROVIDED MATLAB AND SIMULINK FILES. (YOU ARE REQUIRED TO SET F TO THE VALUE GIVEN IN TABLE 4 OF THE INSTRUCTION MANUAL.)
f=500;       % FUNDAMENTAL FREQUENCY
w=2*pi*f;   % FUNDAMENTAL FREQUENCY BASED ON RADIANT, FOR FUTURE CODES
T=20/f;     % TIME OF THE SIMULATION (SET TO TAKE 20 FULL-CYCLES)
T_SMP=2e-6; % STEP SIZE OF THE SIMULATION

%% TO RUN THE SIMULINK FILE
%Overwrite variables with those found in instruction manual for testing
if exist('testing', 'var') && testing
    VDC = 8.33; %Testing
    V1 = 5.893; %Testing
    V3 = 1.572; %Testing
    V5 = 0; %Testing
    a = [21.99, 36.14, 50.44]; %Testing
    RL = 3.335; %Testing
    LL = 6.670e-3; %Testing
    f = 50; %Testing
    w=2*pi*f;   % FUNDAMENTAL FREQUENCY BASED ON RADIANT, FOR FUTURE CODES
    T=20/f;     % TIME OF THE SIMULATION (SET TO TAKE 20 FULL-CYCLES)
    T_SMP=2e-6; % STEP SIZE OF THE SIMULATION
end

run('SHM_SIM.slx')  % THIS CODE ONLY OPENS THE SIMULINK FILE
sim('SHM_SIM.slx')  % THE NECESSARY CODE TO RUN THE SIMULINK FILE


%% USER VARIABLES FOR FOLLOWING SECTIONS
harmonic_num = [1:2:15]';
a_rad = pi/180 .* a;

%% OUTPUT VOLTAGE HARMONICS SHOULD BE CALCULATED HERE
t = linspace(0, T, T/T_SMP+1);
V0_FT = 4/pi * VDC * (1./harmonic_num .* (sin(harmonic_num.*a_rad(1)) - sin(harmonic_num.*a_rad(2)) + sin(harmonic_num.*a_rad(3))));
V0_FT_ANG = -sign(V0_FT) * 180;
V0_FT_ANG(V0_FT_ANG==-180) = 0;

%% OUTPUT VOLTAGE TIME-DOMAIN HARMONIC APPROXIMATION CODES SHOULD BE WRITTEN HERE
V0_FT_TIME = 4/pi * VDC * (1./harmonic_num .* (sin(harmonic_num.*a_rad(1)) - sin(harmonic_num.*a_rad(2)) + sin(harmonic_num.*a_rad(3)))) .* cos(harmonic_num*2*pi*f*t);
V0_FT_TIME_1 = sum(V0_FT_TIME([1],:), 1)';
V0_FT_TIME_3 = sum(V0_FT_TIME([1:3],:), 1)';
V0_FT_TIME_8 = sum(V0_FT_TIME([1:8],:), 1)';

hold off;
plot(t-t(end)+1/f,vo);
hold on;
plot(t-t(end)+1/f,V0_FT_TIME_1);
plot(t-t(end)+1/f,V0_FT_TIME_3);
plot(t-t(end)+1/f,V0_FT_TIME_8);
xlim([0 1/f]);

%% OUTPUT VOLTAGE RMS AND THD CALCULATIONS SHOULD BE WRITTEN HERE
VRMS = VDC*sqrt(2/pi * (a_rad(1) - a_rad(2) + a_rad(3)));
OUTPUT_VOLTAGE_THD = sqrt(2*abs(VRMS)^2 - abs(V0_FT(1))^2)/abs(V0_FT(1)) * 100;

%% OUTPUT CURRENT HARMONICS SHOULD BE CALCULATED HERE
Z = sqrt(RL^2 + (harmonic_num*2*pi*f*LL).^2);
Z_ANG_RAD = atan2((harmonic_num*2*pi*f*LL),RL);
Z_ANG = Z_ANG_RAD * 180/pi;

I = V0_FT./Z;
I_ANG_RAD = atan2((harmonic_num*2*pi*f*LL), RL);
I_ANG = I_ANG_RAD * 180/pi;

%% OUTPUT CURRENT TIME-DOMAIN HARMONIC APPROXIMATION CODES SHOULD BE WRITTEN HERE
PHI = -sign(V0_FT_TIME) * 180;
PHI(PHI<0) = 0;
I_TIME = (V0_FT_TIME./Z) .* cos(harmonic_num*2*pi*f*t + (PHI*pi/180) - atan((harmonic_num*2*pi*f*LL)/(RL)));
I_TIME_1 = sum(I_TIME([1], :), 1);
I_TIME_3 = sum(I_TIME([1:3], :), 1);
I_TIME_8 = sum(I_TIME([1:8], :), 1);

%% OUTPUT CURRENT RMS AND THD CALCULATIONS SHOULD BE WRITTEN HERE
IRMS = sqrt(sum(abs(I).^2/2));
OUTPUT_CURRENT_THD = sqrt(2*abs(IRMS)^2 - abs(I(1))^2)/abs(I(1)) * 100;

hold off;
plot(t-t(end)+1/f,io);
hold on;
plot(t-t(end)+1/f,I_TIME_1);
plot(t-t(end)+1/f,I_TIME_3);
plot(t-t(end)+1/f,I_TIME_8);
xlim([0 1/f]);

%% POWER CALCULATIONS SHOULD BE WRITTEN HERE
PL_CAL = RL*(IRMS^2);
SL_CAL = VRMS*IRMS;
PL_SIM = RL*(rms(VI_O.signals(2).values)^2);
SL_SIM = rms(VI_O.signals(1).values)*rms(VI_O.signals(2).values);

PF_CAL = PL_CAL/SL_CAL;
PF_SIM = PL_SIM/SL_SIM;

%% SOME USEFUL CODES
% EXTRACTING TIME-DOMAIN INFORMATION FROM THE STRUCTURE VARIABLE RECORDED BY THE SIMULINK FILE (SCOPE) INTO THE WORKSPACE
t=VI_O.time;				% EXTRACTING TIME ARRAY
vo=VI_O.signals(1).values;		% EXTRACTING OUTPUT VOLTAGE
io=VI_O.signals(2).values;		% EXTRACTING OUTPUT CURRENT

% FFT CALCULATIONS FOR vo
SAMP=t(2)-t(1);
fs=1/SAMP;
n=length(t);
fr=[0:n-1]*fs/(n);
VO=2*(fft(vo))/n;

% PLOTTING THE OBTAINED FFT FIGURE
figure(4481)
plot(fr,abs(VO),'LineWidth',2)
xlim([0 f*15+150])
grid on
xlabel('Frequency (Hz)')
ylabel('|V_O|')
title('Frequency Content')


% PLOTTING THE OBTAINED TIME-DOMAIN OUTPUT VOLTAGE 
figure(4482)
plot(t-t(end)+1/f,vo,'LineWidth',2)
grid on
xlim([0 1/f])
xlabel('Time (sec)')
ylabel('Voltage (V)')
title('Output Voltage (v_o)')