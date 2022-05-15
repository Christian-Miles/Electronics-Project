%% EGH448-2021 SIMULINK M-FILE
clc
clear all
close all

%% SETTINGS (THE FOLLOWING VALUES MUST BE SET TO THE VALUES WHICH ARE GIVEN FOR THE STUDENT GROUPS IN TABLE 4 OF THE INSTRUCTION MANUAL.)
VDC=8.33;   % DC-LINK VOLTAGE
V1=5.892;	% THE FIRST HARMONIC OF THE OUTPUT VOLTAGE
V3=1.572;	% THE THIRD HARMONIC OF THE OUTPUT VOLTAGE
V5=0;		% THE FIFTH HARMONIC OF THE OUTPUT VOLTAGE
RL=3.335;   % RESISTANCE OF THE SERIES R-L LOAD
LL=6.67e-3; % INDUCTANCE OF THE SERIES R-L LOAD
%% SETTINGS (HERE THE CODES RELATED TO THE NON-LINEAR SET OF EQUATIONS SHOULD BE WRITTEN AND a1, a2, AND a3 SHOULD BE OBTAINED FROM IT USING THE GIVEN V1, V3, V5 AND VDC VALUES FOR EACH GROUP)
a(1)=21.99; % α1 (THE FIRST ANLGE OF THE THREE-ANLGE CONTROL SCHEME)
a(2)=36.14; % α2 (THE SECOND ANLGE OF THE THREE-ANLGE CONTROL SCHEME)
a(3)=50.44; % α3 (THE THIRD ANLGE OF THE THREE-ANLGE CONTROL SCHEME)

%% THESE VALUES IS GIVEN IN THE PROVIDED MATLAB AND SIMULINK FILES. (YOU ARE REQUIRED TO SET F TO THE VALUE GIVEN IN TABLE 4 OF THE INSTRUCTION MANUAL.)
f=50;       % FUNDAMENTAL FREQUENCY
w=2*pi*f;   % FUNDAMENTAL FREQUENCY BASED ON RADIANT, FOR FUTURE CODES
T=20/f;     % TIME OF THE SIMULATION (SET TO TAKE 20 FULL-CYCLES)
T_SMP=2e-6; % STEP SIZE OF THE SIMULATION

%% TO RUN THE SIMULINK FILE
run('SHM_SIM.slx')  % THIS CODE ONLY OPENS THE SIMULINK FILE
sim('SHM_SIM.slx')  % THE NECESSARY CODE TO RUN THE SIMULINK FILE


%% OUTPUT VOLTAGE HARMONICS SHOULD BE CALCULATED HERE

%% OUTPUT VOLTAGE TIME-DOMAIN HARMONIC APPROXIMATION CODES SHOULD BE WRITTEN HERE

%% OUTPUT VOLTAGE RMS AND THD CALCULATIONS SHOULD BE WRITTEN HERE

%% OUTPUT CURRENT HARMONICS SHOULD BE CALCULATED HERE

%% OUTPUT CURRENT TIME-DOMAIN HARMONIC APPROXIMATION CODES SHOULD BE WRITTEN HERE

%% OUTPUT CURRENT RMS AND THD CALCULATIONS SHOULD BE WRITTEN HERE

%% POWER CALCULATIONS SHOULD BE WRITTEN HERE

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