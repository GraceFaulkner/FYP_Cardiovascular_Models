%% Compute values from 2EWM using parameter values from Kind et al. 
% Assign global definitions    

clear all
close all

% Assign parameter values
global driver ddriver test tRange dt 
    
test =  0;
I0   = 530; 
HR   = 65;
tau  = 60/ HR;
dt   = 0.001;
   
% Parameters from Kind et al. 2010
% param = [Rp,    C,   Rd,   L     I0     tau  ]
hypo    = [0.02,  2,   0.6,  0.005, 530,  tau ]; 
normo   = [0.033, 1.5, 0.95, 1.4,   530,  tau ];
hyper   = [0.05,  0.7, 1.4,  0.02,  530,  tau ];
    
tRange = [0:dt:20];
Y0     = [50; 50];                             % mcfp

i       = sin(2*pi*tRange/tau);
i       = i.*(i>0);
i       = i.*i;
driver  = I0 *i;
ddriver = diff(driver)/dt; ddriver= [ddriver,ddriver(end)];

param   = normo;

%Compute ODEs
options  = odeset('Reltol',1e-6);
[time,p ] = ode45(@(t,Y,param) Pode(t,Y, normo ), tRange, Y0, options);

%plot pressure
figure("Name", "3EWM Pressure-Time Waveform")
plot(time, p)
title("3EWM Pressure-Time Waveform")
xlabel("Time (s)")
ylabel("Pressure (mmHg)")
legend("Aortic Pressure", "Arterial pressure")


% derived metrics
Nt        = length(tRange);
Nt_upper  = round(0.75 * Nt); 
t_prime   = time (Nt_upper:end  );
P         = p    (Nt_upper:end,1);
PP        = p    (Nt_upper:end,2);
Psys      = max  (P) ;
Pdias     = min  (P) ;
Pmean     = mean (P) ;
PPsys     = max  (PP);
PPdias    = min  (PP);
PPmean    = mean (PP);

ss        = [];
for i = 1:length(time)
    if p(i)/Psys > 0.95
        tester = p(i)/Psys
        ss(end+1) = time(i);
    end
end
time2ss   = ss(1);

figure("Name", "3EWM Pressure-Time Waveform")
plot(t_prime, P)
hold on;
plot(t_prime, PP)
title("Steady State 3EWM Pressure-Time Waveform")
xlabel("Time (s)")
ylabel("Pressure (mmHg)")
legend("P", "PP")



%%Function
 function dYdt = Pode(t, Y, param)
   global driver ddriver test tRange dt
    % split the initial conditions
    P  = Y(1);
    PP = Y(2);

    %split the parameters
    Rp = param(1);
    C  = param(2);
    Rd = param(3);
    L  = param(4);

    test = test + 1;
    %find the input at this time
    index      = min(find(t<=tRange));
    driverval  = driver (index);
    ddriverval = ddriver(index);
    
    % calculate the output
    dPdt  = 1/C * driverval - 1/(Rd*C) * PP + Rp * ddriverval;
    dPPdt = 1/C * driverval - 1/(Rd*C) * PP ;
    
    %report the output
    dYdt = [dPdt; dPPdt];
 end
     
     
     

