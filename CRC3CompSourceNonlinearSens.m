%% 3 Compartment CRC Windkessel (Rod's current option) with Coupled LV and resistive valves 
%  Analytic Shi double cosine elastance and derivative elastance, to avoid
%  numerical instability 

clear all
close all

%Model parameters
global C_v_h C_v_l C_v_t C_a_h C_a_l C_a_t C_a_tot R_tot rmv rav HR tauS1 tauS2 ELVmax ELVmin t1 t2 t3 t4 grad_l grad_h grad_t C0 Vmax tspan Sol0

%Initial values and component params
mcfp        = 50;
P0          = [mcfp mcfp mcfp];                                            % mean circulatory filling pressure 
R_h         = 3.9;                                                         % Heldt upper body arterial resistance
C_a_h       = 2/3;                                                         % Heldt arterial resistance over 3 compartments 
C_v_h       = 8;                                                           % Heldt
R_l         = 2.4;                                                         % From Heldt total arterial resistance = 1 PRU
C_a_l       = 2/3;                                                         % Heldt arterial resistance over 3 compartments 
C_v_l       = 19;                                                          % Heldt
R_t         = 3;                                                           % Heldt splanchnic arterial resistance
C_a_t       = 2/3;                                                         % Heldt arterial resistance over 3 compartments 
rav         = 0.01;                                                        % Heldt Rlo
rmv         = 0.03;                                                        % Combined resistance of veins and lungs

HR          = 60;
tau         = HR/60;
tauS1       = 1/ 3 *tau;                                                   % Time parameters set to match normalised Heldt values
tauS2       = 1.5 / 3 * tau;
ELVmax      = 0.15;
ELVmin      = 0.1;  

%Artial capacitance and proximal resistance
C_a_tot     = C_a_h + C_a_l + C_a_t;
R_tot       = (1/R_h + 1/R_l +1/R_t)^(-1);

%Source term parameters                                      
T_max       = 30;
tspan       = [0 T_max];
t1          = 10;                                                           % Bias parameters 
t2          = 15;
t3          = 25;
t4          = 30;
LBNP        = 0;%-60;
MBNP        = 0;%-20;
UBPP        = 0;%20;
grad_l      = LBNP/(t2-t1);                                                 % mmHg / sec
grad_h      = UBPP/(t2-t1);
grad_t      = MBNP/(t2-t1);


C0          = 70;                                                              % Time varying capacitance values from heldt
Vmax        = 1450;

para        = [R_h, R_l, R_t, C_a_h, C_a_l, C_a_t, C_v_h, C_v_l, C0, Vmax, rav, rmv, T_max, HR, tauS1, tauS2, ELVmax, ELVmin];

%% Solve

Sol0        = Compute_Metrics(para, P0);
Sens5       = Sensitivity(para, P0, 1.05);
Sens10      = Sensitivity(para, P0, 1.1);
Sens30      = Sensitivity(para, P0, 1.3);

%plot sensitvity
figure('name','sensitivity')
xvalues     = {'R_h', 'R_l', 'R_t', 'C_a_h', 'C_a_l', 'C_a_t', 'C_v_h', ...
               'C_v_l', 'C0', 'Vmax', 'rav', 'rmv', 'T_max', 'HR',  ...
               'tauS1', 'tauS2', 'Elvmax', 'Elvmin'};
yvalues     = {'PLVsys', 'PLVdias', 'PLVmean', 'Partsys', 'Partdias', ...
               'Partmean', 'Pvensys', 'Pvendias', 'Pvenmean', 'Volmax', 'Volmin', 'Volmean'};
subplot(3,1,1)
h           = heatmap(xvalues, yvalues, abs(Sens5));
title('Sensitivity: 5%')
subplot(3,1,2)
i           = heatmap(xvalues, yvalues, abs(Sens10));
title('Sensitivity: 10%')
subplot(3,1,3)
j           = heatmap(xvalues, yvalues, abs(Sens30));
title('Sensitivity: 30%')


figure('name', 'orthogonality')
Orth        = Orthogonality(para, P0);
values      = {'R_h', 'R_l', 'R_t', 'C_a_h', 'C_a_l', 'C_a_t', 'C_v_h', ...
               'C_v_l', 'C0', 'Vmax', 'rav', 'rmv', 'T_max', 'HR',  ...
               'tauS1', 'tauS2', 'Elvmax', 'Elvmin'};
k           = heatmap(values, values, Orth);



%% Functions
function Orth = Orthogonality(para, IC)
    theta0 = para;
    Y0     = IC;

    Sens30 = Sensitivity(theta0, Y0, 1.3)
    for i = 1:length(theta0)
        for j = 1:length(theta0)
            a         = Sens30(:,i);
            b         = Sens30(:,j);
            calc      = dot(a,b)/(norm(a)*norm(b));
            rcalc     = round(calc,4);
            angle     = acosd(rcalc);
            Orth(i,j) = angle;
        end
    end
end

function Sens = Sensitivity(para, IC, percent)
    global Sol0

    theta0 = para;
    Y0     = IC;
    Sens = zeros(length(Sol0), length(theta0));

    for i = 1:length(theta0)
        theta     = theta0;
        theta(i)  = theta0(i) * percent;
        Sol       = Compute_Metrics(theta, Y0);
        Change    = (Sol - Sol0) ./ (Sol0) * percent;
        Sens(:,i) = Change;
    end
end

function Sol = Compute_Metrics(para, IC)
    tspan = [0 para(13)];
    P0    = IC;

    options       = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 0.05);
    [t,P]=ode45(@(t,P)Pode(t,P,para),tspan,P0,options);

    Nt          = length(t);
    Nt_upper    = round(0.8*Nt);
    t_prime     = t (Nt_upper:end);
    PLV         = P (Nt_upper:end, 1);
    Part        = P (Nt_upper:end, 2);
    Pven        = P (Nt_upper:end, 3);

    PLVsys      = max(PLV);
    PLVdias     = min(PLV);
    PLVmean     = mean(PLV);

    Partsys     = max(Part);
    Partdias    = min(Part);
    Partmean    = mean(Part);

    Pvensys     = max(Pven);
    Pvendias    = min(Pven);
    Pvenmean    = mean(Pven);

    %compute volumes
    Vmax      = para(10);
    C0        = para(9);
    C_v_h     = para(7);
    C_v_l     = para(8);
    R_h       = para(1);
    R_l       = para(2);
    R_t       = para(3);
    C_a_h     = para(4);
    C_a_l     = para(5);
    C_a_t     = para(6);
    rav       = para(11);
    rmv       = para(12);
    epara       = para(14:18);
    [Elas, dEdt]= elastance2(t_prime,epara);
    V_lv    = PLV./Elas;

    V_a_h   = C_a_h .* Part;                                                  %head compartment
    V_v_h   = C_v_h .* Pven;

    V_a_l   = C_a_l .* Part;                                                  %leg compartment
    V_v_l   = C_v_l .* Pven;


    Ptrans  = Pven;
    alpha   = 2* Vmax / pi;
    V0      = 0;                                                                     %assumption
    C_v_s   = 1./Ptrans .* (V0 + alpha .* atan(C0/alpha .*Ptrans) );

    V_a_t   = C_a_t .* Part;                                                  %thoracic compartment
    V_v_t   = C_v_s .* Pven;

    V_tot   = V_lv + V_v_h + V_a_h  + V_a_l + V_v_l + V_a_t + V_v_t;            %total
    
    V_max   = max(V_tot);
    V_min   = min(V_tot);
    V_mean  = mean(V_tot);


    Sol     = [PLVsys; PLVdias; PLVmean; Partsys; Partdias; Pvenmean; Pvensys; Pvendias; Pvenmean; V_max; V_min; V_mean];
end



function [dPdt] = Pode(t,P, para)
global t1 t2 t3 t4 grad_l grad_h grad_t
Vmax      = para(10);
C0        = para(9);
C_v_h     = para(7);
C_v_l     = para(8);
R_h       = para(1);
R_l       = para(2);
R_t       = para(3);
C_a_h     = para(4);
C_a_l     = para(5);
C_a_t     = para(6);
rav       = para(11);
rmv       = para(12);
epara     = para(14:18);

%Artial capacitance and proximal resistance
C_a_tot     = C_a_h + C_a_l + C_a_t;
R_tot       = (1/R_h + 1/R_l +1/R_t)^(-1);

t=t
[Elocal, dEdtlocal] =  elastance2 (t, epara);                                      %elastance values
iAV       = valve(P(1)-P(2),rav);                                           %valve equations
iMV       = valve(P(3)-P(1),rmv);

sourceh   = ( heaviside(t-t1   ) - heaviside(t-t2   ) ) * grad_h +...       %Head compartment pressure source
          - ( heaviside(t-t3) - heaviside(t-t4) ) * grad_h; 
sourcel   = ( heaviside(t-t1   ) - heaviside(t-t2   ) ) * grad_l +...       %Leg compartment pressure source
          - ( heaviside(t-t3) - heaviside(t-t4) ) * grad_l; 
sourcet   = ( heaviside(t-t1   ) - heaviside(t-t2   ) ) * grad_t +...       %Thoracic compartment pressure source
         - ( heaviside(t-t3) - heaviside(t-t4) ) * grad_t; 

dPdt(1)   = (iMV-iAV)*Elocal + P(1)/Elocal*dEdtlocal;                       % LV 
dPdt(2)   = (1/C_a_tot)*(iAV...                                             % Artery pressure
          - (1/R_tot)*P(2) + (1/R_tot)*P(3));

alpha     = 2* Vmax / pi;
Ptrans    = P(3) - sourcet;
beta      = C0/(1+(C0/alpha*Ptrans)^2);

dPdt(3)   = 1/(C_v_h + C_v_l + beta)*(-iMV...                                    % Venous pressure
          + (1/R_tot)*P(2) - (1/R_tot)*P(3)...
          + C_v_l*sourcel + C_v_h*sourceh + beta*sourcet);
dPdt      = dPdt';
end

function q = valve (dP,Res)
    %Valve function
    if dP > 0
        wk1 = dP / Res; 
    else
        % reverse bias resistance very large
        wk1 = dP / 10000 / Res; 
    end
    q = wk1;
end

function  [wk1, wk2] = elastance2(t, epara)
% Shi double cosine elastance and derivative elastance 
HR     = epara(1);
tauS1  = epara(2);
tauS2  = epara(3);
ELVmax = epara(4);
ELVmin = epara(5);


t      = rem(t,60/HR);
tau    = t * HR / 60;
e      = 0.5* (1-cos( tau       /tauS1        *pi)).*(tau<=tauS1)              +...
         0.5* (1+cos((tau-tauS1)/(tauS2-tauS1)*pi)).*(tau> tauS1).*(tau<=tauS2)+...
          0                                                  .*(tau> tauS2);
wk1    = ELVmax.*e + ELVmin;
dedtau = 0.5* sin( tau       /tauS1        *pi)*pi/tauS1        .*(tau<=tauS1)              +...
        -0.5* sin((tau-tauS1)/(tauS2-tauS1)*pi)*pi/(tauS2-tauS1).*(tau> tauS1).*(tau<=tauS2)+...
         0                                                               .*(tau> tauS2);
wk2    = ELVmax * dedtau * HR / 60;
end



