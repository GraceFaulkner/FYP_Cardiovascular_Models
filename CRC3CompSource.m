%% 3 Compartment CRC Windkessel (Rod's current option) with Coupled LV and resistive valves 
%  Analytic Shi double cosine elastance and derivative elastance, to avoid
%  numerical instability 

clear all
close all

%Model parameters
global C_v_h C_v_l C_v_t C_a_tot R_tot rmv rav HR tauS1 tauS2 ELVmax ELVmin t1 t2 t3 t4 grad_l grad_h grad_t

%Initial values and component params
mcfp        = 50;
P0          = [mcfp mcfp mcfp];                                            % mean circulatory filling pressure 
R_h         = 1.11;                                                        % Systemmic windkessel
C_a_h       = 1.1;                                                         % NB, Harry / Nicholai values R = 1.11, C = 1.1  
C_v_h       = 10; 
R_l         = 1.11;                                                        % Systemmic windkessel
C_a_l       = 1.1;                                                         % NB, Harry / Nicholai values R = 1.11, C = 1.1  
C_v_l       = 10;
R_t         = 1.11;
C_a_t       = 1.1;
C_v_t       = 10;
rav         = 0.01;                                                        % Valves 
rmv         = 0.3; 

tauS1       = 0.35;                                                        % Shi elastance 
tauS2       = 0.50;
ELVmax      = 0.10;
ELVmin      = 0.05;  
HR          = 60;

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
LBNP        = -60;
MBNP        = -20;
UBPP        = 20;
grad_l      = LBNP/(t2-t1);                                                 % mmHg / sec
grad_h      = UBPP/(t2-t1);
grad_t      = MBNP/(t2-t1);

%% Solve
options       = odeset('RelTol', 1e-06, 'AbsTol', 1e-06, 'MaxStep', 0.05);
[t,P]=ode45(@(t,P)Pode(t,P),tspan,P0,options);

%compute absolute pressure
P_static_l   = 0                     .* (t<=t1)            +...             %Leg compartment
              (t-t1) .* grad_l        .* (t> t1) .* (t<=t2) +...
              LBNP                  .* (t> t2) .* (t<=t3) +...
              (t-t4) .*-grad_l        .* (t> t3) .* (t<=t4) +...
              0                     .* (t> t4);
P_static_l   = smooth(smooth(P_static_l));                              
P_static_h   = 0                     .* (t<=t1)            +...             %Head compartment
              (t-t1) .* grad_h        .* (t> t1) .* (t<=t2) +...
              UBPP                  .* (t> t2) .* (t<=t3) +...
              (t-t4) .*-grad_h        .* (t> t3) .* (t<=t4) +...
              0                     .* (t> t4);
P_static_h   = smooth(smooth(P_static_h));
P_static_t   = 0                     .* (t<=t1)            +...             %Thoracic compartment
              (t-t1) .* grad_t        .* (t> t1) .* (t<=t2) +...
              MBNP                  .* (t> t2) .* (t<=t3) +...
              (t-t4) .*-grad_t        .* (t> t3) .* (t<=t4) +...
              0                     .* (t> t4);
P_static_t   = smooth(smooth(P_static_t));

%%plot the source presures
% figure()
% hold on; grid on;
% plot(t, P_static_h)
% plot(t, P_static_l)
% plot(t, P_static_t)


%Absolute pressures
P_abs_h    = P(:,3) -P_static_h;
P_abs_l    = P(:,3) -P_static_l;
P_abs_t    = P(:,3) -P_static_t;

%%plot abs pressures
% figure()
% hold on; grid on;
% plot(t, P_abs_h)
% plot(t, P_abs_l)
% plot(t, P_abs_t)

%compute volumes
Elas    = elastance2(t);
V_lv    = P(:,1)./Elas;

V_a_h   = C_a_h .* P(:,2);                                                  %head compartment
V_v_h   = C_v_h .* P_abs_h;
V_s_h   = P_static_h * C_v_h;

V_a_l   = C_a_l .* P(:,2);                                                  %leg compartment
V_v_l   = C_v_l .* P_abs_l;
V_s_l   = P_static_l * C_v_l;

V_a_t   = C_a_t .* P(:,2);                                                  %thoracic compartment
V_v_t   = C_v_t .* P_abs_t;
V_s_t   = P_static_t * C_v_t;

V_tot   = V_lv + V_v_h + V_a_h  + V_a_l + V_v_l + V_a_t + V_v_t;            %total

%plots
figure('name','signals1'); hold on; grid on; 
    plot(t,P(:,1)); 
    plot(t,P(:,2)); 
    plot(t,P(:,3)); 
    plot(t,P_abs_h); 
    plot(t,P_abs_l); 
    plot(t,P_abs_t); 
    xlabel('t (secs)'); ylabel('p(t) (mmHg)'); 
    legend('P_L_V','P_A_O','P_V_E_N','P_A_B_S_h','P_A_B_S_l', 'P_A_B_S_t'); 
    title('Pressures')

figure('name', 'volumes');hold on; grid on;
    plot(t, V_lv); 
    plot(t, V_a_h);
    plot(t, V_v_h);
    plot(t, V_a_l);
    plot(t, V_v_l);
    plot(t, V_a_t);
    plot(t, V_v_t);
    plot(t, V_tot);
    xlabel('t (secs)'); ylabel('V(t) (m^3)');
    legend('V_L_V','V_A_O_h','V_V_E_N_h','V_A_O_l','V_V_E_N_l','V_A_O_t','V_V_E_N_t', 'V_T_O_T');
    title('Volumes')


%% Functions

function [dPdt] = Pode(t,P)
global C_v_h C_v_l C_v_t C_a_tot R_tot rmv rav t1 t2 t3 t4 grad_l grad_h grad_t
t=t
[Elocal, dEdtlocal] =  elastance2 (t);                                      %elastance values
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
dPdt(3)   = (1/(C_v_h + C_v_l + C_v_t))*(-iMV...                                    % Venous pressure
          + (1/R_tot)*P(2) - (1/R_tot)*P(3)...
          + C_v_l*sourcel + C_v_h*sourceh + C_v_t*sourcet);
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

function  [wk1, wk2] = elastance2(t)
% Shi double cosine elastance and derivative elastance 
global HR tauS1 tauS2 ELVmax ELVmin
t      = rem(t,60/HR);
tau    = t * HR / 60;
e      = (1-cos( tau       /tauS1        *pi)).*(tau<=tauS1)              +...
         (1+cos((tau-tauS1)/(tauS2-tauS1)*pi)).*(tau> tauS1).*(tau<=tauS2)+...
          0                                                  .*(tau> tauS2);
wk1    = ELVmax.*e + ELVmin;
dedtau = sin( tau       /tauS1        *pi)*pi/tauS1        .*(tau<=tauS1)              +...
        -sin((tau-tauS1)/(tauS2-tauS1)*pi)*pi/(tauS2-tauS1).*(tau> tauS1).*(tau<=tauS2)+...
         0                                                               .*(tau> tauS2);
wk2    = ELVmax * dedtau * HR / 60;
end



