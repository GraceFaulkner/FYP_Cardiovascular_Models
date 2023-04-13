% Assign global definitions    

clear all
close all

global driver ddriver test tRange
    
    test =  0;
    
    % Parameters from Kind et al. 2010
    % param = [Rp,    C,   Rd,   L    ]
    hypo    = [0.02,  2,   0.6,  0.005]; 
    normo   = [0.033, 1.5, 0.95, 1.4  ];
    hyper   = [0.05,  0.7, 1.4,  0.02 ];
    
    %local variables
    I0   = 530; 
    HR   = 60;
    tau  = 60/ HR;
    dt   = 0.001;
    
    % assign ODE45 input
    tRange = [0:dt:10];
    Y0     = [0; 0];                                                       % mcfp
    
    %definition of variable flow
    driver = [];

    i       = sin(10/3*pi*tRange/tau);
    i       = i.*(i>0);
    i       = i.*i;
    driver  = I0 *i;
    ddriver = diff(driver)/dt; ddriver= [ddriver,ddriver(end)];

    %define options for ode solver
    options        = odeset('Reltol',1e-6);

    %solve ODE for states from Kind et al. 2010
    [tHypoSol,  hypoSol ] = ode45(@(t,Y,param) FourEWM(t,Y, hypo ), tRange, Y0, options);     
    [tNormoSol, normoSol] = ode45(@(t,Y,param) FourEWM(t,Y, normo), tRange, Y0, options);   
    [tHyperSol, hyperSol] = ode45(@(t,Y,param) FourEWM(t,Y, hyper), tRange, Y0, options);   
    
    %determine P values fr each state (for comparison with literature)
    PHypo   = hypoSol (:,1);
    PNormo  = normoSol(:,1);
    PHyper  = hyperSol(:,1);
    
    %plot
    figure()
     subplot(2,1,1)
        plot(tRange, driver)
        hold on; grid on;
        plot(tRange, ddriver)
        hold off
        legend("Flow", "Flow derivative")
        xlabel("Time (s)")
        ylabel("Flow (ml/s)")
        title("Flow driver and derivative")
     subplot(2,1,2)
        hold on; grid on;
        plot(tHypoSol, PHypo)
        plot(tNormoSol, PNormo)
        plot(tHyperSol, PHyper)
        legend("Hypotensive", "Normotensive", "Hypertensive")
        xlabel("Time (s)")
        ylabel("Pressure (mmHg)")
        title("Pressure response")
        hold off

        
     %% functions
     function dYdt = FourEWM(t, Y, param)
        global driver ddriver test tRange
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
        dPdt  = (Rp/L - 1/(Rd*C))* PP - Rp/L * P + 1/C * driverval + Rp * ddriverval;
        dPPdt = - 1/(Rd*C)       * PP + 0    * P + 1/C * driverval + 0  * ddriverval;
        
        %report the output
        dYdt = [dPdt; dPPdt];
     end
