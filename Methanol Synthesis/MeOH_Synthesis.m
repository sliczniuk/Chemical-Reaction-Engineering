% ICO - ese02

% MeOH synthesis Reactor Design: Adiabatic Stages with 1D-PFR pseudohom

clc
close all
clearvars

global R NS NR Tc Pc om MM nu cp_mix DhR P
global eps rho_cat Tlim xlim
global Mtot

%% Data Declaration

% We have the universal gas constant:
R  = 8.3144621      ;%[Pa*m^3/(mol*K) = J/(mol*K)]

% We list the number of species in the reactions:
NS = 6;

% We list the number of reactions:
NR = 2;

% We list the species' properties:
Tc = [132.92; 304.19; 33.18; 647.13; 512.58; 190.56]            ;%[K]
      %CO      %CO2     %H2     %H2O    %CH3OH    %CH4    
Pc = [ 34.99;  73.82; 13.13; 220.55;  80.96;  45.90]*1e5        ;%[Pa]
      %CO      %CO2     %H2     %H2O    %CH3OH    %CH4
om = [ 0.066;  0.228;-0.220;  0.345;  0.566;  0.011]            ;%[-]
      %CO      %CO2     %H2     %H2O    %CH3OH    %CH4
MM = [28.0104; 44.0098; 2.0159; 18.0153; 32.0422; 16.0428]*1e-3 ;%[kg/mol]
      %CO      %CO2     %H2     %H2O    %CH3OH    %CH4
      
% Due to the reactions, we write the stoichometric matrix:
nu = [-1      0      -2       0       1       0         %R1 (MeOH)
       1     -1      -1       1       0       0   ];    %R2 (WGS)
       %CO    %CO2	  %H2     %H2O	 %CH3OH   %CH4

% We express the catalyst void fraction:
eps = 0.4                   ;%[-]

% We have the catalyst density:
rho_cat = 1.98e3            ;%[kg_cat/m3]

% We have the mixture density:
cp_mix = 4.081e3            ;%[J/kg]

% We have the reaction enthalpies:
DhR = [-23460; 9510]*4.184  ;%[J/mol]

% We have the operative pressure:
P = 51e5                    ;%[Pa]

% We have the limit temperature for each bed:
Tlim = 270 + 273.15         ;%[K]

% We have the upper limit of the molar fraction of CH3OH
xlim = 0.061;
     
%% Simulation - Pre Processing

% We have the flow state entering the reactor:
Ftotin = 15000/3.6      ;%[mol/s]
%
Tin = 239.85 + 273.15   ;%[K]
%
xin = [0.1339; 0.0997; 0.6333; 0.0000; 0.0000; 0.1331];
       %CO     %CO2     %H2     %H2O     %CH3OH     %CH4

% We evaluate the component molar flows:
Fin = Ftotin*xin        ;%[mol/s]

% We evaluate the component mass flows:
Min = Fin.*MM           ;%[kg/s]

% We evaluate the total mass flow:
Mtot = sum(Min)         ;%[kg/s]

% We evaluate the component mass fractions:
omin = xj2omj(xin)      ;%[-]

%% Simulation - Processing

% options for the solver:
options = odeset('RelTol',1e-7,'AbsTol',1e-9,'Events',@fun_event);

% We declare the initial bed number:
Nbed = 0;

% We set a dummy molar fraction:
xout = 0;

% In the beginning, the initial volume is zero:
Vstart = 0;

% We store the integration variables in an array:
Vtot = [];
Xtot = [];
Vbed = [];

fprintf('+++ Begin Iteration +++\n\n');

% We set the cycle to obtain the effective plot over the beds:
while true
    
    % Advance the bed:
    Nbed = Nbed + 1;
    
    % Set the initial conditions (1):
    omin_bed = omin;
    Tin_bed  = Tin;
    
    % Set the initial conditions (2):
    X0 = [omin_bed; Tin_bed];
    
    % Set the volume with an interval of 200 m3:
    Vspan = [Vstart; Vstart+200];%[m3]
    
    % Integrate the ODE system for the bed:
    [V, X, Vval, Xval] = ode15s(@fun_PFRbed,Vspan,X0,options);
    
    % We assign the new starting points for the initial conditions:
    omin = Xval(1:NS)';
    Tin  = 240 + 273.15;

    % We evaluate the molar fractions:
    x_jout = omj2xj(omin);
    
    % We exploit the molar fraction of CH3OH:
    xout = x_jout(5);
    wout = omin(5);
           
    % Check if the outlet x_MeOH is high enough:
    if xout >= xlim
        
        fprintf('+++ End iteration +++\n\n');
        
        % We compile the matrix of molar fractions:
        NX = size(X,1);
        xfin = zeros(NX,NS);
        %
        for j=1:NX
            xfin(j,:) = omj2xj(X(j,1:NS)')';
        end
        
        % We check for the min. index that evaluates the abs to true:
        [~, k] = min(abs(xfin(:,5) - xlim));
        
        % We assign the new values to V and X based on the limit:
        V = V(1:k);
        X = X(1:k,:);
        
        % We can calculate the actual bed volume:
        Vbed(Nbed) = V(end) - Vstart; %#ok<SAGROW>
        
        % We assign the new starting points for the bed:
        Vstart = V(end);
        
        % We evaluate the terminating w and x for MeOH:
        xout = xfin(k,5);
        wout = X(k,5);
        
        % We print the message for indications:
        fprintf('Bed no. %d\n'                ,Nbed);
        fprintf('Reached volume = %.4f m3\n'  ,Vstart);
        fprintf('Bed volume     = %.4f m3\n'  ,Vbed(Nbed));
        fprintf('Reached w_j    = %.5f\n'     ,wout);
        fprintf('Reached x_j    = %.5f\n\n'   ,xout);  
        fprintf('limit x_CH3OH  = %.5f\n'     ,xlim);        
        
        % We update the plot variables:
        Vtot{Nbed} = V; %#ok<SAGROW>
        Xtot{Nbed} = X; %#ok<SAGROW>
        
        % After the procedure is done, we can exit the loop.
        break;
        
    else
        % We can calculate the actual bed volume:
        Vbed(Nbed) = Vval - Vstart; %#ok<SAGROW>

        % We assign the new starting points for the bed:
        Vstart = Vval;        % We print the message for indications:
    
        fprintf('Bed no. %d\n'                  ,Nbed);
        fprintf('Reached volume = %.4f m3\n'    ,Vstart);
        fprintf('Bed volume     = %.4f m3\n'    ,Vbed(Nbed));
        fprintf('Reached w_j    = %.5f\n'       ,wout);
        fprintf('Reached x_j    = %.5f\n\n'     ,xout);     
    
        % We update the plot variables:
        Vtot{Nbed} = V; %#ok<SAGROW>
        Xtot{Nbed} = X; %#ok<SAGROW>
        
    end    
end


%% Post Processing 
 
% We concatenate all the variables, so the complete reactor is obtained:
Vtot = vertcat(Vtot{:});
Xtot = vertcat(Xtot{:});
 
% We obtain the number of rows (integration steps):
NV = size(Vtot,1);
% 
% We create the reactor matrices:
V_vec = Vtot;

x_mat = zeros(NV,NS);
w_mat = zeros(NV,NS);
T_vec = zeros(NV,1);
 
% We allocate all the variables:
for j=1:NV
%     
    om_vec = Xtot(j,1:NS);          %row vector    
    x_vec  = omj2xj(om_vec');       %column vector
    
    x_mat(j,:) = x_vec';            %row vector   
    w_mat(j,:) = om_vec;            %row vector    
    T_vec(j)   = Xtot(j,NS+1);
    
end
    
figure
plot(V_vec,w_mat(:,5),'LineWidth',1.5);
grid on
xlabel('V [m^3]');
ylabel('CH_3OH mass fraction [-]');
title('Reactor massfrac profile');
saveas(gcf,'02_masfracMeOH_ori.emf');

figure
plot(V_vec,x_mat(:,5),'LineWidth',1.5);
grid on
xlabel('V [m^3]');
ylabel('CH_3OH molar fraction [-]');
title('Reactor molfrac profile');
saveas(gcf,'02_molfracMeOH_ori.emf');

figure
plot(V_vec,T_vec - 273.15,'LineWidth',1.5);
grid on
ylim([235 275]);
xlabel('V [m^3]');
ylabel('T [°C]');
title('Reactor temperature profile');
saveas(gcf,'02_T_ori.emf');

% We evaluate the mass flows:
W_mat = Mtot*w_mat;     %[kg/s]

% We evaluate the molar flows:
F_mat = zeros(NV,NS);
%
for j=1:NV
    F_mat(j,:) = W_mat(j,:).*MM';
end

% We write the CO conversion:
CO_cnv = abs( (F_mat(1,1) - F_mat(:,1))/F_mat(1,1) );

figure
plot(V_vec,CO_cnv,'LineWidth',1.5);
grid on
xlabel('V [m^3]');
ylabel('CO conversion \xi [-]');
title('Reactor CO \xi profile');
saveas(gcf,'02_COcnv_ori.emf');

%% Function Declarations

function f = calc2_DgR(T)% T in [K], dGR in [J/mol]

% This function evaluates the gibbs' free energy of reaction at different
% temperatures.

    global NR
 
    % Coefficient Matrix:
    A = [-22.858  5.602E-2      %R1
           9.418 -0.907E-2];    %R2
    %      %A1    %A2      
    
    I = eye(NR);

    % Temperature Vector:
    Tt = [1; T]; 
    
    % Function as allocated vector:
    f = zeros(NR,1);
    
    for j=1:NR
        
        % Coefficient Vector per reaction:
        At = A'*I(:,j);
        
        % Function value
        f(j) = (At'*Tt)*1e3*4.184; %[J/mol]
        
    end

end

function f = calc_phiPR(T,P,tipo)% T in [K], P in [Pa], phi in [-]

% This function evaluates the activity coefficients per pure specie in the 
% vapour phase with the Peng-Robinson EOS.

    global R Tc Pc om
    
    Tcs = Tc(tipo); Pcs = Pc(tipo); oms = om(tipo);
    
    % Physical variables:
    Som  = 0.37464+1.54226*oms-0.26992*oms^2;
    ktom = (1+Som*(1-sqrt(T/Tcs)))^2;
    
    % Coefficients for the cubic equation:
    a = 0.45724*(R*Tcs)^2*ktom/Pcs;
    b = 0.07780*(R*Tcs)/Pcs;
    A = a*P./(R*T).^2;
    B = b*P./(R*T);
    
    % Solving the cubic equation for the vapour phase:
    Vec    = [1 -1+B A-2*B-3*B.^2 -A*B+B.^2+B.^3];
    VecRis = roots(Vec);
    VecRt  = VecRis(imag(VecRis)==0);
    zliq = min(VecRt); %#ok<NASGU>
    zvap = max(VecRt);
    
    % We declared the vapour phase:
    Z = zvap;
    
    % Residual Gibbs Free Energy:
    gR = R*T*(Z-1-A/(2*sqrt(2)*B)*log((Z+B*(1+sqrt(2)))/(Z+B*(1-sqrt(2))))-log(Z-B));
    
    % Coefficiente di fugacità:
    f = exp(gR/(R*T));

end

function f = omj2xj(om_j)% om_j in [-], x_j in [-]

% This function evaluates the molar fractions given the mass fractions
    
    global MM
    
    % We declare an auxiliary vector:
    vec = om_j./MM;
    
    % We have the molar fractions:
    f = vec/sum(vec);

end

function f = xj2omj(x_j)% x_j in [-], om_j in [-]

% This function evaluates the mass fractions given the molar fractions
    
    global MM
    
    % We declare an auxiliary vector:
    vec = x_j.*MM;
    
    % We have the molar fractions:
    f = vec/sum(vec);

end

function f = calc_rates(om_j,T)% om_j in [-], T in [K], f in [mol/kg_cat/s]

% This function evaluates the reaction rates. Data available from source
% P.L. Villa et al., I&EC (1985), 12.

    global R NS NR P

    % Definition of reaction constants (1):
    AA = [3.49; 2.53; 3.70; 1.54; 5.18];
    
    % Definition of reaction constants (2):
    BB = [4883; -39060; 15948; 8229; 938];
    
    % Explicit coefficients for the rates:
    C = exp(AA + BB*(1/T - 1/506));
    
    % We trasnsform mass fractions in molar fractions:
    x_j = omj2xj(om_j);
    
    % We evaluate the fugacity vector, in [atm]:
    f_j = zeros(NS,1);
    %
    for j=1:NS
        
        f_j(j) = calc_phiPR(T,P,j)*P*x_j(j)/101325;
        
    end

    % This is the order of the species
    % (1)CO, (2) CO2, (3) H2, (4) H2O, (5) CH3OH, (6) CH4

    % We evaluate the equilibrium constant vector:
    K_r = exp(-calc2_DgR(T)/R/T);
    
    % We now write the rate laws:
    NUM1 = f_j(1)*f_j(3)^2 - f_j(5)/K_r(1);
    DEN1 = (C(1) + C(2)*f_j(1) + C(3)*f_j(2) + C(4)*f_j(3))^2;
    %
    NUM2 = f_j(2)*f_j(3) - (f_j(1)*f_j(4))/K_r(2);
    DEN2 = C(5);
    
    % Complete rates with the dimensions [mol/g_cat/min]:
    f = zeros(NR,1);
    %
    f(1) = NUM1/DEN1;
    f(2) = NUM2/DEN2;
    
    % Effective dimensions for the rates:
    f = f*1000/60       ;%[mol/kg_cat/s]  

end

function f = fun_PFRbed(~,X)% V in m3, X in [-], fun in [-]

% This function evaluates the bed mass and energy balances.

    global NS MM nu cp_mix DhR
    global eps rho_cat
    global Mtot

    % We declare the value of the objective value:
    f = zeros(NS+1,1);
    
    % We assign the names to the vector X = [om; T];
    om_j = X(1:NS);
    T    = X(NS+1);%[K]
    
    % We express the reaction rates:
    Rk = calc_rates(om_j,T); %[mol/kg_cat/s]
       
    % We express the mass balances:
    f(1:NS) = 1/Mtot * nu'*Rk .* MM * rho_cat*(1-eps);
    
    % We express the energy balance:
    f(NS+1) = - 1/(Mtot*cp_mix) * DhR'*Rk * rho_cat*(1-eps);
   
end

function [value, isterminal, direction] = fun_event(~,X)% V in m3, X in [-] 

    % This function evaluates the stopping criteria for the integrator:
    % 1. value      = function that when zero, stops the solver
    % 2. isterminal = set to 1 to stop the solver at the event
    % 3. direction  = set to 0 if all zeroes are computed
    
    global Tlim
    
    value      = X(end) - Tlim;
    isterminal = 1;
    direction  = 0;
end




