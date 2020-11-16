% ICO - ese05

% Written Exam Test (reference 11-07-17)

clc
close all
clearvars

global R NS NR nu Dh0f Cp Prif Trif Pval
global Tfed yfed
global cnv Q

%% Data Declaration

% We have the universal gas constant:
R  = 8.3144621/4.184              ;%[cal/(mol*K)]

% We list the number of species in the reactor:
NS = 5;

% We build the atom-element matrix:
A = [2    0    0    4    2      %H
     0    1    1    1    0      %C
     0    1    2    0    1];    %O
     %H2  %CO  %CO2 %CH4 %H2O
     
% We check the rank and the number of reactions:
NR = NS - rank(A);

% We have verified NR = 2, so we use the SAB and WGS reactions:
nu = [-4    0   -1    1    2      %SAB
       1   -1    1    0   -1];    %WGS
       %H2  %CO  %CO2 %CH4 %H2O
       
% We have the formation enthalpies:
Dh0f = [0; -26420; -94050; -17890; -57800]  ;%[cal/mol]

% We have the specific heats:
Cp = [2207; 2253; 3284; 3438; 2662]        ;%[cal/mol/K]       

% We have the reference state:
Prif = 101325   ;%[Pa]
Trif = 298      ;%[K]

% The feed temperature is known:
Tfed = 500      ;%[K]

% The feed composition is known:
yfed = [4/5; 0.0; 1/5; 0.0; 0.0];%[-]

% The operating pressure is known:
Pval = 8*101325 ;%[Pa]

% The CO2 total conversion is known:
cnv = 0.95      ;%[-]

% The heat loss rate is known:
Q = 7000        ;%[cal/mol]

%% Solution

% We have the initial conditions:
X0 = [-0.01; 600];

% We solve the problem:
Xf = fsolve(@fun_r101,X0);

% Solutions in terms of lam and Tout:
lam = zeros(NR,1);
%
lam(1) = Xf(1) + yfed(3)*cnv;
lam(2) = Xf(1);
Tout   = Xf(2);

% Outlet molar fractions:
yout = (yfed + nu'*lam)/(1 + sum(nu'*lam));

% Calculation of Tin:
DhR = calc_DhR(Tout);
%
Tin = Tout + (sum(lam.*DhR) + Q)/sum(yfed.*Cp);

% Calculation of Tsc:
Tsc = Tout - (Tin-Tfed)*sum(yfed.*Cp)/sum(yout.*Cp);

%% Results

fprintf('Exam Results\n');
fprintf('\nTemperatures:\n');
fprintf('\nTin  = %.3f K\nTout = %.3f K\nTsc  = %.3f K\n',Tin,Tout,Tsc);
fprintf('\nOutlet molar fractions from reactor (H2 CO CO2 CH4 H2O):\n');
fprintf('%.5f\n',yout);
fprintf('\n');

% Function Declaration:

function f = calc_Keq(T)% T in [K], f in [-]

    global R NR
    
    f = zeros(NR,1);
    
    % SAB
    f(1) = exp(1/R*(56000/T^2 + 34633/T - 16.4*log(T) + 0.00557*T) + 33.165);
    % WGS
    f(2) = exp(-(-8514 + 7.71*T)/R/T);

end

function f = calc_DhR(T)% T in [K], f in [cal/mol]

    global nu Dh0f Cp Trif
    
    f = nu*Dh0f + nu*Cp*(T - Trif);
      
end

function f = fun_r101(x)

    global NR nu Prif Pval
    global yfed
    global cnv
    
    % We define the degree of advance:
    lam = zeros(NR,1);
    
    % We define the problem's unknowns:
    lam(1) = x(1) + yfed(3)*cnv;
    lam(2) = x(1);
    Tout = x(2);
    
    % We define the problem's equations:
    f = zeros(2,1);
    
    % We have the reacted molar fraction vector:
    yrc = nu'*lam;
    
    % We have the equilibrium molar fractions:
    yeq = (yfed + yrc)/(1 + sum(yrc));
    
    % We write the activities:
    aeq = yeq*Pval/Prif;
    
    % We have the 2 equilibrium conditions:
    Keq = calc_Keq(Tout);
    
    for j=1:NR
        nuvec = nu(j,:)';
        f(j) = Keq(j)/prod(aeq.^nuvec) - 1;
    end
    
end

