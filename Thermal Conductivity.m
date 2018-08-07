clc; 
clearvars;
% Basic Constants 
eps0 = 8.854e-12;
hbar = 1.05457173e-34;
e = 1.60217657e-19;
m0 = 9.10938291e-31;
x = 0; % Alloy composition  (InxGayAl(1-x-y)N)
y = 1;  
xdash = x;
amu = 1.6605e-27;
kb = 1.38e-23;
% Material Constants
% GaN
%=============
GaN_epss = 9.5*eps0;
GaN_epsinf = 5.35*eps0;
GaN_mstar = 0.27*m0;
GaN_alpha = 0.189/e;
GaN_hwlo = 91.2*1e-3*e;
GaN_Tdebye = 600;
GaN_gamma = 0.5;
GaN_a = 3.16e-10;
GaN_c = 5.18e-10;
GaN_c11 = 390e9;
GaN_c12 = 145e9;
GaN_c44 = 105e9;
GaN_Dac = 8.3 *e;
GaN_Uab = 1.2*e;
GaN_delta = GaN_a/2;
GaN_M = amu*83.7;
GaN_Dens = 6150;
GaN_vol_atom= GaN_M/(GaN_Dens);
GaN_ion_r = GaN_vol_atom^(1/3);
%===============
% InN
%===============
InN_epss = 15.3*eps0;
InN_epsinf = 8.4*eps0;
InN_mstar = 0.115*m0;
InN_alpha = 1.43/e;
InN_hwlo = 89.0*1e-3*e;
InN_Tdebye = 660;
InN_gamma = 0.3;
InN_a = 3.54e-10;
InN_c = 5.76e-10;
InN_c11 = 223e9;
InN_c12 = 115e9;
InN_c44 = 48e9;
InN_Dac = 7.1*e;
InN_Uab = 1.8*e;
InN_delta = InN_a/2;
InN_M = amu*128.8;
InN_Dens = 6810;
InN_vol_atom= InN_M/(InN_Dens);
InN_ion_r = InN_vol_atom^(1/3);
%===============
% AlN
%===============
AlN_epss = 9.14*eps0;
AlN_epsinf = 4.6*eps0;
AlN_mstar = 0.35*m0;
AlN_alpha = 0.044/e;
AlN_hwlo = 99.2*1e-3*e;
AlN_Tdebye = 1150;
AlN_gamma = 0.5;
AlN_a = 3.11e-10;
AlN_c = 4.98e-10;
AlN_c11 = 410e9;
AlN_c12 = 149e9;
AlN_c44 = 125e9;
AlN_Dac = 9.5*e;
AlN_Uab = 2.8*e;
AlN_delta = AlN_a/2;
AlN_M = 40.99*amu;
AlN_Dens = 3266;
AlN_vol_atom= AlN_M/(AlN_Dens);
AlN_ion_r = AlN_vol_atom^(1/3); 
%===============
% Overall Properties
%===============
epss = x*InN_epss+y*GaN_epss+(1-x-y)*AlN_epss;
epsinf = x*InN_epsinf+y*GaN_epsinf+(1-x-y)*AlN_epsinf;
mstar = x*InN_mstar+y*GaN_mstar+(1-x-y)*AlN_mstar;
alpha = x*InN_alpha+y*GaN_alpha+(1-x-y)*AlN_alpha;
hwlo = x*InN_hwlo+y*GaN_hwlo+(1-x-y)*AlN_hwlo; %Optical Phonon Energy 
c = x*InN_c+y*GaN_c+(1-x-y)*AlN_c;
ion_r = x*InN_ion_r+y*GaN_ion_r+(1-x-y)*AlN_ion_r;
Dac = x*InN_Dac+y*GaN_Dac+(1-x-y)*AlN_Dac;
Uab = x*InN_Uab+y*GaN_Uab+(1-x-y)*AlN_Uab;
delta =  x*InN_delta+y*GaN_delta+(1-x-y)*AlN_delta;
M = x*InN_M+y*GaN_M+(1-x-y)*AlN_M;
c11 = (GaN_c11*(GaN_delta)^4)/(delta)^4;
c12 = (GaN_c12*(GaN_delta)^4)/(delta)^4;
c44 = (GaN_c44*(GaN_delta)^4)/(delta)^4;
cL = c11 + (2/5)*(c12 + 2*c44 - c11);
cT = c44 - (1/5)*(c12 + 2*c44 - c11); 
Tdebye = GaN_Tdebye*(GaN_M)^0.5*(GaN_delta)^1.5/((M)^0.5*(delta)^1.5);
Dens = x*InN_Dens+y*GaN_Dens+(1-x-y)*AlN_Dens;
v_L  = sqrt(c11/Dens);
v_T1 = sqrt(c44/Dens);
v_T2 = sqrt ((c11-c12)/(2*Dens));
v = 1/((1/3)*(1/v_L+1/v_T1+1/v_T2));
vol_atom = x*InN_vol_atom+y*GaN_vol_atom+(1-x-y)*AlN_vol_atom;
%=================================================================

%=================================================================
d = 1e-6; %Film thickness   
%=================================================================


%Dislocation Parameters and Alloy Scattering 
%======================================================================
ND = 1e13; %Dislocation density 
ND2 = ND/3;
wt_factor = 0.55; %along the direction of transport 
nu = c12/(c11+c12); 
v_T = (v_T1 + v_T2)/2; 
b_s = 5.18e-10; 
a = 3.16e-10; 
b_e = a*sqrt(2)/3; 
gamma = 0.5; 
psi = 1; 
gamma_mass_scatter = y*(((GaN_M-M)/(M))^2 + 39*((GaN_ion_r-ion_r)/(ion_r))^2)+(1-x-y)*( ((AlN_M-M)/(M))^2 + 39*((AlN_ion_r-ion_r)/(ion_r))^2);
% ======================================================================
Temp = 25:2:400;
%For GaN
P = 1.82e-19;  
Cu = 132;

% For AlN
% P = 3.3e-19;  
% Cu = 382;
for i = 1:length(Temp)
    tau_N_inv    =  (@(X)(P*kb^2*Temp(i)^3*exp(-Cu/Temp(i))*X.^2)/(hbar^2)); %Phonon N Process Scattering
    tau_DC_inv   =  (@(X) (wt_factor*ND*vol_atom^(4/3)*X.^(3)*kb^3*Temp(i)^3)/(hbar^3*v^2)); %Core Dislocation scattering
    tau_S_inv    =  (@(X) (2^(3/2)/3^(7/2))*wt_factor*ND2*b_s^2*gamma^2*X*kb*Temp(i)/hbar);  %Screw Dislocation scattering 
    tau_E_inv    =  (@(X) ((2^(3/2)/3^(7/2))*wt_factor*ND2*b_e^2*gamma^2*X*kb*Temp(i)/hbar)*(0.5 +(1/24)*((1-2*nu)/(1-nu))^2*(1+sqrt(2)*(v_L/v_T)^2)^2)); %Edge Dislocation Scattering 
    tau_M_inv    =  (@(X) ((2^(3/2)/3^(7/2))*wt_factor*ND2*gamma^2*X*kb*Temp(i)/hbar)*(b_s^2 + b_e^2*(0.5 +(1/24)*((1-2*nu)/(1-nu))^2*(1+sqrt(2)*(v_L/v_T)^2)^2))); %Mixed Dislocation Scattering 
    tau_UM_inv    = (@(X) ((vol_atom*kb^4*Temp(i)^4)/(4*pi*v^3*hbar^4))*gamma_mass_scatter*X.^4); %Mass defect scattering                            
    tau_bound_inv =  v/(2.38*d); %Boundary scattering 
    tau_C_inv    =  (@(X) tau_N_inv(X)+tau_UM_inv(X)+tau_bound_inv+psi*(tau_DC_inv(X)+tau_S_inv(X)+tau_E_inv(X)+tau_M_inv(X))); %All scattering processes
    Coeff = kb^4*(Temp(i)^3)/(2*pi^2*hbar^3*v);
    integrand = (@(X)(X.^4.*exp(X))./((exp(X)-1).^2.*tau_C_inv(X)));  
    k1(i) = Coeff*integral(integrand,0,Tdebye/Temp(i)); %Thermal conductivity 
end

plot(Temp,k1,'r');
hold on;


