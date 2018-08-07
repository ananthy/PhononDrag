clc; 
clearvars; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%General Parameters 

e = 1.60217657e-19;
k_B=1.3806504*10^-23; % J/K
m0=9.10938188*10^-31; % Free elctron mass in kg
Df = 8.5*e; %Deformation potential constant  
mGaN = 0.22*m0; %Effective mass is very small, GaAs 
GaN_Dens = 6150; %GaN Density 
GaN_c11 = 390e9; %Elastic constant 
GaN_c12 = 145e9; %Elastic constant 
GaN_c44 = 105e9; %Elastic constant 
v_L  = sqrt(GaN_c11/GaN_Dens); %Longitudinal velocity 
v_T1 = sqrt(GaN_c44/GaN_Dens); %Transverse velocity 
v_T2 = sqrt ((GaN_c11-GaN_c12)/(2*GaN_Dens)); %Transverse velocity 
v_s = 1/((1/3)*(1/v_L+1/v_T1+1/v_T2)); %Average speed 
h_bar=6.62606896*10^-34/(2*pi); % in Joule*seconds
h14 = 4.3e9; %Piezoelectric coefficient 
Tdebye = 600; %GaN Debye temperature 
cL = 2.65e11; %Elastic constant 
eps0 = 8.854e-12; %Dielectric constants 
GaN_epss = 10.4*eps0;
GaN_epsinf = 5.47*eps0; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Wavefunction Parameters 
n2D = 1e17; %2DEG sheet density 
n = (11/32)*n2D; %Fang-Forward expression for the wavefunction  
b = (12*mGaN*e^2*n/(GaN_epss*h_bar^2))^(1/3);
g2D = mGaN/(pi*h_bar^2); %2D Density of states 
% plot(z/1e-9,psi); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Fermi-Level 
EF = (n2D/g2D); %Fermi-level 
hwlo = 91.2*1e-3*e; %Optical phonon energy 

%let's look at the coefficient for optical phonon scattering 
coeff_opt = (e^2*mGaN*hwlo*(1/GaN_epsinf-1/GaN_epss))/(8*pi^2*h_bar^3); %Optical Phonon scattering rate 
qz_lim = k_B*Tdebye/(h_bar*v_s); %Debye limit for phonons (in GaN) 

% coefficients for roughness scattering 
Delta = 1e-9; %Average Displacement of the interface (set for thick GaN)  
Lambda = 7.5e-9; %Auto-Correlation Length 
Field = e^2*n2D/(2*GaN_epss);
coeff_rough = (mGaN*Delta^2*Lambda^2*Field^2)/(2*h_bar^3*pi); %roughness scattering coefficient 
H = (@(E,Th)(8+9*((2*sqrt(2*mGaN*E*h_bar^-2).*sin(Th/2))/b)+3*((2*sqrt(2*mGaN*E*h_bar^-2).*sin(Th/2))/b).^2)./(8*(1+((2*sqrt(2*mGaN*E*h_bar^-2).*sin(Th/2))/b)).^3));
Iq_plus_minus = @(qpar,qz)abs(b^6./(b^2+qz.^2).^3)./(qpar.^2+qz.^2);
Iq_plus = (@(qparv)arrayfun(@(qpar)integral(@(qz)Iq_plus_minus(qpar,qz),0,qz_lim/10),qparv)); 


Temp = 50:25:325; 
S_diff = zeros(1,length(Temp)); 
for j = 1:length(Temp)
    T = Temp(j); 
    coeff_piez = ((e*h14)^2*mGaN*k_B*T)/(4*pi*h_bar^3); %Piezoelectric Scattering 
    coeff = (3*b*mGaN*Df^2*k_B*T)/(32*pi*h_bar^3*cL); %Deformation Potential Scattering 
    kf = sqrt(2*mGaN*EF*h_bar^-2); %Fermi wavevector  
    f = (@(E,Th)(mGaN/(4*pi*h_bar^2*k_B*T))*((1-heaviside((2*sqrt(2*mGaN*E*h_bar^-2).*sin(Th/2))-2*kf).*(1-(2*kf./(2*sqrt(2*mGaN*E*h_bar^-2).*sin(Th/2))).^2 ).^0.5)./(cosh((EF-E)/(2*k_B*T))).^2)); %Polarizability Function (Integrand) 
    Pol = (@(T)integral(@(E)f(E,T),0,10*EF,'ArrayValued',true)); %Polarizability Function 
    S = (@(E,T)1+(e^2.*H(E,T).*Pol(T))./(2*GaN_epss*(2*sqrt(2*mGaN*E*h_bar^-2).*sin(T/2)))); %Screening function 
    integrand = (@(E,T)(1./(S(E,T)).^2).*(1-cos(T)));%Integrand for deformation potential scattering  
    FDir = (@(E)1./(1+exp((E-EF)/(k_B*T)))); %Fermi-Dirac function  
    DerFDir = (@(E)(exp((E-EF)/(k_B*T))./(1+exp((E-EF)/(k_B*T))).^2)*(1/(k_B*T))); %Derivative of Fermi-Dirac function 
    FDirDenom = (@(E)1./(1-FDir(E))); 
    Nq = 1/(exp(hwlo/(k_B*T))-1); %Number of optical phonons 
    integrand_opt = (@(E,T)FDirDenom(E).*(1-cos(T)).*((1-FDir(E+hwlo))*Nq.*Iq_plus((2*sqrt(2*mGaN*E*h_bar^-2).*sin(T/2))) +  (1-FDir(E-hwlo))*(Nq+1).*heaviside(E-hwlo).*Iq_plus((2*sqrt(2*mGaN*E*h_bar^-2).*sin(T/2))))); %Integrand for Optical Phonon Scattering 
    integrand_rough = (@(E,T)(1-cos(T)).*(1./(S(E,T))).^2.*exp(-Lambda.^2.*(2*sqrt(2*mGaN*E*h_bar^-2).*sin(T/2)).^2/4)); %Integrand for roughness scattering 
    E_low  = 1e-5*k_B*T; %Energy limits   
    E_high = EF+10*k_B*T; 
    
    angle_low = pi/1e5; %Angle limits 
    angle_high = 2*pi;
   
    ScatrateDef(j) = coeff*integral2(@(E,T)integrand(E,T).*DerFDir(E).*E,E_low,E_high,angle_low,angle_high)/integral(@(E)DerFDir(E).*E,E_low,E_high); %Energy, angle averaged deformation potential scattering rate 
    ScatrateOpt(j) = coeff_opt*integral2(@(E,T)integrand_opt(E,T).*DerFDir(E).*E,E_low,E_high,angle_low,angle_high)/integral(@(E)DerFDir(E).*E,E_low,E_high);  %Energy, angle averaged optical phonon scattering rate 
    ScatrateRough(j) = coeff_rough*integral2(@(E,T)integrand_rough(E,T).*DerFDir(E).*E,E_low,E_high,angle_low,angle_high)/integral(@(E)DerFDir(E).*E,E_low,E_high); %Energy, angle averaged roughness scattering rate 
    Scatrate(j) = ScatrateDef(j)+ScatrateOpt(j)+ScatrateRough(j); %Total scattering rate
   
    
    E_vec = linspace(E_low,E_high,101);
    S_num = zeros(1,length(E_vec));
    S_denom = zeros(1,length(E_vec)); 
    
    for u=1:1:length(E_vec) %Energy 
        ScatrateDefE =   coeff*integral(@(T)integrand(E_vec(u),T),angle_low,angle_high,'ArrayValued',true);
        ScatrateOptE =   coeff_opt*integral(@(T)integrand_opt(E_vec(u),T),angle_low,angle_high,'ArrayValued',true);
        ScatrateRoughE = coeff_rough*integral(@(T)integrand_rough(E_vec(u),T),angle_low,angle_high,'ArrayValued',true); 
        ScatrateE = ScatrateDefE+ ScatrateRoughE + ScatrateOptE;         
        DerFDir1 = ((exp((E_vec(u)-EF)/(k_B*T))/(1+exp((E_vec(u)-EF)/(k_B*T)))^2)*(1/(k_B*T)));         
        S_num(u) = E_vec(u)*DerFDir1*(E_vec(u)-EF)/ScatrateE; %Numerator of Seebeck coefficient  
        S_denom(u) = E_vec(u)*DerFDir1/ScatrateE; %Denominator of Seebeck coefficient                    
    end
    
    S_diff(j) = (trapz(E_vec,S_num))/(trapz(E_vec,S_denom))*(1/(e*T)); %Diffusive Seebeck coefifient 
    
end

%Plots of mobility 
mu_def = e./(mGaN*ScatrateDef); 
mu_Opt = e./(mGaN*ScatrateOpt); 
mu_Rough = e./(mGaN*ScatrateRough); 
mu = e./(Scatrate*mGaN);

%Plot of diffusive Seebeck coefficient 
plot(Temp,S_diff,'--r'); 
