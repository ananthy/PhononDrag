clc; 
clearvars; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%General Parameters 

e = 1.60217657e-19;
k_B=1.3806504*10^-23; % J/K
m0=9.10938188*10^-31; % Free elctron mass in kg 
Df = 8.5*e; %Deformation potential constant  
mGaN = 0.22*m0; %Effective mass is very small, GaAs 
d = 100e-9; %Thickness of GaN layer 
amu = 1.6605e-27;
GaN_M = amu*83.7;
GaN_Dens = 6150;
vol_atom= GaN_M/(GaN_Dens);
ND = 1e13;  %Estimated Dislocation Density
ND2 = ND/3; 
GaN_c11 = 390e9; %Elastic constant 
GaN_c12 = 145e9; %Elastic constant 
GaN_c44 = 105e9; %Elastic constant 
h14 = 4.3e9; %Piezoelectric constant 
v_L  = sqrt(GaN_c11/GaN_Dens);
v_T1 = sqrt(GaN_c44/GaN_Dens);
v_T2 = sqrt ((GaN_c11-GaN_c12)/(2*GaN_Dens));
v_T = (v_T1 + v_T2)/2; 
v_s = 1/((1/3)*(1/v_L+1/v_T1+1/v_T2));
rho = 5.3*1e3; %Density, units of kg/m^3
h_bar=6.62606896*10^-34/(2*pi); % in Joule*seconds
Tdebye = 600; 
eps0 = 8.854e-12;
GaN_epss = 10.4*eps0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Wavefunction Parameters 
n2D = 1e17; %2DEG charge Density 
n = (11/32)*n2D; %Fang-Forward Exapression 
b = (12*mGaN*e^2*n/(GaN_epss*h_bar^2))^(1/3); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Fermi-Level 
EF = 0.107*e; %Where did we get this from? Lowermost Band-Bottom? 
qz_lim = k_B*Tdebye/(h_bar*v_s); %Debye limit for phonons (in GaN) 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Thermal Conductivity Constants 
wt_factor = 0.55; %along the direction of transport 
nu = GaN_c12/(GaN_c11+GaN_c12); 
gamma_mass_scatter = 0;  
b_s = 5.18e-10; 
a = 3.16e-10; 
b_e = a*sqrt(2)/3; 
gamma = 0.5; 
psi2 = 1; 
P = 1.82e-19;  
Cu = 132;
Temp = 50:25:325; 
for j = 1:length(Temp)
    T = Temp(j); 
    %%%%%%%%%%%%%%%%Phonon Scattering Rates 
    %tau_N_inv    =  (@(X)(k_B^2*T^2*X.^2)/(h_bar^2*A_const)); %Phonon N Process Scattering
    tau_N_inv    =  (@(X)(P*k_B^2*T^3*exp(-Cu/T)*X.^2)/(h_bar^2)); %Phonon N Process Scattering
    tau_DC_inv   =  (@(X) (wt_factor*ND*vol_atom^(4/3)*X.^(3)*k_B^3*T^3)/(h_bar^3*v_s^2)); %Core Dislocation scattering
    tau_S_inv    =  (@(X) (2^(3/2)/3^(7/2))*wt_factor*ND2*b_s^2*gamma^2*X*k_B*T/h_bar);  %Screw Dislocation scattering 
    tau_E_inv    =  (@(X) ((2^(3/2)/3^(7/2))*wt_factor*ND2*b_e^2*gamma^2*X*k_B*T/h_bar)*(0.5 +(1/24)*((1-2*nu)/(1-nu))^2*(1+sqrt(2)*(v_L/v_T)^2)^2)); %Edge Dislocation Scattering 
    tau_M_inv    =  (@(X) ((2^(3/2)/3^(7/2))*wt_factor*ND2*gamma^2*X*k_B*T/h_bar)*(b_s^2 + b_e^2*(0.5 +(1/24)*((1-2*nu)/(1-nu))^2*(1+sqrt(2)*(v_L/v_T)^2)^2))); %Mixed Dislocation Scattering 
    tau_UM_inv    = (@(X) ((vol_atom*k_B^4*T^4)/(4*pi*v_s^3*h_bar^4))*gamma_mass_scatter*X.^4); %Mass defect scattering                            
    tau_bound_inv =  v_s/(2.38*d); %Phonon boundary scattering 
    tau_C_inv    =  (@(X) tau_N_inv(X)+tau_UM_inv(X)+tau_bound_inv+psi2*(tau_DC_inv(X)+tau_S_inv(X)+tau_E_inv(X)+tau_M_inv(X))); %All scattering processes

  
    Sg_coeff = ((2*mGaN)^(3/2)*v_s*v_s)/(4*(2*pi)^3*k_B*T^2*e*rho*n2D); %Drag Coefficient      
    omegaq = (@(q,qz)v_s*sqrt(q.^2+qz.^2)); %Phonon frequency 
    Iqz_sq = (@(qz)abs(b^6./(b^2+qz.^2).^3)); %Electron-Phonon momentum conservation, in the z direction 
    kf = sqrt(2*mGaN*EF*h_bar^-2); %Fermi wave vector  
    Pi_0 = (@(E,q)(mGaN/(4*pi*h_bar^2*k_B*T))*((1-heaviside(q-2*kf).*(1-(2*kf./(q)).^2 ).^0.5)./(cosh((EF-E)/(2*k_B*T))).^2)); %Polarizability Function (Integrand) 
    Pi = (@(q)integral(@(E)Pi_0(E,q),1e-2*k_B*T,EF+15*k_B*T,'ArrayValued',true)); %Polarizability Function (Full), integral limits from 0 to 10*Fermi Energy 
    Form = (@(q)(8+9*(q/b)+3*(q/b).^2)./(8*(1+(q/b)).^3)); %Form-factor 
    Screen = (@(q)1+(2*pi*e^2.*Pi(q).*Form(q))./(GaN_epss*q)); %Screening function 
    FD_1 = (@(E)1./(1+exp((E-EF)./(k_B*T)))); %Fermi-Dirac function  
    FD_2 = (@(E,q,qz)1-1./(1+exp((E+h_bar*omegaq(q,qz)-EF)./(k_B*T))));
    Eq = (@(q)h_bar^2*q.^2/(2*mGaN)); 
    G_coeff = (@(q,qz)((1-exp(-h_bar*omegaq(q,qz)/(k_B*T)))./(h_bar*omegaq(q,qz)))); %Energy integral coefficient 
    G = (@(q,qz)(integral(@(E)FD_1(E).*FD_2(E,q,qz)./(sqrt(E-((h_bar*omegaq(q,qz)-Eq(q)).^2./(4*Eq(q))))),1e-2*k_B*T,EF+15*k_B*T,'ArrayValued',true))); %Energy integral      
    %integrand2 = (@(q,qz) (Df^2*Iqz_sq(qz).*q.^2.*(q.^2+qz.^2).*G(q,qz))./((sinh(h_bar*omegaq(q,qz)/(2*k_B*T))).^2.*(Screen(q)).^2.*tau_ph_inv(q,qz))); 
    tau_ph_inv =(@(q,qz)tau_C_inv(h_bar*omegaq(q,qz)/(k_B*T))); %Phonon scattering rate 
    piez_pot = (@(q,qz)(e*h14)^2*(8*qz.^2.*q.^2+q.^4)./(2*(q.^2+qz.^2).^3));%Piezoelectric component of electron-phonon scattering 
    integrand2 =  (@(q,qz) ((piez_pot(q,qz)+Df^2).*Iqz_sq(qz).*q.^2.*(q.^2+qz.^2).*G_coeff(q,qz).*G(q,qz))./((sinh(h_bar*omegaq(q,qz)/(2*k_B*T))).^2.*(Screen(q)).^2.*tau_ph_inv(q,qz))); %Phonon drag integrand 
    Drag(j) = real(2*Sg_coeff*quad2d(@(q,qz)integrand2(q,qz),0,qz_lim,0,qz_lim,'AbsTol',1e-8,'MaxFunEvals',80000)); 
end

plot(Temp,real(Drag),'--r'); 

