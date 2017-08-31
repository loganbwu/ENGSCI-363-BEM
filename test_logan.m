% Logan's test script

% known
v_U = 5.14444;  % 10 knots
P_S = 50;
rho = 1.225;
B=6;
hub_radius = 0.1;

% define aerofoil
C_L_vec = [0.0001	0.11	0.22	0.33	0.44	0.55	0.66	0.746	0.8247	0.8527	0.1325	0.1095	0.1533	0.203	0.2546	0.3082	0.362	0.42	0.477	0.5322	0.587];
C_D_vec = [0.0103	0.0104	0.0108	0.0114	0.0124	0.014	0.0152	0.017	0.0185	0.0203	0.0188	0.076	0.134	0.152	0.171	0.19	0.21	0.231	0.252	0.274	0.297];
alpha_vec = linspace(0, 20, 21);

% guesses
C_P = 0.3;
eta = 1;

R = sqrt(2*P_S/(C_P*eta*rho*pi*v_U^3));
lambda = Omega*R/v_U;

% define blade elements
n_elements = 10;
r =      linspace(r,R,n_elements);
a =     [nan(1,n_elements-1) 0.33];
adash = [nan(1,n_elements-1) 0   ];
alpha =  nan(1,npoints);
beta =   zeros(1,npoints);
% solidity = nan(1,npoints);

% initialise iterative guesses
for i = n_elements:-1:1
    lambda_r = Omega*r(i)/v_U;
    phi = atan(2/(3*lambda_r));
    c = 8*pi*r/(B*C_L) * sin(phi)/(3*lambda_r);
end
