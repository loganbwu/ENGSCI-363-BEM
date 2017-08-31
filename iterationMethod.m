%Iteration Method Aims to follow the Design Steps
%The end script will iterate through steps 1 to 6 until optimised
%parameters
%% Step 1 - Establish Design Specifications
v_U = 5.14444;                      % 10 knots
Omega = 14.6608;                    % Blade rotation/angular speed
Torque = 3.25;                      % Torque, but is this right
eta = 1;                            % Initial Guess for eta
P_E = Omega*Torque;                 % Power Extracted
P_S = P_E * eta;                    % Power Required 

%% Step 2 - Calculate Radius
R = sqrt(2*P_S/(C_P*eta*rho*pi*V_U^3));     %Calculate Radius of Blade 

%% Step 3 - Define Tip Speed Ratio
lambda = Omega*R/v_U;                       %Tip Speed Ratio not local ratio

%% Step 4 - Number of Blades
B = 6;

%% Step 5 - Select Airfoils (basically defines alpha at radius r)
% Airfoil Parameters
C_L = [0.0001	0.11	0.22	0.33	0.44	0.55	0.66	0.746	0.8247	0.8527	0.1325	0.1095	0.1533	0.203	0.2546	0.3082	0.362	0.42	0.477	0.5322	0.587];
C_D = [0.0103	0.0104	0.0108	0.0114	0.0124	0.014	0.0152	0.017	0.0185	0.0203	0.0188	0.076	0.134	0.152	0.171	0.19	0.21	0.231	0.252	0.274	0.297];
Alpha = linspace(0, 20, 21);
rho = 1.225;

% Construction Parameters
% Unsure if c - chord length and if it should be deterministic vs
% iterative?
c = [0.596	0.413	0.303	0.236	0.192	0.162	0.139	0.122	0.109	0.098];
beta = [0.79412481	0.455530935	0.286233997	0.188495559	0.127409035	0.085521133	0.055850536	0.033161256	0.013962634	0];
n_points = length(beta);

% Environment Parameters
% Defined in previous steps, v_U, Omega, B, R, lambda

%% Step 6)a) Initial Guess for a and a' and other initialisation
% Initialisation
a = ones(1,n_points)*0.33;          % axial induction factor, initial 1/3
adash = zeros(1,n_points);          % angular induction factor, intial 0
r = linspace(0.11, R, n_points);    % sample radii
lambda_r = r/R * lambda;            % local wind speed ratio
phi = 2/3*atan(1./lambda_r);        % wind angle
alpha = phi - beta;                 % angle of attack
solidity = B * c ./ (2*pi*r);       % (sigma-dash)

%% Step 6)b - c) 
%Calculate: wind angle, local angle of attack - beta, chord length - c,
%axial and angular induction factors and tip loss - F.
for i = n_points:-1:1     % start at tip
    for j = 1:10
        c_l = interp1(Alpha, C_L, alpha(i)*180/pi);
        c_d = interp1(Alpha, C_D, alpha(i)*180/pi);
        c_n = c_l * cos(phi(i)) + c_d * sin(phi(i));
        c_t = c_l * sin(phi(i)) - c_d * cos(phi(i));
        
        % tip loss factor
        f = B*(R-r(i))/(2*r(i)*sin(phi(i)));
        F = 2/pi*acos(exp(-f));
        
        %took f out for tut2
        
        % guess induction factors
        a(i) = 1./(4.*(sin(phi(i)).^2./(solidity(i).*c_n)) + 1);
        adash(i) = 1./(4.*sin(phi(i)).*cos(phi(i))./(solidity(i)*c_t)-1);
    end
    
    % use factors as guess for next segment (makes no difference but w/e)
    if i >= 2
        a(i-1) = a(i);
        adash(i-1) = adash(i);
    end
end

%% Step 6)d) Calculate p_T - tangential load at each cross-section
p_T = 1/2*rho * v_U^2*(1-a).^2 ./ sin(phi).^2; 
%% Step 6)e) Integrate numerically for incremental torque
Torque = B * trapz(r(~isnan(p_T)), r(~isnan(p_T)).*p_T(~isnan(p_T)));
%% Step 6)f) Power Calculations and Updates
P_E = Torque * Omega;           % extracted
P_T = 1/2*rho*pi*R^2*v_U^3;     % total
C_P = P_E/P_T;                  % power coefficient

