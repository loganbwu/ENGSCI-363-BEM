% Airfoil Parameters
C_L = [0.0001	0.11	0.22	0.33	0.44	0.55	0.66	0.746	0.8247	0.8527	0.1325	0.1095	0.1533	0.203	0.2546	0.3082	0.362	0.42	0.477	0.5322	0.587];
C_D = [0.0103	0.0104	0.0108	0.0114	0.0124	0.014	0.0152	0.017	0.0185	0.0203	0.0188	0.076	0.134	0.152	0.171	0.19	0.21	0.231	0.252	0.274	0.297];
Alpha = linspace(0, 20, 21);
rho = 1.225;

% Construction Parameters
c = [0.596	0.413	0.303	0.236	0.192	0.162	0.139	0.122	0.109	0.098];
beta = [0.79412481	0.455530935	0.286233997	0.188495559	0.127409035	0.085521133	0.055850536	0.033161256	0.013962634	0];
npoints = length(beta);

% Environment Parameters
B = 4;
R = 1.1;
Omega = 23.3838;
v_U = 5.14444;                      % 10 knots
lambda = Omega*R/v_U;               % = 5

% Initialisation
r = linspace(0.1, R, npoints);     % sample radii
lambda_r = Omega*r/v_U;            % local wind speed ratio

a = [nan(1,npoints-1) 0.33];          % initial guess
adash = [nan(1,npoints-1) 0];          % initial guess
phi = nan(1,npoints);        % wind angle
alpha = nan(1,npoints);                 % angle of attack
solidity = nan(1,npoints);       % (sigma-dash)

for i = npoints:-1:1     % start at tip
    for j = 1:10
        % get angles
        phi(i) = atan((1-a(i))/(lambda_r(i)*(1+adash(i))));
        phi(i) = 2/3*atan(1/lambda_r(i));
        alpha(i) = phi(i)-beta(i);
        
        % get lift/drag
        c_l = interp1(Alpha, C_L, alpha(i)*180/pi);
        c_d = interp1(Alpha, C_D, alpha(i)*180/pi);
        c_n = c_l * cos(phi(i)) + c_d * sin(phi(i));
        c_t = c_l * sin(phi(i)) - c_d * cos(phi(i));
        
        c(i) = 8*pi*r(i)/(B*c_l) * (1-cos(phi(i)));
        solidity(i) = B*c(i)/(2*pi*r(i));
        
        % tip loss factor
        f = B*(R-r(i))/(2*r(i)*sin(phi(i)));
        F = 2/pi*acos(exp(-f));
        F=1;

        % induction factors
        a(i) = 1/(4*F*(sin(phi(i))^2/(solidity(i)*c_n)) + 1);
        adash(i) = 1/(4*F*sin(phi(i))*cos(phi(i))/(solidity(i)*c_t)-1);
    end
    
    % use factors as guess for next segment (makes no difference but w/e)
    if i >= 2
        a(i-1) = a(i);
        adash(i-1) = adash(i);
    end
end

p_T = 1/2*rho * v_U^2*(1-a).^2 ./ sin(phi).^2;  % tangential force/length

% integrate forces to find torque
% Q = B * int r.*p_T.*dr
torque = B * trapz(r(~isnan(p_T)), r(~isnan(p_T)).*p_T(~isnan(p_T)));

P_E = torque * Omega;           % extracted
P_T = 1/2*rho*pi*R^2*v_U^3;     % total
P_C = P_E/P_T                  % power coefficient

% plot(r/R, a, r/R, adash); legend('a', 'a''')
% xlabel('r/R');
% ylabel('Induction Factor')