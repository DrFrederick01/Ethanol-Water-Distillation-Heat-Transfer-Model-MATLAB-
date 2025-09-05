% This creates a table for mixture-side properties
mixtureData = table( ...
    [15; 60; 80; 95.6], ...                             % Weight percentage (%)
    [88.413;  80.574; 78.812; 77.955], ...            % T_mix [kg/m^3]
    [935.677; 827.946; 775.789; 747.046], ...              % rho_l_mix [kg/m^3]
    [0.653; 0.964; 1.226; 1.523], ...                    % rho_v_mix [kg/m^3]
    [2080.934e3; 1432.974e3; 1135.451e3; 911.405e3], ...      % h_fg_mix [J/kg]
    [0.662; 0.778; 0.453; 0.217], ...                    % k_l_mix [W/(m·K)]
    [510.665e-6; 418.556e-6; 353.306e-6; 416.551e-6], ...        % mu_l_mix [Pa·s]
    [2.681;  2.094;  2.722;  5.905], ...                    % Pr_l (dimensionless)
    [34.034e-3; 19.913e-3; 17.980e-3;  16.962e-3], ...         % sigma_mix [N/m]
     'VariableNames', {'wtPercent', 'T_mix','rho_l_mix','rho_v_mix','h_fg_mix',...
                      'k_l_mix','mu_l_mix','Pr_l','sigma_mix'});
                  
% Display the table
disp(mixtureData);


% Create a struct array for water-side properties at two temperatures.

% First element: T = 150 °C
waterData(1).T_steam       = 150;               % Temperature in °C
waterData(1).rho_l     = 917.008;           % [kg/m^3]
waterData(1).rho_v     = 2.548;             % [kg/m^3]
waterData(1).h_fg      = 2113.746e3;        % [J/kg] (converted from 2113.746 kJ/kg)
waterData(1).k_l       = 0.681;             % [W/(m·K)]
waterData(1).mu_l      = 182.611e-6;        % [Pa·s] (converted from 182.611 µPa·s)
waterData(1).sigma     = 48.646e-3;         % [N/m] (converted from 48.646 mN/m)

% Second element: T = 180 °C
waterData(2).T_steam       = 180;               % Temperature in °C
waterData(2).rho_l     = 888.999;           % [kg/m^3]
waterData(2).rho_v     = 5.159;             % [kg/m^3]
waterData(2).h_fg      = 2014.161e3;        % [J/kg] (converted from 2014.161 kJ/kg)
waterData(2).k_l       = 0.677;             % [W/(m·K)]
waterData(2).mu_l      = 150.381e-6;        % [Pa·s] (converted from 150.381 µPa·s)
waterData(2).sigma     = 42.037e-3;         % [N/m] (converted from 42.037 mN/m)

% Third element: T = 210 °C
waterData(3).T_steam = 210;          % [°C]
waterData(3).rho_l   = 852.718;      % [kg/m^3]
waterData(3).rho_v   = 9.588;        % [kg/m^3]
waterData(3).h_fg    = 1899.639e3;   % [J/kg] (converted from 1899.639 kJ/kg)
waterData(3).k_l     = 0.653;        % [W/(m·K)]
waterData(3).mu_l    = 127.865e-6;   % [Pa·s] (converted from 127.856 µPa·s)
waterData(3).sigma   = 35.186e-3;    % [N/m] (converted from 35.186 mN/m)


% Fourth element: T = 220 °C
waterData(4).T_steam = 220;          % [°C]
waterData(4).rho_l   = 840.219;      % [kg/m^3]
waterData(4).rho_v   = 11.615;        % [kg/m^3]
waterData(4).h_fg    = 1857.368e3;   % [J/kg] (converted from 1857.368 kJ/kg)
waterData(4).k_l     = 0.645;        % [W/(m·K)]
waterData(4).mu_l    = 121.768e-6;   % [Pa·s] (converted from 121.768 µPa·s)
waterData(4).sigma   = 32.863e-3;    % [N/m] (converted from 32.683 mN/m)




%% --- Global Constants and Other Parameters ---
g            = 9.81;        % [m/s^2]
theta        = 0.0;         % [degrees]
delta        = 1e-3;        % [m] thickness of conduction through the plate
K_conductivity = 16.7;      % [W/(m·K)] plate conductivity
L            = 1;           % [m] characteristic length (water side)
tol          = 1e-4;        % convergence tolerance
maxIter      = 100;         % maximum number of iterations





% compute_Ra_water  Computes the water-side Rayleigh-type number.

function Ra_water = compute_Ra_water(g, theta, rho_l_water, rho_v_water, h_fg_water, L, k_l_water, mu_l_water, DeltaT_sub, sigma_water)

term1 = (g * cos(theta) * rho_l_water * (rho_l_water - rho_v_water) * h_fg_water * L) ...
        / (k_l_water * mu_l_water * DeltaT_sub);
term2 = (sigma_water / (g * (rho_l_water - rho_v_water) * cos(theta)))^(3/2);
Ra_water = term1 * term2;
end


%-----------------------------------------------------------------------
%compute nusslet using specific range
function Nu = compute_Nu(Ra)

    if (Ra >= 1.0e6) && (Ra < 1.0e8)
        % Use Nu = 0.69 * Ra^0.20
        Nu = 0.69 * Ra^0.20;
    elseif (Ra >= 1.0e8) && (Ra < 1.0e10)
        % Use Nu = 0.81 * Ra^0.193
        Nu = 0.81 * Ra^0.193;
    else
        % If Ra is outside the given ranges, handle it here:
        % Option 1: Return an error or warning
        warning('Rayleigh number (%.4e) is outside the known range (1e6 to 1e10).', Ra);
        Nu = NaN;
    end
end

%-----------------------------------------------------------------------------
%calculate h water 
function h_water = compute_h_water(Nu, k_l_water, sigma_water, g, rho_l_water, rho_v_water, theta)
% compute_h_water  Computes the water-side heat transfer coefficient.
%
% Uses:
%   Nu = (h_water/k_l_water) * sqrt( sigma_water/(g*(rho_l_water - rho_v_water)*cos(theta)) )
%
% Rearranged for h_water:
%   h_water = Nu * k_l_water / sqrt( sigma_water/(g*(rho_l_water - rho_v_water)*cos(theta)) )
%

sqrtTerm = sqrt(sigma_water/(g*(rho_l_water - rho_v_water)*cos(theta)));
h_water = Nu * k_l_water / sqrtTerm;
end

function q = compute_heat_flux(h_water,DeltaT_sub)
q = h_water*DeltaT_sub;
end


%----------------------------------------------------------------------------------------------
%computing the kinematic viscosity V variable
function nu = compute_kinViscosity(mu, rho)

%   nu = mu / rho
nu = mu/rho;
end

%----------------------------------------------------------------------------------------------
function la = compute_la_mix(sigma_mix, g, rho_l_mix, rho_v_mix)
% compute_la_mix.

la = sqrt(sigma_mix / ( g*(rho_l_mix - rho_v_mix) ));
end

%----------------------------------------------------------------------------------------------
function h_mix = compute_h_mix(k_l_mix, la_mix, Pr_l, q_dblprime, rho_v_mix, h_fg_mix, nu_mix, pressure, sigma_mix)
% compute_h_mix  Computes the mixture-side heat transfer coefficient.

h_mix = (k_l_mix / la_mix) * 7.0e-4 * (Pr_l^0.35) ...
    * ((q_dblprime * la_mix)/(rho_v_mix * h_fg_mix * nu_mix))^0.7 ...
    * ((pressure * la_mix / sigma_mix))^0.7;
end






%% Set desired temperature for interpolation
% In this case, we want the water properties at T_desired = 217 °C.
T_steam = 217;

% Known temperatures from waterData:
T_low  = waterData(3).T_steam;   % 210 °C
T_high = waterData(4).T_steam;    % 220 °C

% Compute linear interpolation factor
factor = (T_steam - T_low) / (T_high - T_low);

% % Interpolate each property and assign to water-side variables
% rho_l_water = waterData(3).rho_l + factor * (waterData(4).rho_l - waterData(3).rho_l);
% rho_v_water = waterData(3).rho_v + factor * (waterData(4).rho_v - waterData(3).rho_v);
% h_fg_water  = waterData(3).h_fg  + factor * (waterData(4).h_fg  - waterData(3).h_fg);
% k_l_water   = waterData(3).k_l   + factor * (waterData(4).k_l   - waterData(3).k_l);
% mu_l_water  = waterData(3).mu_l  + factor * (waterData(4).mu_l  - waterData(3).mu_l);
% sigma_water = waterData(3).sigma + factor * (waterData(4).sigma - waterData(3).sigma);



% Choose water condition index (here, waterData(1) for 150 °C)
T_steam      = waterData(1).T_steam;
rho_l_water  = waterData(1).rho_l;
rho_v_water  = waterData(1).rho_v;
h_fg_water   = waterData(1).h_fg;
k_l_water    = waterData(1).k_l;
mu_l_water   = waterData(1).mu_l;
sigma_water  = waterData(1).sigma;

%% Display the interpolated water properties at T_desired (217 °C)
fprintf('Interpolated water properties at T = %d °C:\n', T_steam);
fprintf('rho_l_water = %.6f kg/m^3\n', rho_l_water);
fprintf('rho_v_water = %.6f kg/m^3\n', rho_v_water);
fprintf('h_fg_water  = %.2f J/kg\n', h_fg_water);
fprintf('k_l_water   = %.4f W/(m·K)\n', k_l_water);
fprintf('mu_l_water  = %.8e Pa·s\n', mu_l_water);
fprintf('sigma_water = %.6f N/m\n', sigma_water);
%%

% Mixture property 'pressure' (assumed constant for all mixtures)
pressure = 0.1e6;  % [Pa]

%% --- Prepare to Store Results ---
nMix = height(mixtureData);
results = struct('wtPercent',[], 'finalTw_water',[], 'finalTw_mix',[], 'q',[], 'h_water',[], 'h_mix',[], 'iterCount',[]);

%% --- Loop Over Mixture Conditions ---
for i = 1:nMix
    % Extract mixture-side properties for the i-th row
    wtPct      = mixtureData.wtPercent(i);
    T_mix      = mixtureData.T_mix(i);         % Mixture side wall temperature [°C] (base value)
    rho_l_mix  = mixtureData.rho_l_mix(i);
    rho_v_mix  = mixtureData.rho_v_mix(i);
    h_fg_mix   = mixtureData.h_fg_mix(i);
    k_l_mix    = mixtureData.k_l_mix(i);
    mu_l_mix   = mixtureData.mu_l_mix(i);
    Pr_l       = mixtureData.Pr_l(i);
    sigma_mix  = mixtureData.sigma_mix(i);
    
    % Initial guess for water-side wall temperature (Tw_water) [°C]
    Tw_water_old = 100;
    iter = 0;
    
    % Begin iterative process for current mixture condition
    while iter < maxIter
        iter = iter + 1;
        
        % For water side: define DeltaT_sub = T_steam - Tw_water_old
        DeltaT_sub = T_steam - Tw_water_old;
        
        % Compute water-side Rayleigh number (using your prebuilt function)
        Ra_water = compute_Ra_water(g, theta, rho_l_water, rho_v_water, h_fg_water, L, k_l_water, mu_l_water, DeltaT_sub, sigma_water);
        
        % Compute water-side Nusselt number
        Nu_water = compute_Nu(Ra_water);
        
        % Compute water-side heat transfer coefficient
        h_water = compute_h_water(Nu_water, k_l_water, sigma_water, g, rho_l_water, rho_v_water, theta);
        
        % Compute water-side heat flux (q_dblprime)
        q_dblprime = compute_heat_flux(h_water, DeltaT_sub);
        
        % Mixture-side calculations:
        nu_mix = compute_kinViscosity(mu_l_mix, rho_l_mix);
        la_mix = compute_la_mix(sigma_mix, g, rho_l_mix, rho_v_mix);
        h_mix = compute_h_mix(k_l_mix, la_mix, Pr_l, q_dblprime, rho_v_mix, h_fg_mix, nu_mix, pressure, sigma_mix);
        
        % Compute mixture-side wall temperature (Tw_mix)
        % (Using convective balance: Tw_mix = (q_dblprime / h_mix) + T_mix)
        Tw_mix = (q_dblprime / h_mix) + T_mix;
        
        % Update water-side wall temperature through conduction:
        % (Tw_water_new = (delta * q_dblprime / K_conductivity) + Tw_mix)
        Tw_water_new = (delta * q_dblprime / K_conductivity) + Tw_mix;
        
        % Check convergence
        if abs(Tw_water_new - Tw_water_old) < tol
            break;
        end
        
        % Update for next iteration
        Tw_water_old = Tw_water_new;
    end
    % --- Compute mass flow rates ---
        % 1) Water-side mass flow [kg/(m^2·s)]
        massFlowWater_s = q_dblprime / h_fg_water;
        % 2) Mixture-side mass flow [kg/(m^2·s)]
        massFlowMix_s   = q_dblprime / h_fg_mix;
        
        % Convert to kg/(m^2·hr) if needed
        massFlowWater_hr = massFlowWater_s * 3600;
        massFlowMix_hr   = massFlowMix_s   * 3600;
        
        % Store results for the current mixture condition
        results(i).wtPercent     = wtPct;
        results(i).finalTw_water = Tw_water_new;
        results(i).finalTw_mix   = Tw_mix;
        results(i).q             = q_dblprime;
        results(i).h_water       = h_water;
        results(i).h_mix         = h_mix;
        results(i).iterCount     = iter;
        results(i).mDotWater     = massFlowWater_hr;  % or in s, up to you
        results(i).mDotMix       = massFlowMix_hr;
        
        % Display results for this mixture condition
        fprintf('For %g wt%% mixture:\n', wtPct);
        fprintf('  Converged after %d iterations.\n', iter);
        fprintf('  Final Tw_water = %.4f °C\n', Tw_water_new);
        fprintf('  Final Tw_mix   = %.4f °C\n', Tw_mix);
        fprintf('  Final q (heat flux) = %.4f W/m^2\n', q_dblprime);
        fprintf('  Water-side h   = %.4f W/(m^2·K)\n', h_water);
        fprintf('  Mixture-side h = %.4f W/(m^2·K)\n', h_mix);
        fprintf('  mDotWater      = %.4f kg/(m^2·hr)\n', massFlowWater_hr);
        fprintf('  mDotMix        = %.4f kg/(m^2·hr)\n\n', massFlowMix_hr);
        fprintf('  Ra = %.4f W/(m^2·K)\n', Ra_water);
        fprintf('  Ra = %.4f W/(m^2·K)\n', la_mix);
        fprintf('  Ra = %.4f W/(m^2·K)\n', nu_mix);
        
end

 %% --- Convert Results to a Table and Display ---
resultsTable = struct2table(results);
disp('Final Results Table:');
disp(resultsTable);