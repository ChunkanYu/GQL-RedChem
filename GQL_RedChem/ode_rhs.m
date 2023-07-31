function dydt = ode_rhs(t, y, gas, mw,reactor_system,PSR_parameter)


% Set the state of the gas, based on the current solution vector.
setMassFractions(gas, y(2:end), 'nonorm');

nsp = nSpecies(gas);

if strcmp(reactor_system,'isochor')
    set(gas, 'T', y(1), 'Rho', density(gas));
    % energy equation
    wdot = netProdRates(gas);
    tdot = - temperature(gas) * gasconstant * (enthalpies_RT(gas) - ones(nsp,1))' ...
        * wdot / (density(gas)*cv_mass(gas));
    
    % species equation
    rrho = 1.0/density(gas);
    for i = 1:nsp
        dwdt(i,1) = rrho*mw(i)*wdot(i);
    end
end

if strcmp(reactor_system,'isobar') | strcmp(reactor_system,'PSR')
    set(gas, 'T', y(1), 'P', pressure(gas));
    % energy equation
    wdot = netProdRates(gas);
    tdot = - temperature(gas) * gasconstant * enthalpies_RT(gas)' ...
        * wdot / (density(gas)*cp_mass(gas));
    
    % species equation
    rrho = 1.0/density(gas);
    for i = 1:nsp
        dwdt(i,1) = rrho*mw(i)*wdot(i);
    end
    
    if strcmp(reactor_system,'PSR')
        hi = (temperature(gas) * gasconstant * enthalpies_RT(gas))./mw;
        tdot = tdot - (1/cp_mass(gas)/PSR_parameter.tau_res) *...
            (PSR_parameter.w_i_unburnt'*hi - PSR_parameter.h_mass_unburnt);
        dwdt = dwdt - (1/PSR_parameter.tau_res) * (y(2:end ) - PSR_parameter.w_i_unburnt) ;
    end
    
end

% set up column vector for dydt
dydt = [ tdot;
    dwdt];

