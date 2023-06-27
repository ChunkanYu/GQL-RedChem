function dydt = ode_rhs(t, y, gas, mw,reactor_system)


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

if strcmp(reactor_system,'isobar')
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


end

% set up column vector for dydt
dydt = [ tdot;
    dwdt];
