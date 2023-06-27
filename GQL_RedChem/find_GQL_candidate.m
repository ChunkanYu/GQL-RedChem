function [GQL_reduced_chemistry]=find_GQL_candidate(find_GQL_type,...
    Temperature_domain,Pressure_domain,w_species_domain,gas,nsp,mw,Nf)

switch find_GQL_type
    case 'local'
        [T_GQL]=local_Jacobian_algorithm(Temperature_domain,Pressure_domain,w_species_domain,gas,nsp,mw);
    case 'global'
end

[Zs,Zf,Zs_tilde,Zf_tilde] = invariant_subspace(T_GQL,Nf);

GQL_reduced_chemistry = eye(nsp+1,nsp+1);

GQL_reduced_chemistry(2:end,2:end) = [Zs Zf]*[Zs_tilde;0*Zf_tilde];

end

function [T_GQL]=local_Jacobian_algorithm(Temperature_domain,Pressure_domain,...
    w_species_domain,gas,nsp,mw)

% step 1: select randomly ONE set of random species mass fractions
% inside the allowed domain
w_species_attempt = w_species_domain(:,1) + (w_species_domain(:,2) - w_species_domain(:,1)) .* rand(nsp,1);
Temperature_attempt = Temperature_domain(1) + (Temperature_domain(2)-Temperature_domain(1)) * rand(1);
Pressure_attempt = Pressure_domain(1) + (Pressure_domain(2)-Pressure_domain(1)) * rand(1);
% step 2: generate the corresponding Jocobian matrix
for i = 1 : nsp
    dw_species = zeros(nsp,1); dw_species(i,1) = 1e-8 * w_species_attempt(i,1)  + 1e-15;
    %
    setMassFractions(gas, w_species_attempt + dw_species, 'nonorm');
    set(gas, 'T', Temperature_attempt, 'P', Pressure_attempt);
    dwdt_right=chem_source(gas, mw);
    %
    setMassFractions(gas, w_species_attempt - dw_species, 'nonorm');
    set(gas, 'T', Temperature_attempt, 'P', Pressure_attempt,'MassFractions');
    dwdt_left=chem_source(gas, mw);
    % 
    T_GQL(:,i) = ( 0.5/dw_species(i,1) ) * ( dwdt_right - dwdt_left );
end
end

function dwdt=chem_source(gas, mw)

nsp = nSpecies(gas);

wdot = netProdRates(gas);

% species equations
rrho = 1.0/density(gas);
for i = 1:nsp
    dwdt(i,1) = rrho*mw(i)*wdot(i);
end

end

function [Zs,Zf,Zs_tilde,Zf_tilde] = invariant_subspace(T_GQL,Nf)

[E_i,LAMBDA]=eig(T_GQL);

if isreal(LAMBDA)==0
    q = 1;
    for j=1:size(LAMBDA,1)
        if isreal(LAMBDA(j,j))==0
            imag_position(q) = j;
            q = q+1;
        end
    end

    number_complex_pair=size(imag_position,2)/2;

    for j=1:number_complex_pair
        E_i_new1 =  .5*( E_i(:,imag_position(2*j-1)) + E_i(:,imag_position(2*j)) );
        E_i_new2 = -.5*imag( E_i(:,imag_position(2*j-1)) - E_i(:,imag_position(2*j)) );
        E_i(:,imag_position(2*j-1)) = E_i_new1;
        E_i(:,imag_position(2*j)) = E_i_new2;
    end

    LAM = real( LAMBDA );

else
    LAM=LAMBDA;
end

LAM=diag(LAM);E_i_inv=inv(E_i);
[~,idxx]=sort(abs((LAM)));LAM=LAM(idxx,1);

Zs=E_i(:,idxx(1:end-Nf)); Zf=E_i(:,idxx(end-Nf+1:end));
Zs_tilde=E_i_inv(idxx(1:end-Nf),:); Zf_tilde=E_i_inv(idxx(end-Nf+1:end),:);



end