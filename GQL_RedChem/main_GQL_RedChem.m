clear all;clc; close all;

reactor_system='isobar';
find_GQL_type = 'local'; % option={'global','local'}
n_GQL_attempt = 1e5; % the max. number of attempt to find out the GQL reduced chemistry
n_GQL_max = 10; % the max. number of GQL candidates for certain dimension
error_GQL = 10; % in percentage
Nf = 2; % dimension of fast invariant subspace

T0=2000; p0=1e5; Phi=1.0;

gas = Solution('./mechanism_H2_Air/Warnatz.cti');
io2 = speciesIndex(gas,'O2');
in2 = speciesIndex(gas,'N2');
ih2 = speciesIndex(gas,'H2');

nsp = nSpecies(gas);

%% calculate the detailed solution for random choice intervall
X0 = zeros(nsp,1);
X0(ih2) = 2 * Phi;
X0(io2) = 1;
X0(in2) = 79/21;

mw = molecularWeights(gas);
set(gas,'Temperature',T0,'Pressure',p0,'MoleFractions',X0);
nsp = nSpecies(gas);
y0 = [temperature(gas)
    massFractions(gas)];
tel = [0 1e+3];

M=eye(nsp+1,nsp+1);

warning('off');

options = odeset('Mass',M,'RelTol',1.e-8,'AbsTol',1.e-10);

out = ode15s(@ode_rhs,tel,y0,options,gas,mw,reactor_system);

%% determine the random choice interval
% the idea is that this interval should be the domain before the beginning
% of the auto-ignition process.
Temp = out.y(1,:);t = out.x;
pos=find(gradient(Temp,t)==max(gradient(Temp,t)));
IDT_detailed = t(pos);
semilogx(out.x,out.y(1,:)); hold on;
w_species_domain = [out.y(2:end,round(pos/pos)) , out.y(2:end,end)];
Temperature_domain = [out.y(1,round(pos/pos)) , out.y(1,end)];
Pressure_domain = [p0 , p0];
%% begin the generation of GQL reduced chemistry
fprintf([num2str(Nf),' dimension(s) will be reduced. A ',...
    num2str(nsp-Nf-3),'-D GQL reduced chemistry will be attempted.\n']);
n_GQL_candidate = 0;
for i = 1 : n_GQL_attempt
    [Ms] =find_GQL_candidate(find_GQL_type,Temperature_domain,Pressure_domain,...
        w_species_domain,gas,nsp,mw,Nf);
    options = odeset('Mass',Ms,'RelTol',1.e-8,'AbsTol',1.e-10);
    out = ode15s(@ode_rhs,tel,y0,options,gas,mw,reactor_system);
    %
    pos=find(gradient(out.y(1,:),out.x)==max(gradient(out.y(1,:),out.x)));
    IDT_GQL_try = out.x(pos);
    if (100*abs(1-IDT_GQL_try/IDT_detailed) < error_GQL) & (out.x(end)==max(tel))
        n_GQL_candidate = n_GQL_candidate + 1;
        semilogx(out.x,out.y(1,:)); hold on; pause(1);
        GQL_candidate(:,:,n_GQL_candidate) = Ms; 
        fprintf(['A possible ',num2str(nsp-Nf-3),'-D GQL reduced chemistry is found. Total: ',...
            num2str(n_GQL_candidate),'\n']);
        if n_GQL_candidate > n_GQL_max
            break;
        end
    end
end

save GQL_Ms.mat GQL_candidate;

