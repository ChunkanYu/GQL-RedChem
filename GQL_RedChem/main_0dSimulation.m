clear all;clc;

reactor_system='isobar';
% chemistry='detailed_chemistry';
chemistry = 'GQL_chemistry';
% chemistry = 'QSSA_chemistry'; 
QSS_species = {'OH','O'};

T0=2000; p0=1e5; Phi=1.0;

gas = Solution('./mechanism_H2_Air/Warnatz.cti');
% gas = Solution('./mechanism_H2_Air/ELTE2014.cti');
% gas = Solution('./mechanism_H2_Air/grimech30.cti');
% gas = Solution('./mechanism_H2_Air/Keromnes.cti');
% gas = Solution('./mechanism_H2_Air/OConaire.cti');
io2 = speciesIndex(gas,'O2');
in2 = speciesIndex(gas,'N2');
ih2 = speciesIndex(gas,'H2');


nsp = nSpecies(gas);
X = zeros(nsp,1);
X(ih2) = 2 * Phi;
X(io2) = 1;
X(in2) = 79/21;



mw = molecularWeights(gas);
set(gas,'Temperature',T0,'Pressure',p0,'MoleFractions',X);
nsp = nSpecies(gas);
y0 = [temperature(gas)
    massFractions(gas)];
tel = [0 1e+3];

Ms=eye(nsp+1,nsp+1);

switch chemistry
    case 'detailed chemistry'
    case 'GQL_chemistry'
        Ms = importdata('GQL_Ms.mat');
    case 'QSSA_chemistry'
        for i = 1 : size(QSS_species,2)
            iQSS_species = speciesIndex(gas,QSS_species{i});
            Ms(iQSS_species+1,iQSS_species+1) = 0;
        end
end



options = odeset('Mass',Ms(:,:,1),'RelTol',1.e-8,'AbsTol',1.e-10);

out = ode15s(@ode_rhs,tel,y0,options,gas,mw,reactor_system);

Temperature_results = out.y(1,:);
time = out.x;
pos=find(gradient(Temperature_results,time)==max(gradient(Temperature_results,time)));
IDT=time(pos(1,1));

% plot the temperature and OH mole fractions.
% figure(1);
% plot(out.x,out.y(1,:));
% xlabel('time');
% ylabel('Temperature');
% title(['Final T = ' num2str(out.y(1,end)) ' K']);
%
% figure(2);
% ioh = speciesIndex(gas,'OH');
% plot(out.x,out.y(1+ioh,:));
% xlabel('time');
% ylabel('Mass Fraction');
% title('OH Mass Fraction');
