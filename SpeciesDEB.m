function [S,pars,chem,mol,yie,etaO,h_O,h_M]=SpeciesDEB (id,Grid)


S=[];

%Basic Parameters (equal for all species)
muX = 525000; %Molar Gibbs energy (chemical potential) of food (J/mol)
muE = 585000; %Molar Gibbs energy (chemical potential) of reserve (J/mol)
muV = 500000; %Molar Gibbs energy (chemical potential) of structure (J/mol)
muP = 480000; %Molar Gibbs energy (chemical potential) of faeces (J/mol)
muN = 244e3/5; %Molar Gibbs energy (chemical potential) of nitrogenous waste (J/mol), synthesis from NH3, Withers page 119
dV = 0.3; %Dry mass fraction of structure
sG = 0.1; %Gompertz stress coefficient
% ha = 2.16e-11; %Weibull ageing acceleration (1/h2)
E_sm = 1116; %Maximum volume-specific energy density of stomach (J/cm3)
% kJ = 0.0009328296; %Maturity maintenance rate coefficient -- rabbit
kapR = 0.95; %Fraction of reproduction energy fixed in eggs
kapX = 0.8; %Digestive efficiency (decimal \%)
n_X = [1,1.8,0.5,.15]; %Chem. indices of C, O, H and N in food
n_E = [1,1.8,0.5,.15]; %Chem. indices of C, O, H and N in reserve
n_V = [1,1.8,0.5,.15]; %Chem. indices of C, O, H and N in structure
n_P = [1,1.8,0.5,.15]; %Chem. indices of C, O, H and N in faeces
E_G = 7843.150; %Cost of structure (J/cm3)
kap_X_P = 0.1; %Faecation efficiency of food to faeces (-)
T_REF = 20+273.15; %Reference temperature for rate correction (deg C)
T_A = 8085; %Arrhenius temperature
T_AL = 18721; %Arrhenius temperature for decrease below lower boundary of tolerance range \code{T_L}
T_AH = 90000; %Arrhenius temperature for decrease above upper boundary of tolerance range \code{T_H}
T_L = 288; %Lower boundary (K) of temperature tolerance range for Arrhenius thermal response
T_H = 315; %Upper boundary (K) of temperature tolerance range for Arrhenius thermal response
T_A2 = 8085; %Arrhenius temperature for maturity maintenance (causes 'Temperature Size Rule' effect)
T_AL2 = 18721; %Arrhenius temperature for decrease below lower boundary of tolerance range \code{T_L}  for maturity maintenance (causes 'Temperature Size Rule' effect)
T_AH2 = 90000; %Arrhenius temperature for decrease above upper boundary of tolerance range \code{T_H} for maturity maintenance (causes 'Temperature Size Rule' effect)
T_L2 = 288; %Lower boundary (K) of temperature tolerance range for Arrhenius thermal response for maturity maintenance (causes 'Temperature Size Rule' effect)
T_H2 = 315; %Upper boundary (K) of temperature tolerance range for Arrhenius thermal response for maturity maintenance (causes 'Temperature Size Rule' effect)
Tb = 37; %Body temperature assumed to be 37 to every individual
d_V = 0.3; %Dry mass fraction of structure
d_E = 0.3; %Dry mass fraction of reserve
d_Egg = 0.3; %Dry mass fraction of egg
fdry = 0.3; %Dry mass fraction of food
s_M = 1;
L_T = 0;

Ehb_ref=10;
Ehx_ref=100;
Ehp_ref=7000;
% ha_ref = 3.907607e-08; %Weibull ageing acceleration (1/d)
ha_ref = 2.03494e-10; %Weibull ageing acceleration (1/d)


%%%%%%%%%%%%%%%%%%%%%%
%%%%% Species ID %%%%%
%%%%%%%%%%%%%%%%%%%%%%

S(1,1) = id;

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DEB PARAMETERS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%

%Species zoom factor
z=lognrnd(0,1);
while z < 0.85 %Minimum species size
    z=lognrnd(0,1);
end
% z=7.7474

%Species kappa (DEB allocation strategy)
S(1,2) = rand;
while S(1,2)<0.4
    S(1,2) = rand;
end
% S(1,2) = 0.43063;

%Species pM (Volume-specific somatic maintenance (J/cm3)) -- Gerado
%com relação a z e kappa com base no padrão de dados empíricos (AmP Tool para as ordens)
%r² = 0.691; a = 8.246040; b(z)=-1.147917; c(kap)=1.986407
pM = 3812.5*(z^-1.149878)*(S(1,2)^1.492113);
% r = (1+1) *rand(1,1) -1;
% pM = 1555.7600 + r*1000;
% pM = 1.41551e+02;

%Species v (Energy conductance (cm))-- Gerado
%com base no padrão de dados empíricos, sem relação com z e kappa(AmP Tool para as ordens)
r = (1+0.6) *rand(1,1) -0.6;
v = 0.0340330 + r*0.02;
% v = 4.44570e-02;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% DEB mass balance-related calculations %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% match H fraction in organics to stated chemical potentials (needed later for heat production)
n_X(2) = ((muX / 10 ^ 5) - 4.3842 * n_X(1) - (-1.8176) * n_X(3) - (0.0593) * n_X(4)) / 0.9823;
n_V(2) = ((muV / 10 ^ 5) - 4.3842 * n_V(1) - (-1.8176) * n_V(3) - (0.0593) * n_V(4)) / 0.9823;
n_E(2) = ((muE / 10 ^ 5) - 4.3842 * n_E(1) - (-1.8176) * n_E(3) - (0.0593) * n_E(4)) / 0.9823;
n_P(2) = ((muP / 10 ^ 5) - 4.3842 * n_P(1) - (-1.8176) * n_P(3) - (0.0593) * n_P(4)) / 0.9823;

% enthalpies (combustion frame)
h_X = 10^5 * (4.3284 * n_X(1) + 1.0994 * n_X(2) + (-2.0915) * n_X(3) + (-0.1510) * n_X(4)); %J mol^(-1)
h_V = 10^5 * (4.3284 * n_V(1) + 1.0994 * n_V(2) + (-2.0915) * n_V(3) + (-0.1510) * n_V(4)); %J mol^(-1)
h_E = 10^5 * (4.3284 * n_E(1) + 1.0994 * n_E(2) + (-2.0915) * n_E(3) + (-0.1510) * n_E(4)); %J mol^(-1)
h_P = 10^5 * (4.3284 * n_P(1) + 1.0994 * n_P(2) + (-2.0915) * n_P(3) + (-0.1510) * n_P(4)); %J mol^(-1)
h_CO2 = 0; %J mol^(-1)
h_O2 = 0; %J mol^(-1)
h_H2O = 0; %J mol^(-1)
n_M_nitro = [1, 2, 1, 2]; % urea
h_N = 631890;
mu_N = 122e3;
h_O = [h_X, h_V, h_E, h_P];
h_M = [h_CO2, h_H2O, h_O2, h_N];
n_O = [n_X; n_V; n_E; n_P]; % matrix of composition of organics, i.e. food, structure, reserve and faeces
CHON = [12, 1, 16, 14];
wO = CHON * n_O;
w_V = wO(2);
M_V = dV / w_V;
y_EX = kapX * muX / muE; % yield of reserve on food
y_XE = 1 / y_EX; % yield of food on reserve
y_VE = muE * M_V / E_G;  % yield of structure on reserve
y_PX = kap_X_P * muX / muP; % yield of faeces on food
y_PE = y_PX / y_EX; % yield of faeces on reserve
nM = [1 0 2 0;0 2 1 0; 0 0 2 0; n_M_nitro];
n_M_nitro_inv = [-1 * n_M_nitro(1) / n_M_nitro(4); (-1 * n_M_nitro(2)) / (2 * n_M_nitro(4)); (4 * n_M_nitro(1) + n_M_nitro(2) - 2 * n_M_nitro(3)) / (4 * n_M_nitro(4)); 1 /n_M_nitro(4)];
n_M_inv = [1 0 -1 0; 0 1/2 -1/4 0; 0 0 1/2 0; n_M_nitro_inv'];
JM_JO = -1 * n_M_inv; %*% n_O
etaO = [y_XE/muE*-1 0 1/muE; y_PE/muE 0 0; -1/muE 0 0; y_VE/muE -1/muE 0];
w_N = CHON * n_M_nitro';
s_M = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Arrhenius temperature correction factor %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tcorr = exp(T_A / T_REF - T_A / (273.15 + Tb)) * (1 + exp(T_AL / T_REF - T_AL / T_L) + exp(T_AH / T_H - T_AH / T_REF)) / (1 + exp(T_AL / (273.15 + Tb) - T_AL / T_L) + exp(T_AH / T_H - T_AH / (273.15 + Tb)));
Tcorr2 = exp(T_A2 / T_REF - T_A2 / (273.15 + Tb)) * (1 + exp(T_AL2 / T_REF - T_AL2 / T_L2) + exp(T_AH2 / T_H2 - T_AH2 / T_REF)) / (1 + exp(T_AL2 / (273.15 + Tb) - T_AL2 / T_L2) + exp(T_AH2 / T_H2 - T_AH2 / (273.15 + Tb)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Compound parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%pM corrigido por temperatura
% S(1,3) = pM * Tcorr;
S(1,3) = pM;

%v corrigido por temperatura
S(1,4) = v;

%pAm
% z=225.7165*(pM^-0.63021)*(S(1,2)^0.78424);
% z=7.7474;
S(1,5) = S(1,3) * z / S(1,2) * 1;

%Em (J/cm³)
S(1,6) = S(1,5) / S(1,4);

%Maximum Structural Volume (V_m, cm³)
S(1,7) = (S(1,2) * S(1,5) / S(1,3)) ^ 3;

%Maximum structural length(L_m, cm)
S(1,8) = S(1,7)^(1/3);

%Maximum length(C_m, cm)
S(1,9) = S(1,8)/0.211; %Tabela 1.1 - Kooijman 2009, DEB 3rd edition - Neotropical Mammals w/tail

%Maximum Reserve (Ef) - Maximum reserve considering the current
%length
S(1,10) = S(1,6)*(S(1,8)^3);

%Neonate Structural Length (L_0 - Derived from R AmP Database relation between L.i &
% %L.b)
%r² = 0.965, a= -0.93619; b=0.92963
S(1,11) = 0.392119*S(1,8)^0.92963;

%Neonate Body Volume (V_0 = L_0^3)
S(1,12) = S(1,11)^3;

%Neonate Initial Reserve (E0) - Maximum reserve considering the current
%length
S(1,13) = S(1,6)*(S(1,11)^3);

%Species EHb (Maturity at birth (J)) -- Ehb = (z/kappa)³Ehb_ref
% S(1,14) = 1206.525*(S(1,11)^2.63992)*(S(1,2)^-4.82910);
% S(1,14) = 71724.14;
S(1,14) = Ehb_ref*((S(1,8)/S(1,2))^3);

%Species Ehx (Maturity at weaning (J)) -- Ehx = (z/kappa)³Ehx_ref
% S(1,15) = 7683.266*(S(1,8)^0.9015)*(S(1,2)^-4.4519)*(S(1,11)^1.3836);
% S(1,15) = 511548.6;
S(1,15) = Ehx_ref*((S(1,8)/S(1,2))^3);

%Species EHp (Maturity at puberty (J)) -- Ehp = (z/kappa)³Ehp_ref
% S(1,16) = 13860.52*(S(1,8)^2.1374)*(S(1,2)^-5.3265);
% S(1,16) = 41124175;
S(1,16) = Ehp_ref*((S(1,8)/S(1,2))^3);

%Species energy invesment ratio (g)
S(1,17) = E_G / (S(1,2) * S(1,6));

%Gestation Time (d) -- ALOMETRICO AmP(Lb) + Kappa + pM (Lika & Kooijman, 2003)
%r² = 0.842; a=2.09751; b(lb)=1.02889;c(kap)=-0.81848;d(pM)=0.21643
S(1,18) = 8.145861*(S(1,11)^1.02889)*(S(1,2)^-0.81848)*(S(1,3)^0.21643);

%Period Between Reproductive Season (Alometrico Panhtheria Interbirth Interval ~ Adult Body Mass)
W_w = S(:,7) * dV + (S(:,6)*S(:,7))*wO(1)/muE;
S(:,19) = round(21.049*(W_w^0.281));
% S(1,19) = 360;
    
%Longevity (d) -- ALOMETRICO AmP! + Kappa + pM (Lika & Kooijman, 2003)
%r²=0.6633; a=6.07502; b(li) = 0.74706; c(kap)= -0.73683; d (pM) = 0.15460
S(1,20) = 434.8582*(S(1,8)^0.74706)*(S(1,2)^-0.73683)*(pM^0.15460);
     
%Species Grid Origin (Grid)
S(1,21)=ceil(rand*Grid);
S(1,22)=ceil(rand*Grid);

%Litter Size -- Gerado
%com relação a tg e kappa com base no padrão de dados empíricos (AmP Tool para as ordens)
%r² = 0.710; a = 5.88702; b(tg)=-1.08402; c(kap)=-0.87643
S(1,23) = 360.33*(S(1,18)^-1.08402)*(S(1,2)^-0.87643);


%kJ(Maturity maintenance rate coefficient) -- Gerado
%com relação a Lm, kappa e pM com base no padrão de dados empíricos (AmP Tool para as ordens)
%r² = 0.784; a = -14.77649; b(kap)=-0.67632; c(pM)=0.98328;d(li)=1.15739
S(1,24) = 3.8252e-07*(S(1,2)^-0.67632)*(S(1,3)^0.98328)*(S(1,8)^1.15739);

%ha (Weibull ageing acceleration (1/h2)) -- ha = z*ha_ref (Kooijman 2009,pag.324)
% S(1,25) = ha_ref*z;
S(1,25) = 1.47864e-05*S(1,8)^-3.5106;


%Parametros calculados diariamente para o Inidividuo
%Species scaled reserve density (e)
% e = E_pres / E_m; 
% 
% %M_V
% M_V = dV / w_V;
% 
% % heating length
% L_T = p_T / p_MT;
% 
% % Relacionado à metabolic heat production
% kap_G = (d_V * mu_V) / (w_V * E_G)

%Parametros imutaveis

%kJ corrigido por temperatura (Maturity maintenance rate coefficient)
% kJ = kJ * Tcorr2;

%ha (Weibull aging acceleration)
% S(1,25) = S(1,25) * Tcorr;

%Relacionados à alimentacao
yEX = kapX * muX / muE;
yXE = 1 / yEX;
yPX = kap_X_P * muX / muP;
mu_AX = muE / yXE;
eta_PA = yPX / mu_AX;

w_X = wO(1);
w_E = wO(3);
w_V = wO(2);
w_P = wO(4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% IMUTABLE PARAMETERS %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Parameters
pars = [muX, muE, muV, muP, muN, dV, sG, E_sm, kapR, kapX, E_G, kap_X_P, T_REF, T_A, T_AL, T_AH, T_L, T_H, T_A2, T_AL2, T_AH2, T_L2, T_H2, Tb, d_V, d_E, d_Egg, fdry, s_M, mu_AX,eta_PA, L_T];
colNames = {'muX', 'muE', 'muV', 'muP', 'muN', 'dV', 'sG', 'E_sm', 'kapR', 'kapX', 'E_G', 'kap_X_P', 'T_REF', 'T_A', 'T_AL', 'T_AH', 'T_L', 'T_H', 'T_A2', 'T_AL2', 'T_AH2', 'T_L2', 'T_H2', 'Tb', 'd_V', 'd_E', 'd_Egg', 'fdry', 's_M', 'mu_AX','eta_PA', 'L_T'};
pars= array2table(pars,'VariableNames',colNames);

%Chemical composition --match H fraction in organics to stated chemical potentials (needed later for heat production)
chem = [n_X; n_E; n_V;  n_P];
rowNames = {'n_X', 'n_E', 'n_V',  'n_P'};
colNames = {'C','O', 'H', 'N'};
chem= array2table(chem,'RowNames',rowNames,'VariableNames',colNames);

%Matrix of composition of organics, i.e. food, structure, reserve and faeces
chem = [n_X; n_E; n_V;  n_P];
rowNames = {'n_X', 'n_E', 'n_V',  'n_P'};
colNames = {'C','O', 'H', 'N'};
chem= array2table(chem,'RowNames',rowNames,'VariableNames',colNames);

%Molecular Weight composition
mol = wO;
colNames = {'C','O', 'H', 'N'};
mol= array2table(mol,'VariableNames',colNames);

%Yield of Food
yie = [y_EX, y_XE, y_VE, y_PX, y_PE];
colNames= {'y_EX','y_XE', 'y_VE', 'y_PX', 'y_PE'};
yie= array2table(yie,'VariableNames',colNames);

%Species metab_mode (Determinate or Indeterminate growth - 0=Indeterminate /
%1 = Determinate)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Life-History Parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%parametros relacionados à reprpdução
% batching = 1 %Modelo de Pecquerie et. al. 2009
% lambda = 1/2 %Relacionado ao batching
% breeding = 0; %Relacionado ao batching
% pXm %Species pXm (Surface area-specific maximum feeding rate (J/cm2)); relacionada ao processo de ingestao
% %kM
% kM = S(1,3) / E_G; %Relacionado ao batching
% 
% %Parametros relacionados à alimentacao
% %pXm
% pXm = pXm * Tcorr * 1;


