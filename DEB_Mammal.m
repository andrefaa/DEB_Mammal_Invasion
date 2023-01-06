function [A,I,G1,CS,AFR] = DEB_Mammal(S,I,A,G1,pars,mol,CS,AFR)

% Unpack variables
V = A(:,3);% cm^3, structural volume
L = A(:,4); % cm, structural length
E = A(:,5);% J/, reserve density
H = A(:,6);% J, maturity
q = A(:,7);% -, aging acceleration
hs = A(:,8);% -, hazard rate
R = A(:,9);% J, reproduction buffer energy

v = S(A(:,2),4); %cm, energy conductance
kap = S(A(:,2),2); %-,allocation fraction to soma
V_m = S(A(:,2),7); % cm ^ 3, maximum structural volume
g = S(A(:,2),17); %energy investment ratio

e = (A(:,5)/A(:,3)) / S(A(:,2),6);  % -, scaled reserve density

if(A(:,6) < S(A(:,2),14)) % embryo (Equations from Desforges et. al. 2019)
    
    %MOM Characteristics
    MOM = I(I(:,2)==A(:,2) & I(:,1)==A(:,19),:);
    if size (MOM,1) == 0
        A(:,17)= 5;
        dL = 0;
        dV = 0;
        dE = 0;
        dH = 0;
        dEs = 0;
        dq = 0;
        dhs = 0;
        dR = 0;
    else
        %Somatic Maintenance
        dL = S(A(:,2),4)/3;
    %     LN = L+dL;
        p_S = S(A(:,2),3)*L^3;

        %Reserve mobilization
        p_C = E * ((v*pars.('E_G')*(L^2)+p_S)/(kap*E+pars.('E_G')*(L^3)));

        %Energy entering the foetus (demand system and related to mom energy
        %density)
        Ed = MOM(1,5)/MOM(1,3);
        p_E = S(A(:,2),4)*Ed*A(:,4)+p_C;

        % structure
        dL = S(A(:,2),4)/3;
        dV = (A(:,4)+dL)^3- A(:,3);
    %     LN = L+dL;
        p_G = pars.('E_G')*v*L^2;

        %Energy Flux Allocated to Soma
        p_K = (p_S + p_G)/kap;

        %Maturity maintenance & change in maturity
        p_J = S(A(:,2),24) * H;

        %Maturance
        p_R = (1 - kap) * p_K - p_J;

        %Reserve Dynamics
        dE = p_E - p_C;

        %Remaining Metabolic Needs
        p_M = p_K - p_C;

        %Total Energy Flux to Foetus
        p_F = p_M + p_E;

        %Cost for Mother (change in reproduction buffer)    
        Gcost = p_F/pars.('kapR');

        if MOM(1,9) > Gcost %Costs extracted from reproduction buffer
            MOM(1,9) = MOM(1,9) - Gcost;
            I(I(:,2)==A(:,2) & I(:,1)==A(:,19),:) = MOM;
        else % If costs are not supplied by mother, the foetus dies
            A(:,17) = 5;
        end

        % no aging or stomach in embryo
        dH = p_R;
        dEs = 0;
        dq = 0;
        dhs = 0;
        dR = 0;
    end

else
        
    %Feeding, Assimilation & Gut Content
    if (A(:,6) >= S(A(:,2),14) & A(:,6) < S(A(:,2),15)) %Lactating
        
        %MOM Characteristics
        MOM = I(I(:,2)==A(:,2) & I(:,1)==A(:,19),:);
        
        %Milk percentage
        if A(:,6)<(0.25*S(A(:,2),15))
%         if A(:,11) < 30 %(PENSAR EM ESTRATEGIA PARA O LIMIAR DE 30!)
            kmilk = 1;
        else
            A(:,21) = min(A(:,21),A(:,11));
            kmilk = exp(-0.0124*(A(:,11)-A(:,21)));
        end
        
        if size(MOM,1) ~= 0
        
            %Mom Energy Density
            Ed = MOM(1,5)/MOM(1,3);              

            %Assimilation
            fmilk = 3;
            p_Amilk = kmilk*fmilk*S(A(:,2),5)*A(:,4)^2;

            %Lactation costs for mother
            krl=0.95;%Lactation eficiency
            kL = 0.8; %milk digestion efficiency
            Lcost = p_Amilk/krl*kL;

            if MOM(1,9) > Lcost %Costs extracted from reproduction buffer
                MOM(1,9) = MOM(1,9) - Lcost;
                I(I(:,2)==A(:,2) & I(:,1)==A(:,19),:) = MOM;
            elseif (MOM(1,5) > Lcost) & ((Ed / S(MOM(1,2),6)) > 0.8) %Costs extracted from reserves
                MOM(1,5) = MOM(1,5) - Lcost;
                I(I(:,2)==A(:,2) & I(:,1)==A(:,19),:) = MOM;
            else %Mother stops milk production
                p_Amilk = 0;        
            end
        else %Mother is dead, no energy from milk
%             strcat('Mom is dead....=(')
            p_Amilk = 0;
        end
        
        %Foraging Supplement
        p_Aforage = (1-kmilk)*A(:,12)*pars.('kapX');
        p_A = p_Amilk + p_Aforage;
        dEs = A(:,12)-p_Aforage;
        
    else %Adults
        p_X = A(:,12);%pXm = pAm/kX (Kooijman 2009, Cap8, p309)
        p_A = p_X*pars.('kapX'); %pA = pX*kX (Kooijman 2009, Eq. 2.3, p.36)
    %     p_A = S(A(:,2),5) * A(:,14) * L ^ 2;
        dEs = p_X - p_A;
    end
    
    %Somatic Maintenance
    p_S = S(A(:,2),3)*(L^3);
    
    %Pregnant females reduce 5-10% of somatic maintenance (reduce activity)
    %Kooijman 2009, p.283
%     if A(:,17) == 4
%         p_S = 0.9*p_S;
%     end
    
    %Reserve Mobilization
    p_C = E * ((v*pars.('E_G')*(L^2)+p_S)/(kap*E+pars.('E_G')*(L^3)));
    
    %Energy Reserves
    dE = p_A - p_C;
        
    %Growth
    p_G = kap*p_C-p_S;
    if p_G < 0 %No shrinking
        p_G = 0;
    end
    dV = p_G/pars.('E_G');
    
    %Maturity & Reproduction Buffer
    p_J = S(A(:,2),24) * H;
    p_R = (1-kap)*p_C - p_J;
       
	% structure and starvation (assimilated energy < energy required for
	% somatic maintenance)
    if (p_A < p_S)
    	dS = p_S - p_A; % J / t, starvation energy to be subtracted from reproduction buffer if necessary
%         dS = V * r * -1 * pars.('muV') * pars.('d_V') / mol.(3); % J / t, starvation energy to be subtracted from reproduction buffer if necessary
        dV = 0;
        % draw from reproduction or reserves under starvation
        if(R >= dS)
            A(:,9) = R - dS;
        elseif (R < dS && p_R >= dS)
            p_R = p_R - dS;
        elseif (R < dS && p_R < dS)
            E = E - dS;
            A(:,5) = E;
        end
        
        %Individuals have a chance of moving to a better location to search
        %for food
        Porigin = A([15,16]);
        P1=Porigin;
        MaxDist=ceil((1.07*(A(:,14)/1000)^0.68));
        z=1;
        [lin,cols] = size(G1);
        if MaxDist == 1
            P1 = Porigin;%Individual does not leave cell
        else
            for i = 1:MaxDist
                z=z+1;
                [movx,movy]=geramovrecursoDEB(P1(i,:),G1);
                P1(z,1)=movx;
                P1(z,2)=movy;

                %Round-earth correction
                if any(P1(z,:)<1) | (P1(z,1)>cols) | (P1(z,2)>lin)
                    if P1(z,1)<1 & P1(z,2)<1
                        P1(z,1) = cols;
                        P1(z,2) = lin;
                    elseif P1(z,1)>cols & P1(z,2)>lin
                        P1(z,1) = 1;
                        P1(z,2) = 1;
                    elseif P1(z,1)<1 & P1(z,2)>lin
                        P1(z,1) = cols;
                        P1(z,2) = 1;
                    elseif P1(z,1)>cols & P1(z,2)<1
                        P1(z,1) = 1;
                        P1(z,2) = lin;
                    elseif P1(z,1)<1
                        P1(z,1) = cols;
                    elseif P1(z,2)<1
                        P1(z,2) = lin;
                    elseif P1(z,1)>cols
                        P1(z,1)=1;
                    elseif P1(z,2)>lin
                        P1(z,2)=1;
                    end
                end
            end
        end
        A(:,15) = P1(1);
        A(:,16) = P1(2);
    end
    
    %Maturity/Reproduction allocation
    if(H < S(A(:,2),16))
        dH = p_R;
        dR = 0;
    else
        dH = 0;
        dR = p_R * pars.('kapR');
    end

	% Ageing (equation 6.2 in Kooijman 2010 (DEB3)
    dL =  ((A(:,3) + dV)^(1/3)) - A(:,4);
    r = 3/L * dL;
	dq = (q * (V / V_m) * pars.('sG') + S(A(:,2),25)) * e * ((v / L) - r) - r * q; % ageing acceleration
	dhs = q - r * hs; % hazard
end


%Update Individual Parameters
A(:,3) = A(:,3) + dV; %change in Volume (cm³)
A(:,4) = A(:,3)^(1/3); %change in length (cm)
A(:,5) = A(:,5) + dE; %Change in energy density (J/cm³)
if A(:,5) > S(A(:,2),6)*A(:,4)^3
    resid = A(:,5) - (S(A(:,2),6)*A(:,4)^3);
    A(:,5) = S(A(:,2),6)*A(:,4)^3;
else
    resid = 0;
end

A(:,6) = A(:,6) + dH; %Change in maturity (J)
A(:,7) = A(:,7) + dq; %Change in aging acceleration
A(:,8) = A(:,8) + dhs; %Change in hazard rate
A(:,9) = A(:,9) + dR; %Change in reproduction energy buffer(J)
if A(:,9) < 0
    A(:,9) = 0;
end

%Residuals (faeces) back to environment (10% of the energy)
G1(A(:,16),A(:,15)) = G1(A(:,16),A(:,15)) + (dEs+resid)*0.1;

% New value for scaled reserve density
e = (A(:,5)/A(:,3)) / S(A(:,2),6);  

%Reproduction
if ((A(:,10) > round((S(A(:,2),19)*0.95))) & (A(:,10) < round((S(A(:,2),19)*1.05))) & A(:,17) == 3) %Check if invidivual is in the breeding season
%     if(A(:,9) > S(A(:,2),13)) %Check if reproduction buffer has enough energy for a single offspring
%     Ucum = ceil(S(A(:,2),12)*((A(:,13)+S(A(:,2),17))/S(A(:,2),4))*(1+3/4*((S(A(:,2),11)/S(A(:,2),8))/A(:,13)))*(0.2*(A(:,4)/S(A(:,2),8))));
    E_ofs = (S(A(:,2),11)^3)*(A(:,5)/A(:,3));
    if A(:,9) >= 0.8*E_ofs
        O = [];
%       clutchsize = ceil(A(:,9)/Ucum);
%         E_ofs = (A(:,5)/A(:,3))*(S(A(:,2),11)^3);
        clutchsize = ceil(A(:,9)/E_ofs);
        
        %Get Species Clutchsize
        CS = [CS ;A(:,2) clutchsize];
        
        %Get species age of first reproduction (AFR)
        if size (AFR,1) > 0 
            if (AFR(:,1) == A(:,1) & AFR(:,2) == A(:,2)) == 0
                AFR = [AFR; A(:,1) A(:,2) A(:,11)];
            end
        else
            AFR = [AFR; A(:,1) A(:,2) A(:,11)];
        end
        
        
%         pause(0.1)
%         clutchsize = ceil(S(A(:,2),23));
%         E_ofs = (0.0001^3)*(A(:,5)/A(:,4));
        
        clutchenergy = clutchsize * ((0.0001^3)*(A(:,5)/A(:,3)));

        A(:,9) = A(:,9) - clutchenergy;
        A(:,17) = 4; %Individual is pregnant (Stage 4)
        A(:,10) = 0; %Reset diapause timer
        
        %Add foetus to individuals matrix
        for i = 1:clutchsize
            ofs_id = max(I(I(:,2)==A(:,2),1))+1;
            mom_id = A(:,1);
            I = [I; ofs_id A(:,2) (0.0001^3) 0.0001 E_ofs 0 0 0 0 0 0 0 1 0 A(:,15) A(:,16) 0 1 mom_id 0 100000];
        end
    end
elseif (A(:,10) > (S(A(:,2),19)*1.05) & A(:,17) == 3)
        A(:,9) = A(:,9)*0.8; %Lose 20% of the energy invested in that gestation period (arbitraty!)
        A(:,10) = 0; %Reset diapause timer
end
   
%WetMass 
%  wetgonad = ((A(:,6)/ pars.('muE'))* mol.(3))/pars.('d_Egg');
%  wetstorage = A(:,5)* mol.(3)/ pars.('muE');
%  A(:,14) = A(:,3) * pars.('dV') + wetgonad + wetstorage;
%  
 A(:,14) = A(:,3) * pars.('dV') + (A(:,5)+A(:,9))*mol.(1)/pars.('muE');

%Survival
dsurvdt = -1 * A(:,18) * A(:,8);
A(:,18) = A(:,18) + dsurvdt;
if (rand-0.4)>A(:,18) || e<0.4 %Death due to ageing process or low reserve density (0.4 threshold from Desforges et. al. 2019
    A(:,17) = 5;
    if any (I(:,2)==A(:,2) & I(:,19)== A(:,1) & I(:,17) == 0)
        I(I(:,2)==A(:,2) & I(:,19) == A(:,1) & I(:,17) == 0,17) = 5;
    end
end

end