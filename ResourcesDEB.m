function [I1,G1] = ResourcesDEB(S,I,G1,prior)
%Funcao que controla a retirada de recurso do ambiente pelos individuos
%Adultos
EH=[];
I1=I;
%Ordem de consumo
if strcmp(prior,'bodysize')
    %Priority order is defined 70% by species body size and 30% by random
    %factors
	ISp=I1(:,[1 2 4]);
    ISp(:,3) = ISp(:,3)/max(ISp(:,3));
	ISp(:,4)=((0.7*rand(size(ISp,1),1))+(0.2*ISp(:,3)))/2;
    [out,idx] = sort(ISp(:,4),'descend');
    I1 = I1(idx,:);
elseif strcmp(prior,'especialists')
    %Specialists have a greater chance to reach food first 
%         [VALEspe,ORDEsp]=sort(ISp(:,5),'descend');
%         ORDEsp=randperm(size(ISp,1));
%         I1=I1(ORDEsp,:);
elseif strcmp(prior,'random')
        I1=I1(randperm(size(ISp,1)),:);
end
           
    %inicia o LOOP de individus para consumo de recursos
    for Ind=1:size(I1,1)
        %Foetus does not forage
        if I1(Ind,17) ~= 0            
        
            %Individual has not yet fed
            if I1(Ind,12)==0

                %Pregnant or lactating females
                if I1(Ind,17) == 4
                    
                    ofs = I(I(:,2) == I1(Ind,2) & I(:,19) == I1(Ind,17),:);
%                     [pup,x] = find(I1(:,2) == ofs(1,2) & I1(:,1) == ofs(1,1));
                    
                    if any(ofs(:,17) == 1) %Females with child
                        
                        %Ingestion Rate
                        IGMax = (0.859*I1(Ind,14)^0.628)*10*1000; %(J/dia)        
%                         IGMax = IGMax*1.5; % Femeas com juvenis se alimentam 150% do normal
                    
                        %DMD(Maximum Daily Movement Distance; m/dia; Garland 1983)
                        DMD = 1.32*((I1(Ind,14)/1000)^0.28)*1000;%Daily Foraging Distance (Garland 1983)
%                         DMD = DMD*1.5;%Femeas com juvenis se alimentam 150% do normals

                        %ICL(Incremental Cost of Locomotion; J/km Garland 1983
                        ICL = (10678*(I1(Ind,14)/1000)^0.7)/1000;
                %         ICL = ((10678*(I1(Ind,15))^0.7)/1000)/1000;

                        %Parametros de Consumo de Recursos
                        Rmax = 0.7483*(I1(Ind,14)/1000)^0.69;%(g/min ; Shipley et al 1994)
                        Ls = 0.54881*(I1(Ind,14)/1000)^0.29;%(Stride length, m; Shipley 1996)
                        Sf = 1.5373*(I1(Ind,14)/1000)^-0.24;%(Stride frequency, steps/s; Shipley 1996)
                        Vmax = (1.03 +0.03*(Ls*Sf));%(m/s; Shipley 1996)
                %         Vmax = 52.16*(I1(Ind,14)/1000)^0.04; %(m/min ; Shipley 1996)
                        hmax = 0.0118*(I1(Ind,14)/1000)^0.03;%(min/bite ; Shipley et al 1994)
                        Smax = 0.0963*(I1(Ind,14)/1000)^0.71;%(g ; Shipley et al 1994)
                        amax = 1.65*(I1(Ind,14)/1000)^-0.17;%(acceleration m/s ; Shipley et al 1996)

                        %Time allocated to foraging (does not scale with Body weigth;
                        %Hudson & white 1985 (Cap 2) & Shipley 2007
                        Tmax = 6; %Maximo de 6h para alimentação.

                        %Consome os recursos
                        [G1,I1(Ind,:),ofs] = moviment_teseDEB(G1,I1(Ind,:),1,DMD,IGMax,ICL,Rmax,Vmax,hmax,Smax,amax,ofs);
                        I(I(:,2) == I1(Ind,2) & I(:,19) == I1(Ind,17),:) = ofs;
                        
                    else %Females with foetus
                        
                        %Ingestion Rate
                        IGMax = (0.859*I1(Ind,14)^0.628)*10*1000; %(J/dia)        
%                         IGMax = IGMax*1.5; % Femeas com juvenis se alimentam 150% do normal
                    
                        %DMD(Maximum Daily Movement Distance; m/dia; Garland 1983)
                        DMD = 1.32*((I1(Ind,14)/1000)^0.28)*1000;%Daily Foraging Distance (Garland 1983)
%                         DMD = DMD*1.5;%Femeas com juvenis se alimentam 150% do normal

                        %ICL(Incremental Cost of Locomotion; J/km Garland 1983
                        ICL = (10678*(I1(Ind,14)/1000)^0.7)/1000;

                        %Parametros de Consumo de Recursos
                        Rmax = 0.7483*(I1(Ind,14)/1000)^0.69;%(g/min ; Shipley et al 1994)
                        Ls = 0.54881*(I1(Ind,14)/1000)^0.29;%(Stride length, m; Shipley 1996)
                        Sf = 1.5373*(I1(Ind,14)/1000)^-0.24;%(Stride frequency, steps/s; Shipley 1996)
                        Vmax = (1.03 +0.03*(Ls*Sf));%(m/s; Shipley 1996)
                %         Vmax = 52.16*(I1(Ind,14)/1000)^0.04; %(m/min ; Shipley 1996)
                        hmax = 0.0118*(I1(Ind,14)/1000)^0.03;%(min/bite ; Shipley et al 1994)
                        Smax = 0.0963*(I1(Ind,14)/1000)^0.71;%(g ; Shipley et al 1994)
                        amax = 1.65*(I1(Ind,14)/1000)^-0.17;%(acceleration m/s ; Shipley et al 1996)

                        %Time allocated to foraging (does not scale with Body weigth;
                        %Hudson & white 1985 (Cap 2) & Shipley 2007
                        Tmax = 6; %Maximo de 6h para alimentação.

                        %Consome os recursos
                        [G1,I1(Ind,:),~] = moviment_teseDEB(G1,I1(Ind,:),1,DMD,IGMax,ICL,Rmax,Vmax,hmax,Smax,amax,0);
                        
                    end
                        
                else
                    %IG_Max (KJ)--> Herbivores 10KJ/g (Nagy 2001; DMI g/day-->Veio do FMR)
                    IGMax = (0.859*I1(Ind,14)^0.628)*10*1000; %(J/dia)        

                    %DMD(Maximum Daily Movement Distance; m/dia; Garland 1983)
                    DMD = 1.32*((I1(Ind,14)/1000)^0.28)*1000;%Daily Foraging Distance (Garland 1983)

                    %ICL(Incremental Cost of Locomotion; J/km Garland 1983
                    ICL = (10678*(I1(Ind,14)/1000)^0.7)/1000;

                    %Parametros de Consumo de Recursos
                    Rmax = 0.7483*(I1(Ind,14)/1000)^0.69;%(g/min ; Shipley et al 1994)
                    Ls = 0.54881*(I1(Ind,14)/1000)^0.29;%(Stride length, m; Shipley 1996)
                    Sf = 1.5373*(I1(Ind,14)/1000)^-0.24;%(Stride frequency, steps/s; Shipley 1996)
                    Vmax = (1.03 +0.03*(Ls*Sf));%(m/s; Shipley 1996)
            %         Vmax = 52.16*(I1(Ind,14)/1000)^0.04; %(m/min ; Shipley 1996)
                    hmax = 0.0118*(I1(Ind,14)/1000)^0.03;%(min/bite ; Shipley et al 1994)
                    Smax = 0.0963*(I1(Ind,14)/1000)^0.71;%(g ; Shipley et al 1994)
                    amax = 1.65*(I1(Ind,14)/1000)^-0.17;%(acceleration m/s ; Shipley et al 1996)

                    %Time allocated to foraging (does not scale with Body weigth;
                    %Hudson & white 1985 (Cap 2) & Shipley 2007
                    Tmax = 6; %Maximo de 6h para alimentação.

                    %Consome os recursos
                    [G1,I1(Ind,:),~] = moviment_teseDEB(G1,I1(Ind,:),1,DMD,IGMax,ICL,Rmax,Vmax,hmax,Smax,amax,0);
%                         if Figures==1
%                             figure(2)
%                             EH(Ind,1) = sum(sum(G1));
%                             plot(EH(:,1))
%                             
%                             figure(3)
%                             imagesc(G1==0)
% 
%                             pause(0.0001)
%                         end
                end
            end
        end
    end   
end%Fecha funcao