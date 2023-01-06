function [R,I,ofs] = moviment_teseDEB(R,I,nograph,DMD,IGMax,ICL,Rmax,Vmax,hmax,Smax,amax,ofs)
% R - é um mapa da distribuição da quantidade de recursos
% I - é o indivíduo que ira forragear
% nograph - plota(0) ou não os scatterplots
% DMD - Total Daily Movement Distance of an individual
% IGMax - Maximum energy that can be ingested in a day
% ICL - Incremental cost of locomotion (energy lost in the process of
% foraging)
% Rmax - Food intake(g) per time(min)
% Vmax - Individual foraging velocity (m/min)
% hmax - Time taken in each bite (min/bite)
% Smax - Maximum bite size (g)

if ofs ~= 0
    
    %Ingestion Rate
    IGOfs = (0.859*ofs(1,14)^0.628)*10*1000; %(J/dia)
    %Proportion of food coming from mother's milk
    if ofs(:,11) < 30 %(PENSAR EM ESTRATEGIA PARA O LIMIAR DE 30!)
            kmilk = 1;
        else
            kmilk = exp(-0.0124*(ofs(:,13)-10));
    end
    IGOfs = IGOfs * (1-kmilk);
    
    %Movement and cost of movement
    DMDOfs = DMD;
    ICLOfs = (10678*(ofs(1,14)/1000)^0.7)/1000;

	%Parametros de Consumo de Recursos
	RmaxOfs = 0.7483*(ofs(1,14)/1000)^0.69;%(g/min ; Shipley et al 1994)
	LsOfs = 0.54881*(ofs(1,14)/1000)^0.29;%(Stride length, m; Shipley 1996)
	SfOfs = 1.5373*(ofs(1,14)/1000)^-0.24;%(Stride frequency, steps/s; Shipley 1996)
	VmaxOfs = (1.03 +0.03*(Ls*Sf));%(m/s; Shipley 1996)
	hmaxOfs = 0.0118*(ofs(1,14)/1000)^0.03;%(min/bite ; Shipley et al 1994)
	SmaxOfs = 0.0963*(ofs(1,14)/1000)^0.71;%(g ; Shipley et al 1994)
	amaxOfs = 1.65*(ofs(1,14)/1000)^-0.17;%(acceleration m/s ; Shipley et al 1996)
end

if nograph==0
    figure(1)
    imagesc(R);
end
[lin,cols]=size(R);

%Definir pontos de origem e destino
Porigin = I([15,16]);
Ptarget=Porigin;

%o ponto de saida é P origin
mov=1;
% vid=VideoWriter('C:\Users\decoa\OneDrive\Doutorado\Apresentacao\Teste_Video\Caminhos12.avi');
% vid.FrameRate=10;vid.Quality=100;
% open(vid)
if nograph==0
    imagesc(R);
    hold on
    scatter(Porigin(1),Porigin(2),90,[1 0 0],'filled')
    scatter(Ptarget(1),Ptarget(2),90,[1 0 0],'s','filled')
%     IM=getframe;
% 	writeVideo(vid,IM);
end

%Primeiro movimento
P1=Porigin;
DMDt=DMD/6;%Distancia percorrida em uma hora
%Pregnant Females eat more
if I(1,17) == 4
    IGMaxt=IGMax*1.25;
else
    IGMaxt=IGMax;
end

if ofs ~= 0
    IGOfst=IGOfs;
end
z=1;
R1=R;%Matriz de recursos espelhada
F=[];%Matrix for resource consumption
if ofs ~= 0
    OF=[];
end

for k = 1:6 %6h foraging activity
    %Temporary R
    R1=R;    
    %Movement in 1h window
    steps = ceil(DMDt/1000);%Number of cells visited in 1 hour
    if steps == 1
        P1 = Porigin;%Individual does not leave cell
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
    else
        for i = 1:steps
            z=z+1;
            [movx,movy]=geramovrecursoDEB(P1(i,:),R1)
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

    %Processo de forrageio entre as duas celuals (distancia percorrida
    %pelo individuo(m) e energia consumida
    [R1,IG_Step] = forage_cellsDEB(R1,P1,I,DMDt,IGMaxt,ICL,Rmax,Vmax,hmax,Smax,amax);
    
    %Offspring feeding
    if ofs ~= 0
        [R1,IG_Step_Ofs] = forage_cellsDEB(R1,P1,ofs,DMDt,IGOfst,ICLOfs,RmaxOfs,VmaxOfs,hmaxOfs,SmaxOfs,amaxOfs);
    end

    if nograph==0 
        scatter(P1(:,1),P1(:,2),[],[1 0 0],'filled')
%         scatter([Porigin(1),P1(z,1)],[Porigin(2),P1(z,2)],[],[1 0 0],'filled')
    %     IM=getframe;
    % 	writeVideo(vid,IM);
    % 	pause(1)
    end
    
    %Accumulated Resources matrix
    F(k,1)=IG_Step;%quantidade de recurso ingerida (J)
    F(k,2)=DMDt;%distancia percorrida em cada tempo (m)
    F(k,3)=sum(F(:,1));%Recursos acumulados (J)
    F(k,4)=sum(F(:,2));%Movimento acumulado (m)
    
    if ofs ~= 0
        OF(k,1)=IG_Step_Ofs;%quantidade de recurso ingerida (J)
        OF(k,2)=DMDt;%distancia percorrida em cada tempo (m)
        OF(k,3)=sum(F(:,1));%Recursos acumulados (J)
        OF(k,4)=sum(F(:,2));%Movimento acumulado (m)
    end
    
    %Reset P1, z and R
    R=R1;
    P1 = P1(z,:);
    z=1;
    
    %Break if Individual ingested its maximum amount of food (IGMax)
    if F(k,3) >= IGMaxt
        E_ind = IGMaxt;
        if ofs ~= 0 
            E_ofs = IGOfs;
        end
        break
    else
        E_ind = F(k,3);
        if ofs ~= 0 
            E_ofs = OF(k,3);
        end
    end
end

%Ingested food & Final Position
Ind_X = P1(z,1);
Ind_Y = P1(z,2);

%Atualiza energia consumida, posicao, budget de movimento e de ingestão
%Energia extra das fêmeas vai direto para o orçamento reprodutivo
if E_ind > IGMax
    I(12) = I(12)+(E_ind-IGMax);
    E_ind = E_ind-IGMax;
end
I(15) = Ind_X;
I(16) = Ind_Y;
I(12) = E_ind;
I(13) = E_ind/IGMax; 
if I(13)<0
    I(13)=0;
end

if ofs ~= 0
    ofs(15) = Ind_X;
    ofs(16) = Ind_Y;
    ofs(12) = E_ofs;
    ofs(13) = E_ofs/IGOfs;
    if ofs(13)<0
        ofs(13)=0;
    end
end

if nograph==0
    plot(P1(:,1),P1(:,2),'r-')
    hold off
    figure(2)
    plot(F(:,4),F(:,3),'-ro')
    xlabel('Distância percorrida');ylabel('Recursos acumulados')
end