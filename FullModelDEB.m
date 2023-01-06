function FullModelDEB(NSP,N,E1,Grid,tempo,tempo_inv,prior,world,Propagule,NProp,IRep,Figures,DirO,K)

%Parametros iniciais:
%NSP  - numero de especies no inicio do modelo
%N - numero de individuos/especie
%E1 - energia inicial do recurso 1 em todo o sistema
%grid - Numero de celulas na regiao (Grid x Grid)
%tempo - tempo para estabilizar a comunidade nativa(anos)
%tempo_inv - tempo após estabilizacao da comunidade nativa (3650 dias de
%estabilizacao)
%prior - Tipo de especialistas(bodysize/especialists/random)
%world - Neutro / niche
%bsmin - tamamho corporal mínimo na comunidade (kg)
%bsmax - tamamho corporal máximo na comunidade (kg)
%propagule - type of propagule pressure (single/multiple)
%NProp - number of propagules introduced (single=once/multiple=daily for a
%year)
%IRep - number of replicates for different invaders meeting the same native
%community
%Figures - display figures (0/1)
%DirO - Diretorio de outputs
%K - replicate(used for metamodelling)


%CONDICOES INICIAIS
rng(K)
K=num2str(K);
%Number format
format short g
%Species Matrix (S)
S=zeros(NSP,25);

%Conversion to Kelvin Temperature
% T=T+273.15;

%Result matrix
RES=zeros((360*(tempo+tempo_inv)),NSP+1);

%S Matrix Columns
%[Sp_ID  Species Maximum Body Size(g)  FMR(KJ/Day)(Nagy et al 1999)  RaioE1  RaioE2  
% Estrategia(Cresc) Estrategia(Sobrev) POSL POSC POSX POSY NeonatoBodySize TempoGestacao(dias) 
% TempoentrePeriodoReprodutivo Longevidade K(Gompertz) I(Gompertz) FatorCorrecao(Gompertz)]

%Filling Species Matrix(S)
strcat('Creating Species')
for i=1:NSP
    [S(i,:),pars,chem,mol,yie,etaO,h_O,h_M]=SpeciesDEB(i,Grid);
end

%Save table as txt
% T=array2table(S);
% T.Properties.VariableNames = {'ID','Kappa','pM','v','pAM','Em','V_m','L_m','C_m','Ef','L_b','V_b','E_b','Ehb','Ehx','Ehp','g','Tg','r_s','am','posX','posY','Lit_sz','kJ','ha'} 
% writetable(T,'C:\Users\decoa\OneDrive\Doutorado\Bloco2-IBM\Cap3-Competicao\Rerun_DEB\Tabelas\Especies_Simuladas.txt','Delimiter','\t')

%Creating the World
strcat('Creating World')
if strcmp(world,'random')
    G1=rand(Grid,Grid);
    sg1=sum(sum(G1));
    G1=G1/sg1;
%     G1=E1*G1;
    G0=G1;
elseif strcmp(world,'checkerboard')
    
elseif strcmp(world,'landscape')
    G1=zeros(Grid,Grid);
    for ipais=1:round(Grid/3)
        mlp=round(rand*(Grid-1)+1);mcp=round(rand*(Grid-1)+1);
        filhos=Grid*3;
        def=Grid/10;
        for ifilhos=1:filhos
            ml=round(randn*(Grid/def)+mlp);mc=round(randn*(Grid/def)+mcp);

            sigma=rand*0.2*Grid+0.2*Grid;
            for l=1:Grid
                for c=1:Grid
                    G1(l,c)=G1(l,c)+10*exp(-(((l-ml)^2)+(((c-mc)^2)))/sigma);
                end
            end
        end
    end  
    G1=G1/sum(sum(G1));
%     G1=G1*E1;
    G0=G1;
end

% G2=randn(Grid,Grid);
% sg2=sum(sum(G2));
% G2=G2/sg2;
% G2=E2*G2;

% G=G1+G2;
figure(1)
imagesc(G1)
% % figure(2)
% % image(G2)
% % figure(3)
% % image(G)
% pause

%Adult Matrix (A)

%Columns Explanation
%[IND_ID  Sp_ID  Age_IND Body Size_IND IND_FMR  Actual Reserve(KJ/g of Fat)  Maximum Reserve(KJ/g of Fat)
% POSX POSY PosROW PosCOL]
strcat('Creating Individuals')
I=zeros((NSP*N),21);
k=0;
for i=1:NSP
    for j=1:N
        k=k+1;
        I(k,:)=IndividualDEB(S,i,j,pars,mol);
    end
end

%Corrigir Posicao
% A=CheckPOS(A,Grid);

% %Plot do sistema
% j=0;
% for u=1:Grid
%     for p=1:Grid
%         j=j+1;
%         figure(j)
%         temp=A(A(:,10)==p & A(:,11)==u,:);
%         scatter(temp(:,8), temp(:,9),5,temp(:,2))
%     end
% end
% pause

%Matriz SIZE & STRAT
SIZE=[];
STRAT=[];
% EH = [];
CS = [];
AFR=[];

% pp=0;

SAZ = [E1*100 E1*50 E1*20 E1 E1/20 E1/50 E1/100 E1/50 E1/20 E1 E1*20 E1*50];
% SAZ = [E1 E1 E1 E1 E1 E1 E1 E1 E1 E1 E1 E1];

% % DINAMICA % %
%tempo
strcat('Starting Model')
t=0;
for ano = 1:tempo
    for mes = 1:12
        SAZ(mes);
        for dia=1:30
            t=t+1
            STRAT(t,1)=mean(S(I(:,2),2));
            SIZE(t,1)=sum(I(:,17)==0);%foetus
            SIZE(t,2)=sum(I(:,17)==1);%offspring
            SIZE(t,3)=sum(I(:,17)==2);%young
            SIZE(t,4)=sum(I(:,17)==3);%adult
            SIZE(t,5)=sum(I(:,17)==4);%pregnant

            %Replenish Resources
            strcat('Replenishing resources...')
            G1 = G0*SAZ(mes);%Divide resources changing level but mantaining pattern
            
            %Resource Consumption
            strcat('Foraging...')
            [I,G1]=ResourcesDEB(S,I,G1,prior);

            %Resource Allocation
            strcat('Allocating resources...')
        for i = 1:size(I,1)
            [A,I,G1,CS,AFR] = DEB_Mammal(S,I,I(i,:),G1,pars,mol,CS,AFR);
            I(i,:) = A;
        end

            %Atualizar posicoes
            strcat('Updating positions...')
        %     I=CheckPOS(I,Grid);

            %Envelhecimento
            strcat('Individuals Ageing...')
            I(:,11) = I(:,11)+1;

            % Update Individual Stages
            % STD MODEL
            strcat('Updating Stages...')
            for p = 1:size(I,1)
                if I(p,6) < S(I(p,2),14) && I(p,17) ~= 5
                    I(p,17) = 0;
                elseif (I(p,6) > S(I(p,2),14)) && (I(p,6) < S(I(p,2),15)) && I(p,17) ~= 5
                    I(p,17) = 1;
                elseif (I(p,6) > S(I(p,2),15)) && (I(p,6) < S(I(p,2),16)) && I(p,17) ~= 5
                    I(p,17) = 2;
                    I(p,19) = 0;
                elseif I(p,6) >= S(I(p,2),16) && (I(p,17) ~= 4) && (I(p,17) ~= 5)
                    I(p,17) = 3;
                end
            end

            %Zerar a dispersao e a alimentacao
            I(:,12)=0; % Zero p_X
            I(:,13)=0; % Zero f

            %Increase diapause timer
            strcat('Decreasing reproduction timer...')
            if any(I(:,17) == 3 | I(:,17) == 4)
                I(I(:,17) == 3,10)=I(I(:,17) == 3,10)+1;
                I(I(:,17) == 4,10)=I(I(:,17) == 4,10)+1;
                if any(I(:,17)==3 & (I(:,10) > round(S(I(:,2),19)*1.05)))
                    I(I(:,17)==3 & (I(:,10) > round(S(I(:,2),19)*1.05)),10)=0;
                    strcat('Resetou')
                end
            end

            %Increase Gestation Timer
            strcat('Increase Gestation...')
            if any(I(:,17) == 0)
                I(I(:,17) == 0,20) = I(I(:,17) == 0,20)+1;
                %Check individuals gestation period and kills or gives birth  to individuals that
                %trespassed the gestation limit
                [I] = GestationCheckDEB(S,I);
            end

            %Reset females that are no longer with offspring
            strcat('Reset females...')
            if any(I(:,17) == 4)
                I(I(:,17) == 4,10) = 0;
                [rows,~]=find(I(:,17)==4);
                for u=1:size(rows,1)
                    if any(I(:,2)==I(rows(u),2) & I(:,19)==I(rows(u),1)) == 0
                        I(rows(u),17) = 3;
                    end       
                end
        %         [MOM,x] = find(I(I(:,17) == 4,1));
%                 N_MOM = setdiff(I(I(:,17) == 4,1),I(:,19));
%                 I(N_MOM,17) = 3;
        %         I(N_MOM,10) = 0;
            end

            %Remove Dead Individuals (stage 5)
            if any(I(:,17)==5)
                strcat('Remove dead...')
                I=I(I(:,17)~=5,:);
            end

            %Replenish resources
            strcat('Replenishing resources...')
        %     G1=G0;

            %Save species size
            strcat('Save species size...')
            for i=1:NSP+1
                RES(t,i)=size(I(I(:,2)==i,:),1);
            end
            
            %Remove errors
            I=I(I(:,1)~=0,:);

            if Figures==1
        %         strcat('Plot Size...')
                figure(3)
                plot(SIZE(:,1))
                hold on
                plot(SIZE(:,2))
                hold on
                plot(SIZE(:,3))
                hold on
                plot(SIZE(:,4))
                hold on
                plot(SIZE(:,5))
                hold off
                grid
                legend('Foetus', 'Child', 'Young', 'Adults', 'Pregnants')
                pause(0.00001)

                %Plot Strategies
                figure(4)
                plot(STRAT(:,1))
                grid
                legend('Mean k')
        %         strcat('Save EH matrix...')
        %         if sum(I(:,17)==0)>0
        %             EH(t,1) = sum(I(:,17)==0)/sum(I(:,17)==4);
        %         else
        %             EH(t,1) = 0;
        %         end
        %         EH(t,1) = I(I(:,1)==1,9);
        %         EH(t,2) = I(I(:,1)==1,5);
        %         EH(t,3) = I(I(:,1)==1,14);
        %         EH(t,4) = I(I(:,1)==1,18);
        %         EH(t,5) = I(I(:,1)==1,17);
        %         EH(t,6) = I(I(:,1)==1,10);
        %         if any(I(:,1)==21)
        %             pp=pp+1;
        %             EH(pp,1) = I(I(:,1)==21,6);
        %             EH(pp,2) = I(I(:,1)==21,4);
        %             EH(pp,3) = I(I(:,1)==21,11);
        %             EH(pp,4) = I(I(:,1)==21,18);
        % 
        %             figure(5)
        %             plot(EH(:,1))
        %             legend('Maturity')
        %             figure(6)
        %             plot(EH(:,2))
        %             legend('Length')
        %             figure(7)
        %             plot(EH(:,3))
        %             legend('Age')
        %             figure(8)
        %             plot(EH(:,4))
        %             legend('Survival')
        %         end
        %         figure(9)
        %         plot(EH(:,5))
        %         legend('Pregnant')
        %         figure(10)
        %         plot(EH(:,6))
        %         legend('Repr Timer')
        % end
            end
        end %Fecha o ciclo diario pre-invasao
%         figure(2)
%         imagesc(G1)
    end%Fecha o mes
end%Fecha todos os anos

%Write XLS
% xlswrite(strcat(DirO,'NativeCommunity_',K,'.xlsx'),I)
% xlswrite(strcat(DirO,'CommunityStrategy_',K,'.xlsx'),STRAT)
% xlswrite(strcat(DirO,'ClutchSize_',K,'.xlsx'),CS)
% AFR2 = AFR(AFR(:,1)>20,:);
% xlswrite(strcat(DirO,'AFR_',K,'.xlsx'),AFR2)

%Write CSV
csvwrite(strcat(DirO,'NativeCommunity_',K,'.csv'),I)
csvwrite(strcat(DirO,'CommunityStrategy_',K,'.csv'),STRAT)
csvwrite(strcat(DirO,'ClutchSize_',K,'.csv'),CS)
AFR2 = AFR(AFR(:,1)>20,:);
csvwrite(strcat(DirO,'AFR_',K,'.csv'),AFR2)
csvwrite(strcat(DirO,'SpeciesCommunity_',K,'.csv'),S)

% %Plot dda comunidade nativa
% j=0;
% for u=1:Grid
%     for p=1:Grid
%         j=j+1;
%         figure(j)
%         temp=A(A(:,10)==p & A(:,11)==u,:);
%         scatter(temp(:,8), temp(:,9),5,temp(:,2))
%     end
% end
% figure(5)
% plot(RES(:,:))

%Invasao
Result=zeros(IRep,9);
strcat('Invasion Replicates')
for m=1:IRep
    m/IRep
    RES1=RES;
    I2=I;
	S1=S;
    flag_a=0;
    flag_m=0;
    flag_d=0;
    t1=t;
    strcat('Creating Invader Species')
    [S1(NSP+1,:),pars,chem,mol,yie,etaO,h_O,h_M]=SpeciesDEB(NSP+1,Grid);
% 	S1(NSP+1,:)=CreateSpecies(BSMax,BSMin,Mundo,Grid,NSP+1);        
    for ano=(tempo+1):tempo_inv
        for mes=1:12
            for dia=1:30
                t1=t1+1
                %Inserir individuos da especie invasora
                strcat('Creating Invader Individuals')
                if Propagule == 1 %Single invasion event
                    TypeProp=1;
                    if ano==(tempo+1) && mes==1 && dia==1
                        k=size(I2,1);
                        for a=1:NProp
                            k=k+1;
                            I2(k,:)=IndividualDEB(S1,NSP+1,a,pars,mol);
        %                     I2(k,:)=CreateIndividual(a,NSP+1,S1,Class);
                        end
                    end
                end

                if Propagule == 2%Multiple invasion events
                    TypeProp=2;
                    if ano<=(ceil(0.5*tempo_inv)) && mes==1 && dia==1
                        k=size(I2,1);
                        for a=1:NProp
                            k=k+1;
                            I2(k,:)=IndividualDEB(S1,NSP+1,a,pars,mol);
        %                     I2(k,:)=CreateIndividual(a,NSP+1,S1,Class);
                        end 
                    end
                end

                %Replenish Resources
                strcat('Replenishing resources...')
                G1 = G0*SAZ(mes);%Divide resources changing level but mantaining pattern

                %Resource Consumption
                strcat('Foraging...')
                [I2,G1]=ResourcesDEB(S1,I2,G1,prior);

                %Resource Allocation
                strcat('Allocating resources...')
                for i = 1:size(I2,1)
                    [A,I2,G1,CS,AFR] = DEB_Mammal(S1,I2,I2(i,:),G1,pars,mol,CS,AFR);
                    I2(i,:) = A;
                end

                %Atualizar posicoes
                %     strcat('Updating positions...')
                %     I=CheckPOS(I,Grid);

                %Envelhecimento
              strcat('Individuals Ageing...')
                I2(:,11) = I2(:,11)+1;

                % Update Individual Stages
                % STD MODEL
              strcat('Updating Stages...')
                for p = 1:size(I2,1)
                    if I2(p,6) < S1(I2(p,2),14) && I2(p,17) ~= 5
                        I2(p,17) = 0;
                    elseif (I2(p,6) > S1(I2(p,2),14)) && (I2(p,6) < S1(I2(p,2),15)) && I2(p,17) ~= 5
                        I2(p,17) = 1;
                    elseif (I2(p,6) > S1(I2(p,2),15)) && (I2(p,6) < S1(I2(p,2),16)) && I2(p,17) ~= 5
                        I2(p,17) = 2;
                        I2(p,19) = 0;
                    elseif I2(p,6) >= S1(I2(p,2),16) && (I2(p,17) ~= 4) && (I2(p,17) ~= 5)
                        I2(p,17) = 3;
                    end
                end

                %Zerar a dispersao e a alimentacao
                I2(:,12)=0; % Zero p_X
                I2(:,13)=0; % Zero f

                %Increase diapause timer
              strcat('Decreasing reproduction timer...')
                if any(I2(:,17) == 3 | I2(:,17) == 4)
                    I2(I2(:,17) == 3,10)=I2(I2(:,17) == 3,10)+1;
                    I2(I2(:,17) == 4,10)=I2(I2(:,17) == 4,10)+1;
                    if any(I2(:,17)==3 & (I2(:,10) > round(S1(I2(:,2),19)*1.05)))
                        I2(I2(:,17)==3 & (I2(:,10) > round(S1(I2(:,2),19)*1.05)),10)=0;
                        strcat('Resetou')
                    end
                end

                %Increase Gestation Timer
                strcat('Increase Gestation...')
                if any(I2(:,17) == 0)
                    I2(I2(:,17) == 0,20) = I2(I2(:,17) == 0,20)+1;
                    %Check individuals gestation period and kills or gives birth  to individuals that
                    %trespassed the gestation limit
                    [I2] = GestationCheckDEB(S1,I2);
                end
                strcat('Gestation increased')

                %Reset females that are no longer with offspring
                strcat('Reset females...')
                if any(I2(:,17) == 4)
                    I2(I2(:,17) == 4,10) = 0;
                    [rows,~]=find(I2(:,17)==4);
                    for u=1:size(rows,1)
                        if any(I2(:,2)==I2(rows(u),2) & I2(:,19)==I2(rows(u),1)) == 0
                            I2(rows(u),17) = 3;
                        end
                        
                    end
                    
                %     [MOM,x] = find(I2(I2(:,17) == 4,1));
%                     N_MOM = setdiff(I2(I2(:,17) == 4,1),I2(:,19));
%                     I2(N_MOM,17) = 3;
                %     I2(N_MOM,10) = 0;
                end

                 %Remove Dead Individuals (stage 5)
                 if any(I2(:,17)==5)
                     strcat('Remove dead...')
                     I2=I2(I2(:,17)~=5,:);
                 end


                 %Save species size
                 strcat('Save species size...')
                 for i=1:(NSP+1)
                     RES1(t1,i)=size(I2(I2(:,2)==i,:),1);
                 end
            
                 %Check if invader is alive
                 if RES1(t1,NSP+1)==0
                     strcat('Daily Flag Activated')
                     flag_d=1;
                     break
                 end
            end %Fecha o ciclo diario da invasao
            if flag_d==1
                strcat('Monthly Flag Activated')
                flag_m=1;
                break
            end
        end%Fecha o ciclo mensal da invasao
        if flag_m==1
            strcat('Yearly Flag Activated')
            flag_a=1;
            break
        end
    end %Fecha o ciclo da invasao

    strcat('NCom')
    NCom=I2(I2(:,2)==(NSP+1),:);
    if NCom==0
        NCom=0;
    else
        NCom=sum(sum(crosstab(NCom(:,15),NCom(:,16))>0));
    end
    %Matriz de Resultados
    %Abundancia , Numero Comunidades Ocupadas ,Volumetric Length , Especialization(k) , OriginX,
    %OriginY, Tipo de Propagulo (Single/Multiple), Numero de Propagulos,
    %Tempo que  especie invasora ficou no sistema
    strcat('Result Matrix')
    Result(m,:)=[RES1(t1,NSP+1), NCom ,S1(NSP+1,8), S1(NSP+1,2), S1(NSP+1,21), S1(NSP+1,22), Propagule, NProp, t1];
    csvwrite(strcat(DirO,'AbundanciaComunidade_',K,'_Invasor_',num2str(m),'.csv'),RES1)
    disp(Result(1:m,:));
end %Fecha as replicas da invasao

%Write XLS
% xlswrite(strcat(DirO,'Result_',K,'.xlsx'),Result)
% % xlswrite(strcat(DirO,'Abundancia_',K,'.xlsx'),RES1)
% xlswrite(strcat(DirO,'SpeciesCommunity_',K,'.xlsx'),S)

%Write CSV
csvwrite(strcat(DirO,'Result_',K,'.csv'),Result)
% xlswrite(strcat(DirO,'Abundancia_',K,'.xlsx'),RES1)
csvwrite(strcat(DirO,'SpeciesCommunity_',K,'.csv'),S)