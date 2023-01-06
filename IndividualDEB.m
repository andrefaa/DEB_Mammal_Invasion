function [I]=IndividualDEB (S,sp,ind,pars,mol)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Individual parameters %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I=[];

%Individual ID
I(1,1) = ind;

%Species ID
I(1,2) = S(sp,1);

%Individual volume (cm³)
I(1,3) = S(sp,7);
% I(1,3) = S(sp,12);

%Individual structural length (cm)
I(1,4) = S(sp,8);
% I(1,4) = S(sp,11);

%Individual reserve (J)
% I(1,5) = S(sp,10);
% I(1,5) = S(sp,13);
% I(1,5) = I(1,5)/I(1,3);
I(1,5) = S(sp,6)*I(1,4)^3;

%Species scaled reserve density (e)
e = I(1,5) / S(sp,6); 

%Individual maturity (J)
I(1,6) = S(sp,16);%Adult
% I(1,6) = S(sp,15);%Young
% I(1,6) = S(sp,14);%Offspring

%Individual aging acceleration(-)
I(1,7) = 0;

%Individual hazard rate(-)
I(1,8) = 0;

%Individual Reproduction Buffer Energy (J)
I(1,9) = 0;

%Individual Reproduction Timer (d)
I(1,10) = 0;

%Individual Age (d)
I(1,11) = round(0.5*S(sp,20));

%Individual Ingested Food (pX; J)
I(1,12) = 0;

%Individual Scaled Food Response (f)
I(1,13) = 0;

%Individual wet weigth (g)
%  wetgonad = ((I(1,9) / pars.('muE')) * mol.(3)) / pars.('d_Egg');
%  wetstorage = ((I(1,3) * I(1,5) / pars.('muE')) * mol.(3)) / pars.('d_E');
%  I(1,14) = I(1,3) * 1 + wetgonad + wetstorage;
%  wetgonad = ((I(1,6)/ pars.('muE'))* mol.(3))/pars.('d_Egg');
%  wetstorage = I(1,5)* mol.(3)/ pars.('muE');
%  I(1,14) = I(1,3) * pars.('dV') + wetgonad + wetstorage;
%  
 I(1,14) = I(1,3) * pars.('dV') + (I(1,5)+I(1,9))*mol.(1)/pars.('muE');

%Individual Grid Position (Grid)
I(1,15)=S(sp,21);
I(1,16)=S(sp,22);

%Individual Stage (0=faetus/ 1=child/ 2=young/ 3=adult/ 4=pregnant)
I(1,17)=3;

%Individual Survival Probability (%)
I(1,18) = 1;

%Individual MOM ID
I(1,19) = 0;

%Gestation Timer
I(1,20) = 0;

%Lactation Exclusive Limit
I(1,21) = 100000;