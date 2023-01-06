function [movx,movy]=geramovrecursoDEB(P,R)
%gera o movimento como uma escolha que maximiza apenas o uso de rescursos
%com uma escolha aleatória entre as quatro opções com maior recurso
   [lin,cols]=size(R);d=zeros(9,1);r=zeros(9,1);
%    imagesc(R==0)
   m=[-1 -1;-1 0; -1 1;1 -1;1 0;1 1;0 -1;0 1;0 0];
   for i=1:9
      O(i,1)=P(1)+m(i,1);
      O(i,2)=P(2)+m(i,2);
      %Corrigir a posicao no grid
      if any(O(i,:)<1) || (O(i,1)>cols) || (O(i,2)>lin)
          if O(i,1)<1 & O(i,2)<1
              O(i,1) = cols;
              O(i,2) = lin;
          elseif O(i,1)>cols & O(i,2)>lin
              O(i,1) = 1;
              O(i,2) = 1;
          elseif O(i,1)<1 & O(i,2)>lin
              O(i,1) = cols;
              O(i,2) = 1;
          elseif O(i,1)>cols & O(i,2)<1
              O(i,1) = 1;
              O(i,2) = lin;
          elseif O(i,1)<1
              O(i,1) = cols;
          elseif O(i,2)<1
              O(i,2) = lin;
          elseif O(i,1)>cols
              O(i,1)=1;
          elseif O(i,2)>lin
              O(i,2)=1;
          end
      end
      r(i)=R(O(i,2),O(i,1));
   end
   r = [O r];
   r1 = sortrows(r,3);
   r1(r1(:,3)<0,3)=0;
   r1(isnan(r1(:,3)),3)=0;
   if sum(r1(:,3))<=0
       targ = randi(size(r1,1),1);
   else
	r1(:,3) = r1(:,3)/sum(r1(:,3));
    r1(:,3) = cumsum(r1(:,3));
    targ = min(find(r1(:,3)>=rand));
   end
   
%    [r,I]=sort(r,'descend');
%    plim=randi(4);
   movx=r1(targ,1);
   movy=r1(targ,2);
   if movx==0 | movy==0
       r1
       targ
       pause
   end