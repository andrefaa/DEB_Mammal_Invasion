function I = GestationCheckDEB (S,I)

%Check if individual surpassed the gestation period and is not
%ready to be born --> Dead
if sum(I(:,17) == 0)>0
	idx = find(I(:,17) == 0);
	for i = 1:size(idx,1)
        if I(idx(i),20) >= S(I(idx(i),2),18)
            if I(idx(i),6) >= 0.8*S(I(idx(i),2),14)
                I(idx(i),17) = 1;
                I(idx(i),20) = 0;
%                 I(idx(i),19) = 0;
                I(idx(i),18) = I(idx(i),6)/S(I(idx(i),2),14); %Decrease survival probability
                
            else
                I(idx(i),17) = 5; %Offspring is not born
%                 I(idx(i),19) = 0;
            end
        end
	end
end