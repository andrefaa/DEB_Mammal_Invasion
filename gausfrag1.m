function gausfrag(aresta,p,pais,filhos,def,filename)

D=zeros(aresta,aresta);
for ipais=1:pais
    mlp=round(rand*(aresta-1)+1);mcp=round(rand*(aresta-1)+1);
    for ifilhos=1:filhos
        ml=round(randn*(aresta/def)+mlp);mc=round(randn*(aresta/def)+mcp);
        
        sigma=rand*0.2*aresta+0.2*aresta;
        for l=1:aresta
            for c=1:aresta
                D(l,c)=D(l,c)+10*exp(-(((l-ml)^2)+(((c-mc)^2)))/sigma);
            end
        end
    end
end
D(D<30)=0;
image(D)
Dor=D(D>=0);
[min(Dor) max(Dor)]
%figure(1)
%image(D)
%figure(2)
D1=D;
for i=1:length(p)
    qp=prctile(Dor,100*p(i))
    D2=(D1>=qp)*1;
    [D,f]=frag_stats(D2);
    f
    %file=[filename num2str(round(100*p(i)))];
    [lin,cols]=size(D);minx=1;miny=1;cellsize=1;empvalue=-9999;
    %writeasc(D,1,file,'.asc',minx,miny,cellsize) ;
    %save(file,'D','cols','lin','minx','miny','cellsize','empvalue');
    subplot (2,ceil(length(p)/2),i)
    imagesc(D)
    title(['PLAND=' num2str(1-p(i),1)])
    axis off
    axis square
    colormap([1 1 1;0 0 0])
end
mean(mean(D))
