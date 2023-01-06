function [R,IG_Step] = forage_cellsDEB(R,P1,I,DMD,IGMax,ICL,Rmax,Vmax,hmax,Smax,amax)
%Function to caculate dispersal movement among cells (in m) and energy
%consumed in the process


%Foraging steps within 1h window
IG_Step = 0;

if size(P1,1) == 1 
    
   %Individual does not leave cell
   
	%Current position
    CellX = P1(1,1);
    CellY = P1(1,2);

    %New position
	NewX = CellX;
    NewY = CellY;

    %Density within each cell
    E_un = 10000; %10000J for each food unit (1g dry mass of C).
    
    %Current Cell
    if R(CellY,CellX) ==0
        DCell=0;
        VCell=0;
        Tcell=0;
        IG_Cell=0;
    else
        DCell = 1/sqrt(R(CellY,CellX)/E_un/1000);%Units of food/cell area (1000 m²).
        %Velocity within cell (Shipley 1996, eq 11)
        VCell = (Vmax*DCell/(DCell+(Vmax^2/amax)))*60;%(m/s --> m/min)
        %Time within cells (30min in each cell)
        Tcell = (60/size(P1,1))/2;
        %Energy ingested in the current cell
        if DCell>(VCell*Smax/Rmax)
            Denscell = sqrt(R(CellY,CellX)/E_un/1000);
            IG_Cell = (VCell*Denscell*Smax)/(1+hmax*VCell*Denscell);%g/min
            IG_Cell = IG_Cell*10*1000;%J/min
            IG_Cell = IG_Cell * Tcell;%J/cell
        else
            IG_Cell = ((Rmax*Smax)/(Rmax*hmax+Smax))*10*1000;%J/min
            IG_Cell = IG_Cell * Tcell;%J/cell
        end
    end
    
    %New Cell
    if R(NewY,NewX)==0
        DNew=0;
        VNew=0;
        TNewcell=0;
        IG_New=0;
    else
        DNew = 1/sqrt(R(NewY,NewX)/E_un/1000);
        %Velocity within cell (Shipley 1996, eq 11)
        VNew = (Vmax*DNew/(DNew+(Vmax^2/amax)))*60;%(m/s --> m/min)
        TNewcell=(60/size(P1,1))/2;
        %Energy ingested in the new cell
        if DNew>(VNew*Smax/Rmax)
            DensNew = sqrt(R(NewY,NewX)/E_un/1000);
            IG_New = (VNew*DensNew*Smax)/(1+hmax*VNew*DensNew);%g/min
            IG_New = IG_New*10*1000;%J/min
            IG_New = IG_New * Tcell;%J/cell
        else
            IG_New = ((Rmax*Smax)/(Rmax*hmax+Smax))*10*1000;%J/min
            IG_New = IG_New * TNewcell;%J/cell
        end
    end

    %Velocity within cell (Shipley 1996, eq 11)
%     VCell = (Vmax*DCell/(DCell+(Vmax^2/amax)))*60;%(m/s --> m/min)
%     VNew = (Vmax*DNew/(DNew+(Vmax^2/amax)))*60;%(m/s --> m/min)

    %Time within cells (30min in each cell)
%     Tcell = (60/size(P1,1))/2;
%     TNewcell = (60/size(P1,1))/2;

    %Energy ingested in the current cell
%     if DCell>(VCell*Smax/Rmax)
%         Denscell = sqrt(R(CellY,CellX)/E_un/1000);
%         IG_Cell = (VCell*Denscell*Smax)/(1+hmax*VCell*Denscell);%g/min
%         IG_Cell = IG_Cell*10*1000;%J/min
%         IG_Cell = IG_Cell * Tcell;%J/cell
%     else
%         IG_Cell = ((Rmax*Smax)/(Rmax*hmax+Smax))*10*1000;%J/min
%         IG_Cell = IG_Cell * Tcell;%J/cell
%     end
    
    %Correct for overhead feeding within a single cell
    if IG_Cell > IGMax/2
        IG_Cell = IGMax/2;
    end

    %Correct consumption overhead
    if IG_Cell > R(CellY,CellX)
        IG_Cell = R(CellY,CellX);
    end

    %Discount for locomotion costs
    Gasto_Cell=ICL*(DMD/2);%move half distance in each cell
    LCell = IG_Cell-Gasto_Cell;

    %Energy ingested in the new cell
%     if DNew>(VNew*Smax/Rmax)
%         DensNew = sqrt(R(NewY,NewX)/E_un/1000);
%         IG_New = (VNew*DensNew*Smax)/(1+hmax*VNew*DensNew);%g/min
%         IG_New = IG_New*10*1000;%J/min
%         IG_New = IG_New * Tcell;%J/cell
%     else
%         IG_New = ((Rmax*Smax)/(Rmax*hmax+Smax))*10*1000;%J/min
%         IG_New = IG_New * TNewcell;%J/cell
%     end
    
    %Correct for overhead feeding within a single cell
    if IG_New > IGMax/2
        IG_New = IGMax/2;
    end

    %Correct consumption overhead
    if IG_New > R(NewY,NewX)
        IG_New = R(NewY,NewX);
    end

    %Discount for locomotion costs
    Gasto_New=ICL*(DMD/2);%move half distance in each cell
    LNew = IG_New-Gasto_New;

    %Maximum Consumption
%     IG_MAX_STEP = (((Rmax*Smax)/(Rmax*hmax+Smax))*10*1000)*60;
%     IG_MAX_STEP = IG_MAX_STEP-(ICL*DMDt);
%     F_food = (LCell+LNew)/IG_MAX_STEP;

    %Remove energy from the system
    R(CellY,CellX) = R(CellY,CellX)-IG_Cell;
    if R(CellY,CellX)<0
        R(CellY,CellX)=0;
    end
    R(NewY,NewX) = R(NewY,NewX)-IG_New;
    if R(NewY,NewX)<0
        R(NewY,NewX)=0;
    end
    
    %Ingested energy in step within 1h
    IG_Step = IG_Step+(LCell+LNew);
    
else

    %Individual leaves cell
    
    for i = 1:(size(P1,1)-1)
        %Current position
        CellX = P1(i,1);
        CellY = P1(i,2);

        %New position
        if size(P1,1)>1
            NewX = P1(i+1,1);
            NewY = P1(1+1,2);
        else
            NewX = P1(i,1);
            NewY = P1(i,2);
        end

        %Density within each cell
        E_un = 10000; %10000J for each food unit (1g dry mass of C).
        if R(CellY,CellX) ==0
            DCell=0;
            VCell=0;
        else
            DCell = 1/sqrt(R(CellY,CellX)/E_un/1000);%Units of food/cell area (1000 m²).
            %Velocity within cell (Shipley 1996, eq 11)
            VCell = (Vmax*DCell/(DCell+(Vmax^2/amax)))*60;%(m/s --> m/min)
        end

        if R(NewY,NewX)==0
            DNew=0;
            VNew=0;
        else
            DNew = 1/sqrt(R(NewY,NewX)/E_un/1000);
            %Velocity within cell (Shipley 1996, eq 11)
            VNew = (Vmax*DNew/(DNew+(Vmax^2/amax)))*60;%(m/s --> m/min)
        end
%         DCell = 1/sqrt(R(CellY,CellX)/E_un/1000);%Units of food/cell area (1000 m²).
%         DNew = 1/sqrt(R(NewY,NewX)/E_un/1000);
% 
%         %Velocity within cell (Shipley 1996, eq 11)
%         VCell = (Vmax*DCell/(DCell+(Vmax^2/amax)))*60;%(m/s --> m/min)
%         VNew = (Vmax*DNew/(DNew+(Vmax^2/amax)))*60;%(m/s --> m/min)

        %Time within cells (30min in each cell)
        Tcell = (60/size(P1,1))/2;
        TNewcell = (60/size(P1,1))/2;

        %Energy ingested in the current cell
        if DCell>(VCell*Smax/Rmax)
            Denscell = sqrt(R(CellY,CellX)/E_un/1000);
            IG_Cell = (VCell*Denscell*Smax)/(1+hmax*VCell*Denscell);%g/min
            IG_Cell = IG_Cell*10*1000;%J/min
            IG_Cell = IG_Cell * Tcell;%J/cell
        else
            IG_Cell = ((Rmax*Smax)/(Rmax*hmax+Smax))*10*1000;%J/min
            IG_Cell = IG_Cell * Tcell;%J/cell
        end

        %Correct consumption overhead
        if IG_Cell > R(CellY,CellX)
            IG_Cell = R(CellY,CellX);
        end

        %Discount for locomotion costs
        Gasto_Cell=ICL*(DMD/2);%move half distance in each cell
        LCell = IG_Cell-Gasto_Cell;

        %Energy ingested in the new cell
        if DNew>(VNew*Smax/Rmax)
            DensNew = sqrt(R(NewY,NewX)/E_un/1000);
            IG_New = (VNew*DensNew*Smax)/(1+hmax*VNew*DensNew);%g/min
            IG_New = IG_New*10*1000;%J/min
            IG_New = IG_New * Tcell;%J/cell
        else
            IG_New = ((Rmax*Smax)/(Rmax*hmax+Smax))*10*1000;%J/min
            IG_New = IG_New * TNewcell;%J/cell
        end

        %Correct consumption overhead
        if IG_New > R(NewY,NewX)
            IG_New = R(NewY,NewX);
        end

        %Discount for locomotion costs
        Gasto_New=ICL*(DMD/2);%move half distance in each cell
        LNew = IG_New-Gasto_New;

        %Maximum Consumption
    %     IG_MAX_STEP = (((Rmax*Smax)/(Rmax*hmax+Smax))*10*1000)*60;
    %     IG_MAX_STEP = IG_MAX_STEP-(ICL*DMDt);
    %     F_food = (LCell+LNew)/IG_MAX_STEP;

        %Remove energy from the system
        R(CellY,CellX) = R(CellY,CellX)-IG_Cell;
        if R(CellY,CellX)<0
            R(CellY,CellX)=0;
        end
        R(NewY,NewX) = R(NewY,NewX)-IG_New;
        if R(NewY,NewX)<0
            R(NewY,NewX)=0;
        end

        %Ingested energy in step within 1h
        IG_Step = IG_Step+(LCell+LNew);
    end
end