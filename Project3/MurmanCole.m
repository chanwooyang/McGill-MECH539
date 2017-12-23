function [ p, x, y, Cp, residual, time, numberOfIteration ] = MurmanCole( Mfree )
%Justin ChanWoo Yang
%260368098

%INPUT
%Mfree: Freestream mach number

%OUTPUT
%p: Solution of potential function
%x: x value at each grid point along x-axis
%y: y value at each grid point along y-axis
%Ufree: Freestream velocity
%numberOfIteration: Number of iteration

Tstart = cputime;

Tol = eps('single');
maxIteration = 10^20;
numberOfIteration = 0;

gamma = 1.4;                            %Specific heat ratio for air
TR = 0.08;                              %Thickness ratio
Tfree = 273.15;                         %Freestream temperature in kelvin (Assumption)
R = 8314.4621;                          %Universal gas constant
MWAir = 28.97; 
densityAir = 1.2922;
Pfree = densityAir*(R/MWAir)*Tfree;                         %Molecular weight of air
aFree = sqrt(gamma*(R/MWAir)*Tfree);    %Freestream speed of sound
Ufree = Mfree*aFree;                    %FreeStream mach number

%Number of grid
GridXBeforeLE = 20;
GridXAlongAF = 20;
GridXAfterTE = 20;
numberOfGridY = 30;

%Grid setup (Polynomial stretching of the grid)
%x-Grid setup
for i=1:GridXBeforeLE-1
    x((GridXBeforeLE+1)-i)=20-((20/(GridXBeforeLE*(GridXBeforeLE-1)))*(i^2)-...
        (20/(GridXBeforeLE*(GridXBeforeLE-1)))*i);
end

dx = 1/GridXAlongAF;
for i=1:GridXAlongAF
    x(i+GridXBeforeLE)=20+dx*i;
end

for i=1:(GridXAfterTE+1)
    x(i+(GridXBeforeLE+GridXAlongAF-1))=21+((29/(GridXAfterTE*(GridXAfterTE+1)))*(i^2)-...
        (29/(GridXAfterTE*(GridXAfterTE+1)))*i);
end

%y-Grid setup
for i=GridXBeforeLE:GridXBeforeLE+GridXAlongAF
    airfoilY(i) = TR*(-2*((x(i))^2)+82*x(i)-840);
    y(1,i)=airfoilY(i);
end

for i=1:GridXBeforeLE+GridXAlongAF+GridXAfterTE
    for j=1:numberOfGridY
        if i<GridXBeforeLE || i>(GridXBeforeLE+GridXAlongAF)
            y(j,i)=(50/(numberOfGridY*(numberOfGridY-1)))*(j^2)+(-50/(numberOfGridY*(numberOfGridY-1)))*j;
        end
    end
end

for i=GridXBeforeLE:GridXBeforeLE+GridXAlongAF      %y(1) starts from the airfoil surface
    for j=2:numberOfGridY           %Polynomial stretching that first grid space of y is (Thickness Ratio)/2 (=0.04)
        y(j,i)=airfoilY(i)+(-((((50-airfoilY(i))/(numberOfGridY-numberOfGridY^2))-...
            ((50-airfoilY(i))*2/(numberOfGridY-numberOfGridY^2)+0.04)*((1+numberOfGridY)/(2-numberOfGridY)))+...
            ((50-airfoilY(i))*2/(numberOfGridY-numberOfGridY^2)+0.04)*(numberOfGridY/(2-numberOfGridY))))*(j^2)+...
            (((50-airfoilY(i))/(numberOfGridY-numberOfGridY^2))-...
            ((50-airfoilY(i))*2/(numberOfGridY-numberOfGridY^2)+0.04)*((1+numberOfGridY)/(2-numberOfGridY)))*j +...
            ((50-airfoilY(i))*2/(numberOfGridY-numberOfGridY^2)+0.04)*(numberOfGridY/(2-numberOfGridY));
        
    end
end

%Initialization
p = zeros(numberOfGridY,GridXBeforeLE+GridXAlongAF+GridXAfterTE);   %Potential function
A = zeros(numberOfGridY,GridXBeforeLE+GridXAlongAF+GridXAfterTE);
APrevious = zeros(numberOfGridY,GridXBeforeLE+GridXAlongAF+GridXAfterTE);
mu = zeros(numberOfGridY,GridXBeforeLE+GridXAlongAF+GridXAfterTE);
muPrevious = zeros(numberOfGridY,GridXBeforeLE+GridXAlongAF+GridXAfterTE);
Cp = zeros(numberOfGridY,GridXBeforeLE+GridXAlongAF+GridXAfterTE);

%Transonic Small Disturbance using Murman-Cole Method
while (numberOfIteration < maxIteration)
    
    %Call Boundary Conditions
    for i=2:GridXBeforeLE-1
        p(1,i) = p(2,i);
    end
    for i=(GridXBeforeLE+GridXAlongAF+1):(GridXBeforeLE+GridXAlongAF+GridXAfterTE-1)
        p(1,i) = p(2,i);
    end
    for i=GridXBeforeLE:(GridXBeforeLE+GridXAlongAF)
        p(1,i) = Ufree*(TR*(-4*x(i)+82))*(y(1,i)-y(2,i))+p(2,i);
    end

    
    pPrevious = p;    %Store solution of previous iteration for Gauss-Seidel and residual calculation
    
    for i=3:(GridXBeforeLE+GridXAlongAF+GridXAfterTE-1)
        for j=2:(numberOfGridY-1)
            
            A(j,i) = (1-(Mfree^2))-(gamma+1)*((Mfree^2)/Ufree)*((p(j,i+1)-p(j,i-1))/(x(i+1)-x(i-1)));    %Evaluate A
            
            if A(j,i)>0         %Evaluate mu
                mu(j,i) = 0;
            elseif A(j,i)<0
                mu(j,i) = 1;
            end
            
            %Compute Coefficients
            a=(-2*(1-muPrevious(j,i))*APrevious(j,i))/((x(i+1)-x(i-1))*(x(i+1)-x(i)))-...
                2*(1-muPrevious(j,i))*APrevious(j,i)/((x(i+1)-x(i-1))*(x(i)-x(i-1)))-...
                2/((y(j+1,i)-y(j-1,i))*(y(j+1,i)-y(j,i)))-2/((y(j+1,i)-y(j-1,i))*(y(j,i)-y(j-1,i)))+...
                2*muPrevious(j,i-1)*APrevious(j,i-1)/((x(i)-x(i-1))*(x(i)-x(i-2)));
            b=2/((y(j+1,i)-y(j-1,i))*(y(j+1,i)-y(j,i)));
            c=2/((y(j+1,i)-y(j-1,i))*(y(j,i)-y(j-1,i)));
            d=(2*(1-muPrevious(j,i))*APrevious(j,i))/((x(i+1)-x(i-1))*(x(i)-x(i-1)))-...
                (2*muPrevious(j,i-1)*APrevious(j,i-1))/((x(i)-x(i-1))*(x(i)-x(i-2)))-...
                (2*muPrevious(j,i-1)*APrevious(j,i-1))/((x(i)-x(i-2))*(x(i-1)-x(i-2)));
            e=(2*(1-muPrevious(j,i))*APrevious(j,i))/((x(i+1)-x(i-1))*(x(i+1)-x(i)));
            g=(2*muPrevious(j,i-1)*APrevious(j,i-1))/((x(i-1)-x(i-2))*(x(i)-x(i-2)));
            
            %Gauss-Seidel
            p(j,i)=(-c*p(j-1,i)-g*p(j,i-2)-d*p(j,i-1)-e*pPrevious(j,i+1)-b*pPrevious(j+1,i))/a;
            

        end
    end
    
    residual(numberOfIteration+1) = norm(abs(p-pPrevious));
    
    if residual(numberOfIteration+1)<Tol
        break;
    end
    
    numberOfIteration = numberOfIteration + 1;
    
    APrevious = A;
    muPrevious = mu;
    
    time(numberOfIteration+1) = cputime - Tstart;
 
end

for i=2:(GridXBeforeLE+GridXAlongAF+GridXAfterTE-1)
    for j = 1:numberOfGridY
        Cp(j,i) = -2*((p(j,i+1)-p(j,i-1))/(x(i+1)-x(i-1)))/Ufree;
    end
end

%Calculate local pressure and local mach number

% px = zeros(numberOfGridY,GridXBeforeLE+GridXAlongAF+GridXAfterTE);
% for j=1:numberOfGridY
%     for i=2:(GridXBeforeLE+GridXAlongAF+GridXAfterTE-1)
%         px(j,i) = (p(j,i+1)-p(j,i-1))/(x(i+1)-x(i-1));
%     end
% end
% 
% for j=1:numberOfGridY
%     for i=1:(GridXBeforeLE+GridXAlongAF+GridXAfterTE)
%         vel(j,i) = px(j,i)+Ufree;
%         a(j,i) = sqrt((aFree^2)-(gamma-1)*Ufree*px(j,i));                 %Local speed of sound
%         Pressure(j,i) = Pfree*(1-gamma*(Mfree^2)*(vel(j,i)/Ufree));       %Local Pressure
%         MLocal(j,i) = vel(j,i)/a(j,i);                                    %Local Mach number
%     end
% end

end
