function [ W,P,PtotalExit,M,densityResidual,error,x,dx,time,numberOfIteration ] = EulerExplicitScheme( gridNumber,PRatio,CFL,epsilon,scheme,option )
%Justin ChanWoo Yang
%260368098
%
%INUPUT
%gridNumber: Number of grids to discritize the channel
%epsilon: value of coefficient for the Scalar Dissipation scheme
%option: Choose between 1-Scalar Dissipation and 2-Steger Warming
%
%OUTPUT
%W: W matrix [rho rho*u e]
%P: P matrix that stores local static pressure across the channel
%PtotalExit: Matrix that stores total pressure across the channel
%M: matrix that stores local Mach number across the channel
%densityResidual: matrix that stores density residual at each iterations
%error: Total Pressure loss in % across the channel
%x: Location of the channel at each grids
%ds: Grid size
%time: time taken at each iterations
%numberOfIteration: total number of iterations taken to converge to the solution

Tstart = tic;

%   Variables
Tolerance = eps('single');
maxIteration = 2*10^5;
numberOfIteration = 0;

%   Grid Setting
dx = 1/gridNumber;
for i = 1:gridNumber
    x(i) = dx/2+dx*(i-1);
end

%   Nozzle Parameters
h = 0.15;
t1 = 0.8;
t2 = 3;
S = @(x) 1-h*(sin(pi*(x^t1)))^t2 ;

%   Flow Condition
gamma = 1.4;        %   Specific heat ratio
Tt = 531.2;    %   Inlet Total Temperature [R]
Pt = 2117.0;   %   Inlet Total Pressure [lb/ft^2]
R = 1716;           %   Gas Constant [ft*lb/slug*R]
Cv = R/(gamma-1);   %   Specific heat coefficient

%   Allocating Memory for State Matrix, Pressure vector and Mach Vector
W = zeros(gridNumber,3);
RESIDU = zeros(gridNumber,3);
P = zeros(gridNumber,1);
M = zeros(gridNumber,1);
u = zeros(gridNumber,1);
c = zeros(gridNumber,1);
rho = zeros(gridNumber,1);

%   Initialize intlet flow
Pexit = PRatio*Pt;

M(:) = sqrt(2/(gamma-1)*((Pexit/Pt)^(-(gamma-1)/gamma) -1));
Minlet = M(1);

Pressure = Pt*(1+((gamma-1)/2)*(M(1)^2))^(-gamma/(gamma-1));
P(:) = Pressure;
Tinlet = Tt*(1+((gamma-1)/2)*(M(1)^2))^(-1);
rhoInlet = Pressure/(R*Tinlet);
rho(:) = rhoInlet ;
cInlet = sqrt(gamma*Pressure/rhoInlet);
c(:) = cInlet ;
uInlet = Minlet*cInlet;
u(:) = uInlet ;
eInlet = Pressure/(gamma-1)+(1/2)*rhoInlet*(uInlet^2);
a = sqrt(2*gamma*(gamma-1)/(gamma+1)*Cv*Tt);

for i = 1:gridNumber
    W(i,:) = [rhoInlet rhoInlet*uInlet eInlet];
end

%Select Temporal Discretization method
if scheme == 1
    Cycle = 1 ; %Euler Explicit
elseif scheme == 2
    Cycle = 4 ; %Jameson's Modified Runge-Kutta
end

while (numberOfIteration < maxIteration)
    W_storage = W;
    
    %   -----  SELECT SPATIAL DISCRETIZATION SCHEME  -----  %
    
    for k = 1:Cycle
        for i=2:gridNumber-1

            SLeft = (1/2)*(S(dx*(i-1/2))+S(dx*(i-3/2)));
            SRight = (1/2)*(S(dx*(i-1/2))+S(dx*(i+1/2)));
            Q = [0 P(i)*(SRight-SLeft) 0] ;
            V = (1/2)*(SRight+SLeft)*dx;

            rhoCurrent = W(i,1);
            uCurrent = W(i,2)/rhoCurrent;
            cCurrent = sqrt(gamma*P(i)/rhoCurrent);

 
            if scheme == 1              %Euler Explicit
                dt = CFL*dx/uCurrent ;  
            elseif scheme == 2          %Modified 4th Order Runge Kutta
                dt = CFL*((2*sqrt(2))/(abs(uCurrent)/dx+cCurrent*sqrt(1/(dx^2))));
            end
            
            %Select Flux Evaluation Method
            if option == 1
                FluxRight = ScalarDissipation(W(i,:),W(i+1,:),epsilon);
                FluxLeft = ScalarDissipation(W(i-1,:),W(i,:),epsilon);
            elseif option == 2
                FluxRight = StegerWarming(W(i,:)',W(i+1,:)',epsilon);
                FluxLeft = StegerWarming(W(i-1,:)',W(i,:)',epsilon);
            elseif option == 3
                FluxRight = CorrectedModifiedStegerWarming(W(i,:),W(i+1,:),epsilon);
                FluxLeft = CorrectedModifiedStegerWarming(W(i-1,:),W(i,:),epsilon);
            elseif option == 4
                FluxRight = Roe(W(i,:),W(i+1,:));
                FluxLeft = Roe(W(i-1,:),W(i,:));
            end
            
            alpha = 1/(Cycle+1-k);
            RESIDU(i,:) = (FluxRight*SRight-FluxLeft*SLeft)-Q ;
            W(i,:) = W_storage(i,:)-alpha*(dt/V)*RESIDU(i,:);
        end
    end
    %   ---------   Updating Flow Parameters   ----------   %
    
    for i = 2:gridNumber-1
        P(i) = (gamma-1)*(W(i,3)-(1/2)*(W(i,2)^2)./W(i,1));
        rho(i) = W(i,1);
        u(i) = W(i,2)/rho(i);
        c(i) = sqrt((gamma*P(i))/rho(i));
        M(i) = u(i)/c(i);
    end
    
    
    %Characteristic Boundary Conditions

    %-------------Inlet Boundary---------------
    if Minlet<1
        dp_du = Pt*(gamma/(gamma-1))*(1-(gamma-1)/(gamma+1)*u(1)^2/a^2)^(1/(gamma-1))*(-2*(gamma-1)/(gamma+1)*u(1)/a^2);

        dt_inlet = CFL*dx/(u(1)+c(1));
        Lambda = ((u(2)+u(1))/2 - (c(2)+c(1))/2)*(dt_inlet/dx);
        du_inlet = -Lambda*( P(2)-P(1) - rho(1)*c(1)*(u(2)-u(1)) ) / (dp_du - rho(1)*c(1));
        
        %Update Flow Properties
        u(1) = u(1) + du_inlet;
        T_inlet = Tt*(1-(gamma-1)/(gamma+1)*u(1)^2/a^2);
        P(1) = Pt*(T_inlet/Tt)^(gamma/(gamma-1));
        rho(1) = P(1)/(R*T_inlet);
        Energy_inlet = rho(1)*(Cv*T_inlet+1/2*u(1)^2);
        c(1) = sqrt(gamma*P(1)/rho(1));
        M(1) = u(1)/c(1);
        
        W(1,:) = [rho(1) rho(1)*u(1) Energy_inlet];
    end
    
    %-------------Exit Boundary---------------
    dt_exit = CFL*dx/(u(gridNumber)+c(gridNumber));
    
    %Compute eigenvalues
    Lambda_1 = (u(gridNumber)+u(gridNumber-1))/2 * dt_exit/dx;
    Lambda_2 = ((u(gridNumber)+u(gridNumber-1))/2 + (c(gridNumber)+c(gridNumber-1))/2) * dt_exit/dx;
    Lambda_3 = ((u(gridNumber)+u(gridNumber-1))/2 - (c(gridNumber)+c(gridNumber-1))/2) * dt_exit/dx;
    
    %Comput exit Mach number
    R1 = -Lambda_1*(rho(gridNumber)-rho(gridNumber-1) - 1/c(gridNumber)^2*(P(gridNumber)-P(gridNumber-1)));
    R2 = -Lambda_2*(P(gridNumber)-P(gridNumber-1) + rho(gridNumber)*c(gridNumber)*(u(gridNumber)-u(gridNumber-1)));
    R3 = -Lambda_3*(P(gridNumber)-P(gridNumber-1) - rho(gridNumber)*c(gridNumber)*(u(gridNumber)-u(gridNumber-1)));
    
    %Compute dP based on either a subsonic or supersonic exit
    M(gridNumber) = ((u(gridNumber)+u(gridNumber-1))/2) / ((c(gridNumber)+c(gridNumber-1))/2);
    
    %Compute dP based on either a subsonic or supersonic exit
    if M(gridNumber)>1
        d_p = (R2+R3)/2;
    else
        d_p = 0;
    end
    d_rho = R1 + 1/c(gridNumber)^2*(d_p);
    d_u = (R2 - d_p)/(rho(gridNumber)*c(gridNumber));
    
    %Update flow properties
    rho(gridNumber) = rho(gridNumber) + d_rho;
    u(gridNumber) = u(gridNumber) + d_u;
    P(gridNumber) = P(gridNumber) + d_p;
    T_exit = P(gridNumber)/(R*rho(gridNumber));
    Energy_exit = rho(gridNumber)*(Cv*T_exit+1/2*u(gridNumber)^2);
    c(gridNumber) = sqrt(gamma*P(gridNumber)/rho(gridNumber));
    M(gridNumber) = u(gridNumber)/c(gridNumber);
    
    W(gridNumber,:) = [rho(gridNumber) rho(gridNumber)*u(gridNumber) Energy_exit];
    %---------------------------------------------------
    
    numberOfIteration = numberOfIteration + 1;
%     numberOfIterations(numberOfIteration) = numberOfIteration;
    time(numberOfIteration) = toc(Tstart);
    densityResidual(numberOfIteration) = norm(W(1,:)-W_storage(1,:));
    P_ratio = P./Pt;
    
    Convergence = norm(RESIDU);
    if (Convergence < Tolerance || numberOfIteration == maxIteration)
        break;
    end
    
    fprintf('%.17f\n',Convergence);
end

for i=1:gridNumber
    PtotalExit(i) = P(i)/((1+((gamma-1)/2)*(M(i)^2))^(-gamma/(gamma-1)));
    error(i) = 100*(PtotalExit(i)-Pt)/Pt;
end

%   --------------------------    THE END     -------------------------
%
% figure(1)
% plot(x,P)
% title('Pressure Distribution Along the Channel')
% xlabel('Channel location')
% ylabel('Pressure [lb/ft^2]')
%
% figure(2)
% plot(x,PtotalExit)
% title('Total Pressure Distribution Along the Channel')
% xlabel('Channel location')
% ylabel('Total Pressure [lb/ft^2]')
%
% figure(3)
% plot(x,M)
% title('Mach Number Distribution Along the Channel')
% xlabel('Channel location')
% ylabel('Mach Number')
%
% figure(4)
% plot(densityResidual)
% title('Residual Convergence')
% xlabel('Iteration Number')
% ylabel('Residual')
%
% figure(5)
% plot(x,error)
% title('Pressure Loss Across the Channel')
% xlabel('Channel location')
% ylabel('Pressure Loss [%]')

end