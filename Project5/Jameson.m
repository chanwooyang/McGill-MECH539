function [ W,P,PtotalExit,M,densityResidual,error,x,dx,time,numberOfIteration ] = Jameson( gridNumber, epsilon, option )
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

%Variable
Tolerance = eps;
maxIteration = 10^5;
numberOfIteration = 0;
dx = 1/gridNumber;
for i = 1:gridNumber
    x(i) = dx/2+dx*(i-1);
end
CFL = 0.95;

%Nozzle Parameter
h = 0.025;  %Bump height
t1 = 0.8;   %Maximum location of the bump
t2 = 3;     %Controls the width of the bump
S = @(x) 1-h*((sin(pi*(x^t1)))^t2);

%Flow Conditions
gamma = 1.4;    %Specific heat ratio
Tt = 531.2;     %Inlet Total Temperature [R]
Pt = 2117.0;    %Inlet Total Pressure [lb/ft^2]
R = 1716;       %Gas Constant [ft*lb/slug*R]
Minlet = 1.2;   %Inlet Mach Number

%Declare variables
W = zeros(gridNumber,3);
W1 = zeros(gridNumber,3);
W2 = zeros(gridNumber,3);
W3 = zeros(gridNumber,3);
P = zeros(gridNumber,1);      %Pressure

%Initialize intlet flow
Pinlet = Pt*(1+((gamma-1)/2)*(Minlet^2))^(-gamma/(gamma-1));
Tinlet = Tt*(1+((gamma-1)/2)*(Minlet^2))^(-1);
P(:) = Pinlet;
M(1) = Minlet;
rhoInlet = Pinlet/(R*Tinlet);
cInlet = sqrt(gamma*Pinlet/rhoInlet);
uInlet = Minlet*cInlet;
eInlet = Pinlet/(gamma-1)+(1/2)*rhoInlet*(uInlet^2);

for i=1:gridNumber
    W(i,:) = [rhoInlet rhoInlet*uInlet eInlet];
    W1(i,:) = [rhoInlet rhoInlet*uInlet eInlet];
    W2(i,:) = [rhoInlet rhoInlet*uInlet eInlet];
    W3(i,:) = [rhoInlet rhoInlet*uInlet eInlet];
end

previousW = W;

while (numberOfIteration < maxIteration)
    
    %Jameson's fourth order Runge Kutta
    for j=1:4
        for i=2:gridNumber-1
            rho = W(i,1);
            u = W(i,2)/rho;
            c = sqrt(gamma*P(i)/rho);
            dt = CFL*((2*sqrt(2))/(abs(u)/dx+c*sqrt(1/(dx^2))));
            
            SLeft = (1/2)*(S(dx/2+dx*(i-1))+S(dx/2+dx*(i-2)));
            SRight = (1/2)*(S(dx/2+dx*(i-1))+S(dx/2+dx*i));
            V = (1/2)*(SRight+SLeft)*dx;
            Q = [0 P(i)*(SRight-SLeft) 0];
            
            if j==1
                if option==1
                    FluxRight0 = ScalarDissipation(W(i,:),W(i+1,:),epsilon);
                    FluxLeft0 = ScalarDissipation(W(i-1,:),W(i,:),epsilon);
                elseif option==2
                    FluxRight0 = StegerWarming(W(i,:)',W(i+1,:)');
                    FluxLeft0 = StegerWarming(W(i-1,:)',W(i,:)');
                end
                W1(i,:)=W(i,:)+((1/4)*dt/V)*(-(FluxRight0*SRight-FluxLeft0*SLeft)+Q);
            elseif j==2
                if option==1
                    FluxRight1 = ScalarDissipation(W1(i,:),W1(i+1,:),epsilon);
                    FluxLeft1 = ScalarDissipation(W1(i-1,:),W1(i,:),epsilon);
                elseif option==2
                    FluxRight1 = StegerWarming(W1(i,:)',W1(i+1,:)');
                    FluxLeft1 = StegerWarming(W1(i-1,:)',W1(i,:)');
                end
                W2(i,:)=W(i,:)+((1/3)*dt/V)*(-(FluxRight1*SRight-FluxLeft1*SLeft)+Q);
            elseif j==3
                if option==1
                    FluxRight2 = ScalarDissipation(W2(i,:),W2(i+1,:),epsilon);
                    FluxLeft2 = ScalarDissipation(W2(i-1,:),W2(i,:),epsilon);
                elseif option==2
                    FluxRight2 = StegerWarming(W2(i,:)',W2(i+1,:)');
                    FluxLeft2 = StegerWarming(W2(i-1,:)',W2(i,:)');
                end
                W3(i,:)=W(i,:)+((1/2)*dt/V)*(-(FluxRight2*SRight-FluxLeft2*SLeft)+Q);
            elseif j==4
                if option==1
                    FluxRight3 = ScalarDissipation(W3(i,:),W3(i+1,:),epsilon);
                    FluxLeft3 = ScalarDissipation(W3(i-1,:),W3(i,:),epsilon);
                elseif option==2
                    FluxRight3 = StegerWarming(W3(i,:)',W3(i+1,:)');
                    FluxLeft3 = StegerWarming(W3(i-1,:)',W3(i,:)');
                end
                W(i,:) = W(i,:)+(dt/V)*(-(FluxRight3*SRight-FluxLeft3*SLeft)+Q);
            end
        end
    end
    
    for i=2:gridNumber
        P(i) = (gamma-1)*(W(i,3)-(1/2)*(W(i,2)^2)/W(i,1));
        rho = W(i,1);
        u = W(i,2)/rho;
        c = sqrt((gamma*P(i))/rho);
        M(i) = u/c;
    end
    
    numberOfIteration = numberOfIteration + 1;
    densityResidual(numberOfIteration+1) = norm(abs(W(:,1)-previousW(:,1)));
        
    time(numberOfIteration+1) = toc(Tstart);
    
    if (densityResidual(numberOfIteration+1) < Tolerance || numberOfIteration==maxIteration)
        for i=1:gridNumber
            PtotalExit(i) = P(i)/((1+((gamma-1)/2)*(M(i)^2))^(-gamma/(gamma-1)));
            error(i) = 100*abs(Pt-PtotalExit(i))/Pt;
        end
        break;
    end
    
    %Exit Boundary Condition
    W(gridNumber,:) = W(gridNumber-1,:);
    
    previousW = W;

end

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

