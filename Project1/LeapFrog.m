function [ u , numberOfConvergeTestIteration] = LeapFrog( dx, dt )
%INPUT
%dt: Discretization size of time domain
%dx: Discretization size of spatial domain

%OUTPUT
%u: Array of values of wave equation at each grid point

TOL = 10^(-5); %Tolerance of discrepancy between predicted value and corrected value
x=0;    %Spatial domain boundary (start)
X=40;   %Spatial domain boundary (end)
t=0;    %Time domain boundary (start)
T=10;   %Time domain boundary (end)
C=1/2;  %C = wave speed
v=C*dt/dx;    

numberOfConvergeTestIteration=0;    %Initialize number of iteration of converge test

% Initialize the exact Dirichlet boundary conditions
for n=1:((T/dt)+1)
    u(1,n)=0;
    u(((X/dx)+1),n)=1;
end

% Initialize the initial condition
for j=1:((X/dx)+1)
    u(j,1)=( 1/2)*(1+tanh(250*(x-20)));
    x=x+dx;
end


%Since leap-frog scheme is a two-level scheme, requiring two previous time
%steps,'Predictor/Corrector' method is used to obtain values at time step 1
%from the values at time step 0

%Initialize the boundary contitions
halfStepU(1)=0;
tempU(1)=0;
uPredicted(1)=0;
halfStepU(((X/dx)+1))=1;
tempU(((X/dx)+1))=1;
uPredicted(((X/dx)+1))=1;

%Values at time step 1 is predicted using 'Upwind' scheme (Predictor)
for j=2:((X/dx)+1)
    uPredicted(j)=(1-v)*u(j,1)+v*(u(j-1,1));
end

%Then, values at time step 1/2 is estimated by averaging step 0 and
%predicted step 1 to be used in corrector method, which is Leap Frog scheme

%infinite loop
while 1
    
    numberOfConvergeTestIteration = numberOfConvergeTestIteration + 1;  
    %Count how many iteration has been made to converge the corrected values
    convergeTest = 0;       %Initialize converge test index every iteration
    
    %Find values at time step 1/2 by averaging step 0 and predicted step 1
    for j=2:(X/dx)
        halfStepU(j)=(1/2)*(uPredicted(j)+u(j,1));
    end

    %Store Corrected values, obtained by Leap Frog method, in temporary array
    for j=2:(X/dx)
        tempU(j)=u(j,1)-v*(halfStepU(j+1)-halfStepU(j-1));
    end
    
    %Test if the discrepancy of corrected value and predicted value is
    %lower than the tolerance for every element
    for j=1:((X/dx)+1)
        if (abs(tempU(j)-uPredicted(j))<TOL)
            convergeTest = convergeTest + 1;    %If every elements are satisfied the tolerance,
        end                                     %Converge Test index will be number of
    end                                         %elements in spatial domain array
    
    %If every elements are satisfied, store corrected value to solution array, u
    if (convergeTest==((X/dx)+1))
        for j=1:((X/dx)+1)
            u(j,2)=tempU(j);
        end
        break       %Then, exit the infinite loop
    
    %If not, store corrected value in predicted value array, and then
    %repeat the iteration until it converges
    else
        for j=1:((X/dx)+1)
            uPredicted(j)=tempU(j);
        end
    end 
end

%Leap-Frog Method
for n=2:(T/dt)
    for j=2:(X/dx)
        u(j,n+1)=u(j,n-1)-v*(u(j+1,n)-u(j-1,n));
    end
end

end
