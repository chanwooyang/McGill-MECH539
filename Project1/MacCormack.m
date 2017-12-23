function [ u ] = MacCormack( dx, dt )
%INPUT
%dt: Discretization size of time domain
%dx: Discretization size of spatial domain

%OUTPUT
%u: Array of values of wave equation at each grid point


x=0;    %Spatial domain boundary (start)
X=40;   %Spatial domain boundary (end)
t=0;    %Time domain boundary (start)
T=10;   %Time domain boundary (end)
C=1/2;  %C = wave speed
v=C*dt/dx;    

% Initialize the exact Dirichlet boundary conditions
for n=1:((T/dt)+1)
    uPredicted(1,n)=0;
    uPredicted(((X/dx)+1),n)=1;
    u(1,n)=0;
    u(((X/dx)+1),n)=1;
end

% Initialize the initial condition
for j=1:((X/dx)+1)
    u(j,1)=(1/2)*(1+tanh(250*(x-20)));
    x=x+dx;
end

for n=1:(T/dt)
    for j=2:(X/dx)
        uPredicted(j,n+1)=(1+v)*u(j,n)-v*u(j+1,n);      %Predictor Method
        u(j,n+1)=(1/2)*(u(j,n)+(1-v)*uPredicted(j,n+1)+v*uPredicted(j-1,n+1));  %Corrector Method
    end
end


end
