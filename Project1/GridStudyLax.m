function [ u ] = GridStudyLax( dx,dt,t )
%INPUT
%dt: Discretization size of time domain
%dx: Discretization size of spatial domain
%t: Specific time frame at which the values will be plotted (0? t ?10)

%OUTPUT
%u: Array of values of wave equation

%Compute the wave equation using Lax scheme
u = Lax(dx,dt);

T=10;   %Time domain boundary (end)
X=40;   %Spatial domain boundary (end)

dn = ((T/dt))/10;   %Calculate the discrepancy between each time fame
n = 1+dn*t;         %Calculate appropriate index that corresponds to specified time frame

hold on,
xlabel('x')
ylabel('u(x)')
y=linspace(0,40,((X/dx)+1));
plot(y,u(:,n),'-bo');       %Plot the values at specified time frame
hold off

end
