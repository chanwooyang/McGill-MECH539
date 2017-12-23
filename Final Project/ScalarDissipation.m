function [ FluxRight ] = ScalarDissipation( Wcurrent, Wnext,epsilon )
%Justin ChanWoo Yang
%260368098
%
%INUPUT
%Wcurrent: W matrix of grid at i
%Wnext: W matrix of grid at i+1
%epsilon: value of coefficient for the Scalar Dissipation scheme
%
%OUTPUT
%FluxRight: Flux matrix at i+1/2

gamma = 1.4;

uCurrent = Wcurrent(2)/Wcurrent(1);
rhoCurrent = Wcurrent(1);
Pcurrent = (gamma-1)*(Wcurrent(3)-(1/2)*rhoCurrent*(uCurrent^2));
eCurrent = Wcurrent(3);
FluxCurrent = zeros(1,3);
FluxCurrent = [rhoCurrent*uCurrent rhoCurrent*(uCurrent^2)+Pcurrent (eCurrent+Pcurrent)*uCurrent];

uNext = Wnext(2)/Wnext(1);
rhoNext = Wnext(1);
Pnext = (gamma-1)*(Wnext(3)-(1/2)*rhoNext*(uNext^2));
eNext = Wnext(3);
FluxNext = zeros(1,3);
FluxNext = [rhoNext*uNext rhoNext*(uNext^2)+Pnext (eNext+Pnext)*uNext];

rhoMid = (1/2)*(rhoCurrent+rhoNext);
rhoUMid = (1/2)*(Wcurrent(2)+Wnext(2));
Pmid = (1/2)*(Pcurrent+Pnext);

uMid = rhoUMid/rhoMid;
cMid = sqrt((gamma*Pmid)/rhoMid);

if uMid>0
    lambda = uMid+cMid;
else
    lambda = max([uMid uMid+cMid uMid-cMid]);
end

FluxRight = (1/2)*(FluxCurrent+FluxNext)-(1/2)*epsilon*lambda*(Wnext-Wcurrent);

end

