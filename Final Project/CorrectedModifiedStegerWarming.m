function [ FluxRight ] = CorrectedModifiedStegerWarming( Wcurrent,Wnext,epsilon )
%Justin ChanWoo Yang
%260368098
%
%INUPUT
%Wcurrent: W matrix of grid at i
%Wnext: W matrix of grid at i+1
%
%OUTPUT
%FluxRight: Flux matrix at i+1/2

gamma = 1.4;

rhoCurrent = Wcurrent(1);
uCurrent = Wcurrent(2)/rhoCurrent;
Pcurrent = (gamma-1)*(Wcurrent(3)-(1/2)*rhoCurrent*(uCurrent^2));

rhoNext = Wnext(1);
uNext = Wnext(2)/rhoNext;
Pnext = (gamma-1)*(Wnext(3)-(1/2)*rhoNext*(uNext^2));

dPdx = (Pnext - Pcurrent)/min(Pnext,Pcurrent);
omegaAvg = 1/(1+dPdx^2);

FluxRight = omegaAvg*ModifiedStegerWarming(Wcurrent',Wnext',epsilon)+(1-omegaAvg)*StegerWarming(Wcurrent',Wnext',epsilon);

end