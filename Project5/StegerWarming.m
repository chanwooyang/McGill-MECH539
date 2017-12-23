function [ FluxRight ] = StegerWarming( Wcurrent,Wnext )
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
cCurrent = sqrt(gamma*Pcurrent/Wcurrent(1));

alphaCurrent = (1/2)*(uCurrent^2);
betaCurrent = gamma-1;

Smatrix = [1 0 0; -uCurrent/rhoCurrent 1/rhoCurrent 0; alphaCurrent*betaCurrent -uCurrent*betaCurrent betaCurrent];
S1matrix = [1 0 0; uCurrent rhoCurrent 0; alphaCurrent rhoCurrent*uCurrent 1/betaCurrent];
positiveDmatrix = [(uCurrent+sqrt(uCurrent^2))/2 0 0; 0 ((uCurrent+cCurrent)+...
    sqrt((uCurrent+cCurrent)^2))/2 0; 0 0 ((uCurrent-cCurrent)+sqrt((uCurrent-cCurrent)^2))/2];
Cmatrix = [1 0 -1/(cCurrent^2); 0 rhoCurrent*cCurrent 1; 0 -rhoCurrent*cCurrent 1];
C1matrix = [1 1/(2*(cCurrent^2)) 1/(2*(cCurrent^2)); 0 1/(2*rhoCurrent*cCurrent)...
    -1/(2*rhoCurrent*cCurrent); 0 1/2 1/2];

positiveAmatrix = S1matrix*C1matrix*positiveDmatrix*Cmatrix*Smatrix;



rhoNext = Wnext(1);
uNext = Wnext(2)/rhoNext;
Pnext = (gamma-1)*(Wnext(3)-(1/2)*rhoNext*(uNext^2));
cNext = sqrt(gamma*Pnext/Wnext(1));

alpha = (1/2)*(uNext^2);
beta = gamma-1;

Smatrix = [1 0 0; -uNext/rhoNext 1/rhoNext 0; alpha*beta -uNext*beta beta];
S1matrix = [1 0 0; uNext rhoNext 0; alpha rhoNext*uNext 1/beta];
negativeDmatrix = [(uNext-sqrt(uNext^2))/2 0 0; 0 ((uNext+cNext)-...
    sqrt((uNext+cNext)^2))/2 0; 0 0 ((uNext-cNext)-sqrt((uNext-cNext)^2))/2];
Cmatrix = [1 0 -1/(cNext^2); 0 rhoNext*cNext 1; 0 -rhoNext*cNext 1];
C1matrix = [1 1/(2*(cNext^2)) 1/(2*(cNext^2)); 0 1/(2*rhoNext*cNext)...
    -1/(2*rhoNext*cNext); 0 1/2 1/2];

negativeAmatrix = S1matrix*C1matrix*negativeDmatrix*Cmatrix*Smatrix;

Fluxright = positiveAmatrix*Wcurrent+negativeAmatrix*Wnext;
FluxRight = Fluxright';
end