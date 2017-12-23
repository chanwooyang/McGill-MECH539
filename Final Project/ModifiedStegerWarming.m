function [ FluxRight ] = ModifiedStegerWarming( Wcurrent,Wnext,epsilon )
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

Waverage = (Wcurrent+Wnext)/2;

rhoAverage = Waverage(1);
uAverage = Waverage(2)/rhoAverage;
Paverage = (gamma-1)*(Waverage(3)-(1/2)*rhoAverage*(uAverage^2));
cAverage = sqrt(gamma*Paverage/Waverage(1));

alphaAverage = (1/2)*(uAverage^2);
betaAverage = gamma-1;

Smatrix = [1 0 0; -uAverage/rhoAverage 1/rhoAverage 0; alphaAverage*betaAverage -uAverage*betaAverage betaAverage];
S1matrix = [1 0 0; uAverage rhoAverage 0; alphaAverage rhoAverage*uAverage 1/betaAverage];
positiveDmatrix = [(uAverage+sqrt(uAverage^2+epsilon^2))/2 0 0; 0 ((uAverage+cAverage)+...
    sqrt((uAverage+cAverage)^2+epsilon^2))/2 0; 0 0 ((uAverage-cAverage)+sqrt((uAverage-cAverage)^2+epsilon^2))/2];
negativeDmatrix = [(uAverage-sqrt(uAverage^2+epsilon^2))/2 0 0; 0 ((uAverage+cAverage)-...
    sqrt((uAverage+cAverage)^2+epsilon^2))/2 0; 0 0 ((uAverage-cAverage)-sqrt((uAverage-cAverage)^2+epsilon^2))/2];
Cmatrix = [1 0 -1/(cAverage^2); 0 rhoAverage*cAverage 1; 0 -rhoAverage*cAverage 1];
C1matrix = [1 1/(2*(cAverage^2)) 1/(2*(cAverage^2)); 0 1/(2*rhoAverage*cAverage)...
    -1/(2*rhoAverage*cAverage); 0 1/2 1/2];


positiveAmatrix = S1matrix*C1matrix*positiveDmatrix*Cmatrix*Smatrix;

negativeAmatrix = S1matrix*C1matrix*negativeDmatrix*Cmatrix*Smatrix;

Fluxright = positiveAmatrix*Wcurrent+negativeAmatrix*Wnext;
FluxRight = Fluxright';
end