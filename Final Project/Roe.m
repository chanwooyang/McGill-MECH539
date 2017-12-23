function [ FluxRight ] = Roe( Wcurrent, Wnext )
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
cCurrent = sqrt(gamma*Pcurrent/Wcurrent(1));
FluxCurrent = zeros(1,3);
FluxCurrent = [rhoCurrent*uCurrent rhoCurrent*(uCurrent^2)+Pcurrent (eCurrent+Pcurrent)*uCurrent];

uNext = Wnext(2)/Wnext(1);
rhoNext = Wnext(1);
Pnext = (gamma-1)*(Wnext(3)-(1/2)*rhoNext*(uNext^2));
eNext = Wnext(3);
cNext = sqrt(gamma*Pnext/Wnext(1));
FluxNext = zeros(1,3);
FluxNext = [rhoNext*uNext rhoNext*(uNext^2)+Pnext (eNext+Pnext)*uNext];

rhoHat = sqrt(rhoCurrent*rhoNext);
uHat = (sqrt(rhoCurrent)*uCurrent + sqrt(rhoNext)*uNext)/(sqrt(rhoCurrent)+sqrt(rhoNext));
hHat = (sqrt(rhoCurrent)*((eCurrent+rhoCurrent)/rhoCurrent) + sqrt(rhoNext)*((eNext+rhoNext)/rhoNext))/(sqrt(rhoCurrent)+sqrt(rhoNext));
cHat = sqrt((gamma-1)*(hHat-(1/2)*uHat));


alphaHat = (1/2)*(uHat^2);
betaHat = gamma-1;

Smatrix = [1 0 0; -uHat/rhoHat 1/rhoHat 0; alphaHat*betaHat -uHat*betaHat betaHat];
S1matrix = [1 0 0; uHat rhoHat 0; alphaHat rhoHat*uHat 1/betaHat];
Cmatrix = [1 0 -1/(cHat^2); 0 rhoHat*cHat 1; 0 -rhoHat*cHat 1];
C1matrix = [1 1/(2*(cHat^2)) 1/(2*(cHat^2)); 0 1/(2*rhoHat*cHat) -1/(2*rhoHat*cHat); 0 1/2 1/2];

lambdaCurrent = [uCurrent 0 0; 0 uCurrent+cCurrent 0; 0 0 uCurrent-cCurrent];
lambdaNext = [uNext 0 0; 0 uNext+cNext 0; 0 0 uNext-cNext];
lambdaHat = [uHat 0 0; 0 uHat+cHat 0; 0 0 uHat-cHat];

for i=1:3
    epsilon = max([0;lambdaHat(i,i)-lambdaCurrent(i,i);lambdaNext(i,i)-lambdaHat(i,i)]);
    if abs(lambdaHat)<=epsilon
        lambdaHat(i,i) = (1/2)*((lambdaHat^2)/epsilon+epsilon);
    end
end

positiveDmatrix = [(uHat+sqrt(uHat^2+epsilon^2))/2 0 0; 0 ((uHat+cHat)+...
    sqrt((uHat+cHat)^2+epsilon^2))/2 0; 0 0 ((uHat-cHat)+sqrt((uHat-cHat)^2+epsilon^2))/2];
negativeDmatrix = [(uHat-sqrt(uHat^2+epsilon^2))/2 0 0; 0 ((uHat+cHat)-...
    sqrt((uHat+cHat)^2+epsilon^2))/2 0; 0 0 ((uHat-cHat)-sqrt((uHat-cHat)^2+epsilon^2))/2];

positiveAmatrix = S1matrix*C1matrix*positiveDmatrix*Cmatrix*Smatrix;
negativeAmatrix = S1matrix*C1matrix*negativeDmatrix*Cmatrix*Smatrix;

AmatrixHat = positiveAmatrix-negativeAmatrix;

Fluxright = (1/2)*(FluxCurrent'+FluxNext')-(1/2)*AmatrixHat*(Wnext'-Wcurrent');
FluxRight = Fluxright';
end