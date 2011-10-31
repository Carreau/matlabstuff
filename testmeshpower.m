%%
clear

[X,Y]=meshgrid(-1:0.5:1,-1:0.5:1);
[XX,YY]=meshgrid(-1:0.1:1,-1:0.1:1);

Z=cos(X).*sin(Y);
ZZ=cos(XX).*sin(YY);

tpi=trapPowerInterpolation(X,Y,Z);

surf(XX,YY,tpi.interpolatedPower(XX,YY))
%same as 
f=tpi.interpolatedPower;
surf(XX,YY,f(XX,YY))



