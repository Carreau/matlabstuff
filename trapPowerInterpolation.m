classdef (ConstructOnLoad) trapPowerInterpolation < handle
% trapPowerInterpolation
%
% Object to deal with power interpolation
% create it with the calibration like this
% 
% > obj = trapPowerInterpolation(X,Y,Z)
%       X is the list of X postion
%       Y is the list of Y postion
%           X,Y created as meshgrid would do 
%       Z the efalution of the fonction you wish to interpolate between X and Y
%
% > obj.interpolatedPower(x,y)
%       will return the '*linear' power interpolation at x,y
% > obj.interpolatedPower % no argument
%       will return an fontion handle which doeas the same
%       
% %Example 
% [X,Y]=meshgrid(-1:0.5:1,-1:0.5:1);
% [XX,YY]=meshgrid(-1:0.1:1,-1:0.1:1);
% 
% Z=cos(X).*sin(Y);
% ZZ=cos(XX).*sin(YY);
% 
% tpi=trapPowerInterpolation(X,Y,Z);
% 
% surf(XX,YY,tpi.interpolatedPower(XX,YY))
% %same as 
% f=tpi.interpolatedPower;
% surf(XX,YY,f(XX,YY))
%
    properties
        XMeshGrid;
        YMeshGrid;
        ZMeshGrid;
    end
    methods
        function p=interpolatedPower(self,x,y)
            if nargin==1
            p = @(x,y)interp2(self.XMeshGrid,self.YMeshGrid,self.ZMeshGrid,x,y,'*linear');
            else
            p = interp2(self.XMeshGrid,self.YMeshGrid,self.ZMeshGrid,x,y,'*linear');
            end
        end
        function obj=trapPowerInterpolation(X,Y,Z)
            obj.XMeshGrid=X;
            obj.YMeshGrid=Y;
            obj.ZMeshGrid=Z;
        end
    end
end