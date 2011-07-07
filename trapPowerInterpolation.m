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
            if nargin == 1
                q=X;
                %cas où l'onpasse directemetn une structures des données
                %raw enregistrées:
                stp = 2/(q.grid_size-1);

                data=q.data_r(:,1:end);
                %xtrap   = data(1,:);
                %ytrap   = data(2,:);
                %aodpower= data(3,:);
                %qpdx    = data(4,:);
                %sigmax  = data(5,:);
                %qpdy    = data(6,:);
                %sigmay  = data(7,:);
                %qpdsum  = data(8,:);
                %sigmaqpdsum          = data(9,:);
                powerphotodiode      = data(10,:);
                %sigmapowerphotodiode = data(11,:);
                
                mpd = mean(reshape(powerphotodiode,q.number_of_point,[]));
                inpower = mpd /max(mpd);

                [X Y] = meshgrid(-1:stp:1,-1:stp:1);
                Z=reshape(inpower,q.grid_size,[]);
                
                clear q mpd inpower stp;
                
            end
            obj.XMeshGrid=X;
            obj.YMeshGrid=Y;
            obj.ZMeshGrid=Z;
            clear X Y Z;
        end
    end
end