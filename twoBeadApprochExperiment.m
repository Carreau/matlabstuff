classdef (ConstructOnLoad) twoBeadApprochExperiment < handle
    %properties, ie actually store variables
    properties
        rawdata
        still_trap
        moving_trap
        datevec_beggining
    end
    %properties that are recalculed when accessed (ie depend on rawdata)
    properties(Dependent)
        %angle between the x axis and the vector represented by the closest
        %end approch point os the moving trap and the start approch of the
        %moving trap xxxDeg et xxxRad ar convenient methor calling
        %approchAngle
        approachAngleDeg
        approachAngleRad
        approachAngle
        %wrapper to acces rawdata.f(:,:,:) with tangent and perpendicular
        %at the moving direction acessor (we should move moving trap itself
        %to a subclass with lazy accessor to parameters like position and
        %force
        moving_trap_force
        still_trap_force
        %
        time_m
    end
    methods
        %lazy accesor
        function f = get.moving_trap_force(self)
            f.x=squeeze(self.rawdata.f(2,1,:));
            f.y=squeeze(self.rawdata.f(2,2,:));
            f.tangent=f.x*cos(self.approachAngle)+f.y*sin(self.approachAngle);
            f.perp=f.x*sin(self.approachAngle)+f.y*cos(self.approachAngle);
        end
        %lazy accessor copy past from previous: not good
        function f = get.still_trap_force(self)
            f.x=squeeze(self.rawdata.f(1,1,:));
            f.y=squeeze(self.rawdata.f(1,2,:));
            f.tangent=f.x*cos(self.approachAngle)+f.y*sin(self.approachAngle);
            f.perp=f.x*sin(self.approachAngle)+f.y*cos(self.approachAngle);
        end
        %constructor
        function obj=twoBeadApprochExperiment(f)
            obj.rawdata=f;
            obj.datevec_beggining=datevec(0);
            g = f.trap_pos;
            obj.still_trap.x = squeeze(g(1,1,:));
            obj.still_trap.y = squeeze(g(1,2,:));
            obj.moving_trap.x = squeeze(g(2,1,:));
            obj.moving_trap.y = squeeze(g(2,2,:));
        end
        %lazy accesor (convenient)
        function theta=get.approachAngleDeg(self)
            theta=rad2deg(self.approachAngle);
        end
        
        %lazy accesor (convenient)
        function theta=get.approachAngleRad(self)
            theta=self.approachAngle;
        end
        
        %lazy accesor
        function theta=get.approachAngle(self)
            dx = self.moving_trap.x(1)-self.still_trap.x(1);
            dy = self.moving_trap.y(1)-self.still_trap.y(1);
            theta=angle(dx+1i*dy);
        end
        function setDatevec(self,year,month,day,h,m,s)
            for i=1:length(self)
                self(i).datevec_beggining=[year month day h m s];
            end
        end
        function t =get.time_m(self)
            t=datenum(self.rawdata.datevec-self.datevec_beggining)*24*60;
        end;
    end
end
