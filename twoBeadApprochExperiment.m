classdef (ConstructOnLoad) twoBeadApprochExperiment < handle
    %properties, ie actually store variables
    properties
        rawdata
        %traps position
        still_trap
        moving_trap
        datevec_beggining
        fTouchDetection
        fitvalue
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
        % time in minutes after datevec_beginning
        time_m
        %trap distance:
        trap_distance
        bead_distance
        %E
        E_final
        %touch distance
        touch_d
        param
    end
    %% calcul du module d'young 
    % %F=4/3*E/(1-(nu)^2).*d.^(3/2)*sqrt(r)+f0;
    % avec 
    % F force sur la bille
    % E module d'young
    % nu coefficient de poisson 
    % d ?
    % ?
    % f0
    methods
        function ret = description(self)
            ret=[];fprintf('filename\tyyyy-mm-dd hh:mm\tangl\tspeed\trest_time\n');
             for i=1:length(self)
                s=self(i);
                format = '';
                [pathstr, name, ext] = fileparts(s.rawdata.name);
                [~, sub , ~  ] = fileparts(pathstr);
                A=sprintf('%s/%s%s',sub(end-3:end),name(end-2:end),ext);
                format = strcat(format,'%s\t');

                A2 = s.rawdata.datevec(1:5);
                format = strcat(format,'%04i-%02i-%02i %02i:%02i\t');
    
                A3 = s.approachAngleDeg;
                format = strcat(format,'%4.0f\t');

                A4 = s.param.speed.value;
                format = strcat(format,'%3i\t');
     
                A5 = s.param.resting_time.value;
                format = strcat(format,'%3i\t');
    
                ret = [ret;sprintf(format,A,A2,A3,A4,A5)];
             end
        end
        function fit(self)
           for i=1:length(self)
                f = self(i).still_trap.force.r(1:self(i).moving_trap.event.appr.stop);
                d = self(i).bead_distance     (1:self(i).moving_trap.event.appr.stop);
                %plot(d,f);
                [d0,f0,E,err] = youngfit(d,f);
                self(i).fitvalue.d0  = d0;
                self(i).fitvalue.f0  = f0;
                self(i).fitvalue.E   = E*1e12;
                self(i).fitvalue.err = err*1E20;
           end
           
        end
        
        %lazy accessor
        function p = get.param(self)        
            p = self.rawdata.parameters;
        end
        function distance = get.touch_d(self)
            distance = self.trap_distance(self.rawdata.touch_point);
        end
        %lazy accessor
        function e = get.E_final(self)
           for i=1:length(self)
            e(i) = self(i).rawdata.E_final;
           end
        end
        %lazy accessor
        %function d = get.trap_distance(self,point)
        %   d = sqrt((self.moving_trap.x(point)+ self.still_trap.x(point)).^2 + (self.moving_trap.y(point)-self.still_trap.y(point)).^2);
        %end
        function d = get.trap_distance(self)
           %disp('this is the trap distance in aod units'); 
           %d = sqrt((self.moving_trap.x+ self.still_trap.x).^2 +
           %(self.moving_trap.y-self.still_trap.y).^2); 
           d= self.rawdata.d;
        end
        
        function d= get.bead_distance(self)
           mb = self.moving_trap.absolute_bead_pos;
           sb = self.still_trap.absolute_bead_pos;
           dx = sb.x-mb.x;
           dy = sb.y-mb.y;
           d = sqrt( dx .^2 + dy .^2); 
        end
        %lazy accesor
        function x = applyfun(self,fun)
            for i=1:length(self)
                x(i) = fun(self(i));
            end
            
        end;
        function f = get.moving_trap_force(self)
            f.x=squeeze(self.rawdata.f(2,1,:));
            f.y=squeeze(self.rawdata.f(2,2,:));
            f.tangent=(f.x * cos(self.approachAngle)) +(f.y *sin(self.approachAngle));
            f.perp   =(f.x *-sin(self.approachAngle)) +(f.y *cos(self.approachAngle));
        end
        %lazy accessor copy past from previous: not good
        function f = get.still_trap_force(self)
            f.x=squeeze(self.rawdata.f(1,1,:));
            f.y=squeeze(self.rawdata.f(1,2,:));
            f.tangent=f.x*cos(self.approachAngle)+f.y*sin(self.approachAngle);
            f.perp   =-f.x*sin(self.approachAngle)+f.y*cos(self.approachAngle);
        end
        %constructor
        function obj=twoBeadApprochExperiment(f)
            obj.rawdata=f;
            obj.datevec_beggining=datevec(0);
            g = f.trap_pos;
            obj.still_trap = f.still_trap;
            obj.moving_trap = f.moving_trap;
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
            dx = self.moving_trap.pos_um.x(1)-self.still_trap.pos_um.x(1);
            dy = self.moving_trap.pos_um.y(1)-self.still_trap.pos_um.y(1);
            theta=angle(dx+1i*dy);
            %theta=deg2rad(45);
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
