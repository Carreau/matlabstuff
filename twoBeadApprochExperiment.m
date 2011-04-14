classdef (ConstructOnLoad) twoBeadApprochExperiment < handle
    % class to handle two bead approch experiemtn 
    % to create use temptest(f.save_data)
    % the set date by calling self.setDatevec(YYYY,MM,DD,hh,mm,ss)
    % and self.setCommentaire('String');
    % do a self.fit
    %
    
    
    %properties, ie actually store variables
    properties
        rawdata
        %traps position
        still_trap
        moving_trap
        datevec_beggining
        fTouchDetection
        fitvalue
        commentaire
        arp
        cp
        % this will store the result for partial fit on a 
        % sliding two make an educated guess of the young modulus 
        % and the touch distance over the contact of the two bead
        partialfit
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
        %timo's young modulus
        % we assume the touchdistance is ok, and then we inverse the hertz
        % equation to get a young mudulus for each point
        Etimo
        
        
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
        %make the methode to guess E by timo's way
        function ret = get.Etimo(self)
            for i = 1:length(self)
               if(size(self(i).fitvalue)== [0 0])
                   disp('please fit the data first');
                   ret = Nan;
                   return
               end
               % let's find be the force at the touch point
               % then 
               % d0
               % f0
               % drift
               stop  = self(i).moving_trap.event.appr.stop;   
               f_ind = (self(i).fitvalue.drift*self(i).fitvalue.d0+self(i).fitvalue.f0);
               d00= self(i).fitvalue.d0*1.0;
               [~,j] = min(abs(self(i).bead_distance(1:stop)-d00));
               touchevent = j; 
               re=@(f,d) (f-f_ind)*3/4*(1-0.5^2) ./ (sqrt(4.5) .* (d00-d).^3/2);
               a=self(i).still_trap.force.r(touchevent:stop);
               b= self(i).bead_distance(touchevent:stop);
               ret = re(a,b);
               hold on;
               plot(b(1000:end)/d00,ret(1000:end),'r.');
               errorbar(4.5/d00,0,1e-12,'r');
            end
        end
        %set the capting proteine concentration
        function setCp(self,value)
            for i = 1:length(self)
                self(i).cp = value;
            end
        end
        function setArp(self,value)
            for i = 1:length(self)
                self(i).arp = value;
            end
        end
        function setCommentaire(self,str)
            for i = 1:length(self)
               self(i).commentaire=str;
            end
        end
        function ret = description(self)
            ret=[];fprintf('filename\tyyyy-mm-dd hh:mm\tangl\tspeed\tr_time\ttime_m\tkx/akx\tky/aky\n');
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

                A6 = floor(s.time_m*10)/10;
                format = strcat(format,'%04.1f\t');
                
                A8 = s.still_trap.kappa.x/s.still_trap.auto_kappa.x;
                A7 = s.still_trap.kappa.y/s.still_trap.auto_kappa.y;
                format = strcat(format,'%4.1f\t%4.1f\t');
                

                ret = [ret;sprintf(format,A,A2,A3,A4,A5,A6,A7,A8)];
             end
        end
        function doPartialFit(self)
            for i=1:length(self)
                clf
                s= self(i);
                s.partialfit=[];
                f = self(i).still_trap.force.r(1:self(i).moving_trap.event.appr.stop);
                d = self(i).bead_distance     (1:self(i).moving_trap.event.appr.stop);
                n = length(f);
                fprintf(1,'we hav %d value to fit\n',length(f));
                step = 2000;
                sb = 200;
                a = floor(step/sb);
                n = floor(n/step);
                dinit=30;
                Einit=1e-15;
                f0init=0;
                Starting = [dinit,Einit];
                for j=1:(n-1)*a
                    intv = (1+(j-1)*sb:1+(j-1)*sb+step);
                    [d0,E,err,flag] = yf(d(intv),f(intv),Starting);
                    if(flag)
                        
                        Starting = [d0,E];
                        s.partialfit.d0(j)    = d0;
                        s.partialfit.E(j)     = E*1e12;
                        s.partialfit.err(j)   = err*1E20;
                        s.partialfit.d(j) = mean(d(intv));
                    else 
                        disp('max reach /////////');
                    end
                end
            end
        end
        function fit(self)
           for i=1:length(self)
                f = self(i).still_trap.force.r(1:self(i).moving_trap.event.appr.stop);
                d = self(i).bead_distance     (1:self(i).moving_trap.event.appr.stop);
                %l=length(f);
                %n=50;
                %f0 = f(1:floor(l/n):end);
                %d0 = d(1:floor(l/n):end);
                [d0,f0,E,err,drift] = youngfit(d,f);
                self(i).fitvalue.d0  = d0;
                self(i).fitvalue.f0  = f0;
                self(i).fitvalue.E   = E*1e12;
                self(i).fitvalue.err = err*1E20;
                self(i).fitvalue.drift = drift;
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
