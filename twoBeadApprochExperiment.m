classdef (ConstructOnLoad) twoBeadApprochExperiment < handle
    % Class to handle two bead approch experiement
    %
    % CREATING OBJECT
    %      to create use
    %
    %      >> exp = temptest('experimentfile.mat')
    %
    %      if a file named 'experimentfile_meta.mat' il will also be loaded to
    %      add information about actin melange time, concentration of caping
    %      protein, arp...usw
    %
    %      if the _meta file does not exist il will be created after you enter
    %      the information
    %
    % STATIC DATA STRUCTURE
    %      rawdata
    %      still_trap  -> structure to acces info about the still trap
    %      moving_trap -> structure to acces info about the movong trap
    %            cf trap.m
    %      datevec_beggining -> datevec of actin mix, use to caculate time_m
    %      fTouchDetection -> not Used anymore, keep for compatibility
    %      fitvalue        -> store the fitvalue after .fit is called
    %      commentaire     -> misc caracter string
    %      arp             ->concentration in arp
    %      cp              -> concentrationin cp protein
    %      partialfit      -> fit Values for a slindingFit after doPartialFit is
    %                         called
    % DYNAMIC DATA STRUCTURE
    %      % those data type ara calculated each time you call them and are
    %      % dependant on the values of the previous one
    %
    %      approachAngleDeg
    %      approachAngleRad
    %      approachAngle
    %      moving_trap_force -> almost same as still_trap.force but
    %                           convenient acces to .tangent and .perp
    %                           composantes
	%      still_trap_force  -> same
    %                   % note : x and y might be inverted and this one is line
    %                   % bined
    %
    %      time_m            -> time in minute after value in datevec
	%      trap_distance     -> distance between the two traps
    %      bead_distance     -> distance between the two beads
	%      E_final           -> not used
	%      touch_d           -> not used
    %      param             -> convenient acces to rawdata.parameters
    %      Etimo             -> timo's way of findig E by averaging. need .fit
    %                           before
    %      % function %
    %      doPartialFit      -> run an slinding fit
    %      fit               -> fit the curve
    %      description       -> plot a short array description of each exp.
    %      applyfun          -> apply a anonypous function on each object
    %                           usefull to get acces to subparameters
    %                           eg
    %      >> exp.applyfun(@(x) x.still_trap.force.x(1))
    %            will give you an array of the fist value of the force on the
    %            still_trap for each experiement
    %      >> figure(2)
    %      >> hold on
    %      >> g.applyfun(@(x) plot(x.bead_distance,x.moving_trap.force.x,'+'));
    %            will plot the force on the moving trap as a fonction of the
    %            bead distance for all experiement.

    %properties, ie actually store variables
    properties
        rawdata
        index
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
        metafile
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
        
        function showfitTau(self)
           for i = 1:length(self)
                e = self(i).moving_trap.event.appr.stop;
                ee = self(i).moving_trap.event.retr_start;
                d =[0:ee-e]*1/self(i).rawdata.parameters.Sampling_rate.value;
                f = [self(i).still_trap.force.r];
                a = self(i).fitvalue.estimates.a;
                b = self(i).fitvalue.estimates.b;
                tau = self(i).fitvalue.estimates.tau;
                fh = @(x)a*exp(-x./tau)+b;
                figure(2);
                clf;
                hold on
                plot(d,f(e:ee),'+');
                hold on
                plot(d,fh(d),'r+');
           end
        end
        
        function showfitYoung(self)
           for i = 1:length(self)
                e = self(i).moving_trap.event.appr.stop;
                d = [self(i).bead_distance];
                f = [self(i).still_trap.force.r];
                fh = @(x) youngHertz(x,self(i).fitvalue.d0,self(i).fitvalue.f0,self(i).fitvalue.E*1e-12,self(i).fitvalue.drift);
                figure(1);
                hold off;
                plot(d(1:e),f(1:e),'+');
                hold on;
                plot(d(1:e),fh(d(1:e)),'r-');
           end
        end
        %% make the methode to guess E by timo's way
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
               b= b(1000:end);
               ret=ret(1000:end);
               w= (ret/d00>0);

               b=b(w);
               ret=ret(w);
               step=30;


               %plot(b(1:step:end)/d00,ret(1:step:end),'g.');
               %errorbar(4.5/d00,0,1e-12,'r');
               ret = [a b];
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
                k=0;
                for j=1:(n-1)*a
                    intv = (1+(j-1)*sb:1+(j-1)*sb+step);
                    [d0,E,err,flag] = yf(d(intv),f(intv),Starting);
                    if(flag)
                        k=k+1;
                        Starting = [d0,E];
                        s.partialfit.d0(k)    = d0;
                        s.partialfit.E(k)     = E*1e12;
                        s.partialfit.err(k)   = err*1E20;
                        s.partialfit.d(k) = mean(d(intv));
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
                l=length(f);
                n=50;
                f0 = f(1:floor(l/n):end);
                d0 = d(1:floor(l/n):end);
                [d0,f0,E,err,drift] = youngfit(d,f);
                self(i).fitvalue.d0  = d0;
                self(i).fitvalue.f0  = f0;
                self(i).fitvalue.E   = E*1e12;
                self(i).fitvalue.err = err*1E20;
                self(i).fitvalue.drift = drift;
                %% let's calculate the relaxation time
                t0 = self(i).moving_trap.event.appr.stop;
                t1 = self(i).moving_trap.event.retr_start;
                d_t = self(i).still_trap.force.r(t0:t1);
                %%
                d_t = [self(i).still_trap_force.tangent]';
                d_t=self(i).rawdata.d(self(i).rawdata.appr_stop:self(i).rawdata.retr_start);
                d_xxx=self(i).rawdata.d();
                
                figure(4)
                plot(d_xxx,'+');
                title('distance, time');
                
                %   1           -   
                %    t_t=[0:length(d_t)-1]*1/self(i).rawdata.parameters.Effective_Sampling_Rate.value;
                t_t=[0:length(d_t)-1]*1/self(i).rawdata.parameters.Sampling_rate.value;%average per trap ?    
                start_point = [(d_t(1)-d_t(end)),d_t(end),t_t(end)/40];
                %options=optimset('iter');
                
                [a,err] = fminsearch(@exp_dec, start_point,[],t_t,d_t);
                figure(3);
                hold off;
                plot(t_t,d_t,'+');
                title('décroissance de la force lors du repos des pièges.')
                xlabel('temps');
                ylabel('force');
                hold on;
                aa= @(x) a(1)*exp(-x./a(3))+a(2);
                bb =@(x) arrayfun(aa,x);
                plot(t_t,bb(t_t),'r','LineWidth',2);
                self(i).fitvalue.estimates.a = a(1);
                self(i).fitvalue.estimates.b = a(2);
                self(i).fitvalue.estimates.tau = a(3);
                self(i).fitvalue.estimates.err = err;
                fprintf('fit %d/%d\n',i,length(self));
                f= load(self(i).metafile);
                f.submeta(self(i).index).fitvalue= self(i).fitvalue;
                %disp 'save in' self(i).metafile;
                save(self(i).metafile,'-struct','f');
                clear f;
                %self(i).showfitTau;
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
            x=zeros(size(self));
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
            %f.x=squeeze(self.rawdata.f(1,1,:));
            %f.y=squeeze(self.rawdata.f(1,2,:));
            f.x=self.still_trap.force.x;
            f.y=self.still_trap.force.y;
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

function [sse, FittedCurve] =exp_dec(params,t,f)
        %global eta2 fc
        a = params(1);
        b = params(2);
        tau= params(3);
        FittedCurve=a*exp(-t./tau)+b;
        ErrorVector = (FittedCurve - f)./(f);
        %ErrorVector = (FittedCurve - f).*(FittedCurve-b);
        %ErrorVector = (FittedCurve - f)./(FittedCurve);
        ErrorVector(isnan(ErrorVector))=[];
        sse = abs(sum(abs(ErrorVector) .^ 2))*1e10;
end
