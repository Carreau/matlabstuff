classdef trap < handle 
    properties
        kappa       %x et y
        corrected_kappa = false % bool, which return corrected kappa if asked
        pos_aod         %x et y
        cpos_aod
        slopes      %x et y
        % lin_binnig add a variable to return the results linbinned by a
        % certain amount
        trapInterpolateur
        QPD_dx 
        QPD_dy
        QPD_sum
        aod_microm
        linbinvalue = 0; 
    end
    properties(Dependent)
        %sould add a time function to be able to over plot files with
        %différent rates.
        % event.appr.start
        % event.appr.stop
        % event.retr.start
        %% methods posfixed with ld send lin_bined values, cf without nb
        %% for description 
        
        event_lb
        bead_pos_in_trap_lb
        force_lb
        pos_um_lb
        pos_lb
        absolute_bead_pos_lb
        
        %% 
        % event : store approch start, stop, and retr start
        %     
        %kappa
        event
        bead_pos_in_trap    % in um 
        force               % x and y for now (should implement theta phi)
        pos_um              % same
        pos                 % not implemented should decide in which unity
        absolute_bead_pos   % in um (auto depend on kappa.x(.y)...)
        auto_kappa          % for now, this will be a trap stiffness rescaled by the QPD power
    end
    
    methods
        %let's do a kappa rescaled by QPD sum
        function setkappa(self,value)
            self.rkappa = value;
        end
%         function ret=get.pos_aod(self)
%             if isempty(self.cpos_aod)
%                 dig2noru=@(d)(((d/256)*500/2^23)-75)/15;
%                 nor2digu= @(n) (((n*15+75)*2^23)/500+2^23)*256;
%                 correctn2ds = @(x) dig2noru(nor2digu(x));
%                 correctn2da = @(u)arrayfun(correctn2ds,u);
%                 self.cpos_aod.x= correctn2da(self.pos_aod.x);
%                 self.cpos_aod.y= correctn2da(self.pos_aod.y);
%             end
%             ret = self.cpos_aod;
%         end
        function set.trapInterpolateur(self,q)
            self.trapInterpolateur = trapPowerInterpolation(q);
        end
        function ret=ckappa(self)
            if self.corrected_kappa
                pf=self.trapInterpolateur.interpolatedPower(self.pos_aod.x(1),self.pos_aod.y(1));
                ret.x=self.kappa.x*pf;
                ret.y=self.kappa.y*pf;
            else
                ret = self.kappa;
            end
        end
        function ret=cslopes(self)
            if self.corrected_kappa
                pf=self.trapInterpolateur.interpolatedPower(self.pos_aod.x(1),0*self.pos_aod.y(1));
                ret.x=self.slopes.x*pf;
                ret.y=self.slopes.y*pf;
            else
                ret = self.slopes;
            end
        end
        function ret=get.kappa(self)
            if self.corrected_kappa == true 
                ret=self.ckappa;
            else
                ret=self.kappa;
            end
        end
        function ret = get.auto_kappa(self)
            % i'll use 30_mars_50cp_25apr.mat as reference
            % is 
            % for a qpd mean of 3.1136
            % kappa.x will be  3.4466e-05
            % kappa.y wil be 4.1101e-05
            ret.x = mean(self.QPD_sum(1:1000))/3.1136*3.4466e-05;
            ret.y = mean(self.QPD_sum(1:1000))/3.1136*4.1101e-05;
        end
        %fonction to linbin 4parameters array with x,y,r et theta
        function ret = lbxyrt(self,data) %LinBin X Y R Theta
            ret.x = linbin(data.x,self.linbinvalue);
            ret.y = linbin(data.y,self.linbinvalue);
            ret.r = linbin(data.r,self.linbinvalue);
            ret.theta = linbin(data.theta,self.linbinvalue);
        end
        function ret=get.bead_pos_in_trap_lb(self)
            ret = lbxyrt(self,self.bead_pos_in_trap);
        end
        function ret=get.force_lb(self)
            ret = lbxyrt(self,self.force);
        end
        function ret=get.pos_um_lb(self)
            ret = lbxyrt(self,self.pos_um);
        end
        function ret=get.pos_lb(self)
            ret = lbxyrt(self,self.pos_aod);
        end
        
     function ret = get.event_lb(self)
            tr=self.pos_um_lb;
            tr_m=sqrt(   squeeze((tr.x)').^2 + squeeze((tr.y)').^2 );
            tr_d_0=find(diff(tr_m)==0);
            tr_d_no0=find(diff(tr_m)~=0);
            ret.appr.start=tr_d_no0(1);
            ret.retr.start=tr_d_0(end);
            tr_d_0_cut=tr_d_0;
            %tr_d_0_cut(find(tr_d_0_cut<=tr_d_no0(1)))=[];
            tr_d_0_cut(tr_d_0_cut<=tr_d_no0(1))=[];
            ret.appr.stop=tr_d_0_cut(1);
        end


        function ret = get.event(self)
            tr=self.pos_um;
            tr_m=sqrt(   (tr.x).^2 + (tr.y).^2 );
            tr_d_0=find(diff(tr_m)==0);
            tr_d_no0=find(diff(tr_m)~=0);
            ret.appr.start=tr_d_no0(1);
            ret.retr_start=tr_d_0(end);
            tr_d_0_cut=tr_d_0;
            %tr_d_0_cut(find(tr_d_0_cut<=tr_d_no0(1)))=[];
            tr_d_0_cut(tr_d_0_cut<=tr_d_no0(1))=[];
            ret.appr.stop=tr_d_0_cut(1);
        end
        
        
        
        function ret = get.pos_um(self)
            ret.x = self.pos_aod.x ./ self.aod_microm.x; 
            ret.y = self.pos_aod.y ./ self.aod_microm.y;
            [ret.theta,ret.r] = cart2pol(ret.x,ret.y);
        end
        function ret = get.bead_pos_in_trap(self)
            ret.x = -1./self.slopes.x*self.QPD_dx./self.QPD_sum;
            ret.y = -1./self.slopes.y*self.QPD_dy./self.QPD_sum;
            [ret.theta,ret.r] = cart2pol(ret.x,ret.y);
        end
        function ret=get.absolute_bead_pos(self)
            ret.x = self.pos_aod.x./self.aod_microm.x + self.bead_pos_in_trap.x;
            ret.y = self.pos_aod.y./self.aod_microm.y + self.bead_pos_in_trap.y;
            [ret.theta,ret.r] = cart2pol(ret.x,ret.y);
        end
        function ret=get.force(self)
            %1E-6 car kappa doit être en pN/m ou N/mu m
            ret.x = self.kappa.x/self.slopes.x*self.QPD_dx ./ self.QPD_sum * 1E-6;
            ret.y = self.kappa.y/self.slopes.y*self.QPD_dy ./ self.QPD_sum * 1E-6;
            [ret.theta,ret.r] = cart2pol(ret.x,ret.y);
        end
    end

end

