classdef trap < handle 
    properties
        kappa       %x et y
        pos_aod         %x et y
        slopes      %x et y
        % lin_binnig add a variable to return the results linbinned by a
        % certain amount
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
        %% methods posfixed with ld send lin_bined values 
        event_lb
        bead_pos_in_trap_lb
        force_lb
        pos_um_lb
        pos_lb
        absolute_bead_pos_lb
        
        %%
        event
        bead_pos_in_trap
        force
        pos_um
        pos
        absolute_bead_pos
    end
    
    methods
        function ret=get.bead_pos_in_trap_lb(self)
            f= self.bead_pos_in_trap;
            ret.x = linbin(f.x,self.linbinvalue);
            ret.y = linbin(f.y,self.linbinvalue);
        end
        function ret=get.force_lb(self)
            f= self.force;
            ret.x = linbin(f.x,self.linbinvalue);
            ret.y = linbin(f.y,self.linbinvalue);
        end
        function ret=get.pos_um_lb(self)
            f= self.pos_um;
            ret.x = linbin(f.x,self.linbinvalue);
            ret.y = linbin(f.y,self.linbinvalue);
        end
        function ret=get.pos_lb(self)
            f= self.pos_aod;
            ret.x = linbin(f.x,self.linbinvalue);
            ret.y = linbin(f.y,self.linbinvalue);
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
        end
        function ret = get.bead_pos_in_trap(self)
            ret.x = -1./self.slopes.x*self.QPD_dx./self.QPD_sum;
            ret.y = -1./self.slopes.y*self.QPD_dy./self.QPD_sum;
        end
        function ret=get.absolute_bead_pos(self)
            ret.x = self.pos_aod.x./self.aod_microm.x + self.bead_pos_in_trap.x;
            ret.y = self.pos_aod.y./self.aod_microm.y + self.bead_pos_in_trap.y;
        end
        function ret=get.force(self)
            %1E-6 car kappa doit être en pN/m ou N/mu m
            ret.x = self.kappa.x/self.slopes.x*self.QPD_dx ./ self.QPD_sum * 1E-6;
            ret.y = self.kappa.y/self.slopes.y*self.QPD_dy ./ self.QPD_sum * 1E-6;
        end
    end

end

