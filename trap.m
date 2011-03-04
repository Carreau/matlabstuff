classdef trap < handle
    properties
        kappa       %x et y
        pos         %x et y
        bead_pos    %x et y
        slopes      %x et y
        force       % x and y
        QPD_dx 
        QPD_dy
        QPD_sum
        aod_microm
    end
    properties(Dependent)
        absolute_bead_pos
    end
    
    methods
        function ret=get.absolute_bead_pos(self)
            ret.x = self.pos.x./self.aod_microm.x-1./self.slopes.x*self.QPD_dx/self.QPD_sum;
            ret.y = self.pos.y./self.aod_microm.y-1./self.slopes.y*self.QPD_dy/self.QPD_sum;
        end
        function obj=twoBeadApprochExperiment(all)
            obj.kappa.x      
            obj.kappa.y      
            obj.pos.x        
            obj.pos.y
            obj.bead_pos.x
            obj.bead_pos.y
            obj.slopes.x
            obj.slopes.y
            obj.force.x
            obj.force.y
            obj.QPD_dx
            obj.QPD_dy
            obj.QPD_sum
            obj.aod_microm.x
            obj.aod_microm.y
u
        end
        
    end 

end
