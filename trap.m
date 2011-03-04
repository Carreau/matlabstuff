classdef trap < handle 
    properties
        kappa       %x et y
        pos_aod         %x et y
        bead_pos    %x et y
        slopes      %x et y
        QPD_dx 
        QPD_dy
        QPD_sum
        aod_microm
    end
    properties(Dependent)
        force
        pos_um
        pos
        absolute_bead_pos
    end
    
    methods
        function ret = get.pos_um(self)
            ret.x = self.pos_aod.x ./ self.aod_microm.x; 
            ret.y = self.pos_aod.y ./ self.aod_microm.y;
        end
        function ret=get.absolute_bead_pos(self)
            ret.x = self.pos_aod.x./self.aod_microm.x-1./self.slopes.x*self.QPD_dx./self.QPD_sum;
            ret.y = self.pos_aod.y./self.aod_microm.y-1./self.slopes.y*self.QPD_dy./self.QPD_sum;
        end
        function ret=get.force(self)
            %1E-6 car kappa doit être en pN/m ou N/mu m
            ret.x = self.kappa.x/self.slopes.x*self.QPD_dx ./ self.QPD_sum * 1E-6;
            ret.y = self.kappa.y/self.slopes.y*self.QPD_dy ./ self.QPD_sum * 1E-6;
        end
    end

end

