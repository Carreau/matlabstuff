classdef twoBeadApprochExperiment
    properties
        rawdata
        still_trap
        moving_trap
    end
    properties(Dependent)
        aproachAngleDeg
    end
    methods
        function obj=twoBeadApprochExperiment(f)
            obj.rawdata=f;
            g = f.trap_pos;
            obj.still_trap.x = squeeze(g(1,1,:));
            obj.still_trap.y = squeeze(g(1,2,:));
            obj.moving_trap.x = squeeze(g(2,1,:));
            obj.moving_trap.y = squeeze(g(2,2,:));
        end
        function theta=get.aproachAngleDeg(self)
            dx = self.moving_trap.x(1)-self.still_trap.x(1);
            dy = self.moving_trap.y(1)-self.still_trap.y(1);
            theta=rad2deg(angle(dx+i*dy));
        end
    end
end
