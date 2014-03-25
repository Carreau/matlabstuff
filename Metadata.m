classdef (ConstructOnLoad) Metadata < handle
    
    %properties, ie actually store variables
    properties
        nom
        hauteur
        rves
        rpip
    end
    %properties that are recalculed when accessed (ie depend on rawdata)
    properties(Dependent)
        tension

    end
        methods
        
        function t = get.tension(self)
            t = self.rves/self.rpip;
        end


        %constructor
        function self=Metadata(data)
            self.nom = data.nom     
            self.hauteur = data.hauteur
            self.rves = data.rves
            self.rpip = data.rpip
            % si dqtq.UUID existe alors 
            % self.UUID = data.uui
            % sinon 
            % self.uuid = java.util.UUID.randomUUID.toString
        end

        function s = serialize(self)
            s.nom = self.nom    ; 
            s.hauteur = self.hauteur;
            s.rves = self.rves  ;
            s.rpip = self.rpip  ;
            s.uuid = self.uuid
        end

        
    end
end


