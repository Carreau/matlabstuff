classdef (ConstructOnLoad) erasableBuffer < handle
% this class should play a the role of a rewritable buffer which should
% look like a modifing string in the console
% it should retain the number of carracter written and send enough \b, on
% next write

    properties
        currentstring='';
        laststring
    end
    methods
        function clear(self)
           for i=1:length(self.laststring);
                fprintf('\b');
           end
           %self.laststring='';
        end
        function print(self,string)
            self.clear();
            fprintf(string);
            self.laststring=string;
        end 
        
        function counter(self,i,m)
              if i==1
                  fprintf('\n');
              end
              self.print(sprintf('working ... %d/%d',i,m));
              if i==m
                  fprintf('\n');
              end
        end
    end
    
end