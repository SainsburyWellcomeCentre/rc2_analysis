classdef LFPPowerPlot < RC2Axis
    
    properties
        
        ap
        
        
    end
    
    methods
        
        function obj = LFPPowerPlot(h_ax)
            
            VariableDefault('h_ax', []);
            
            obj = obj@RC2Axis(h_ax);
        end
        
        
        
    end
end
