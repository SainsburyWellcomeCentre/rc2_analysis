function val = constants(var)
%%list of parameters to use for analysis
%   Input:  var - name of the constant
%           Possible:   'spiking_class_threshold_ms'

consts.spiking_class_threshold_ms = 0.45;

val = consts.(var);
