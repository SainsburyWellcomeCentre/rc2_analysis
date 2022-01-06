classdef Probe < handle
% Probe Class for handling geometrical information about probe
%
%   Probe Properties:
%       n_shanks                    - number of shanks on the probe
%       electrodes_per_shank        - number of electrodes per shank
%       electrode_size_x            - size of electrode horizontally (um)
%       electrode_size_y            - size of electrode vertically (um)
%       reference_electrode_ids     - ID of the reference electrode
%       tip_length_um               - length of probe tip below lowest electrode (um)
%       electrode_y_separation_um   - separation of electrodes vertically (um)
%       electrode_x_separation_um   - separation of electrodes horizontally (um)
%       shank_x_separation_um       - separation of shanks horizontally (um)
%       n_electrode_columns_per_shank - number of columns of electrodes horizontally
%       bad_electrode_ids           - list of bad electrode IDs 
%
%   Probe Methods:
%       electrode_from_tip_um       - compute distance of electrode from probe tip
%       electrode_from_left_um      - compute distance of electrode from leftmost electrode
%
%   Here we work with electrode IDs rather than channel IDs. To convert you
%   can use the SpikeGLXMetaData class.
%
%   See also: SpikeGLXMetaData, SpikeGLXRecording

    properties
        
        type
    end
    
    properties (Dependent = true)
        
        n_shanks
        electrodes_per_shank
        electrode_size_x
        electrode_size_y
        reference_electrode_ids
        tip_length_um
        electrode_y_separation_um
        electrode_x_separation_um
        shank_x_separation_um
        n_electrode_columns_per_shank
        bad_electrode_ids
    end
    
    
    
    methods
        
        function obj = Probe(type)
        % Probe
        %
        %   Probe(TYPE) class for querying geometrical properties of different probe types
        %   TYPE is either '3A' or '24'
        
            obj.type = type;
        end
        
        
        
        function val = get.n_shanks(obj)
        %%number of shanks for probe
        
            if strcmp(obj.type, '3A')
                val = 1;
            elseif strcmp(obj.type, '24')
                val = 4;
            end
        end
        
        
        
        function val = get.n_electrode_columns_per_shank(obj)
        %%number of electrode columns per shank
        
            if strcmp(obj.type, '3A')
                val = 4;
            elseif strcmp(obj.type, '24')
                val = 2;
            end
        end
        
        
        
        function val = get.electrode_size_x(obj)
        %size of electrode horizontally (in um)
            val = 12;
        end
        
        
        
        function val = get.electrode_size_y(obj)
        %size of electrode vertically (in um)
            val = 12;
        end
        
        
        
        function val = get.electrodes_per_shank(obj)
        %number of electrodes per shank
        
            if strcmp(obj.type, '3A')
                val = 960;
            elseif strcmp(obj.type, '24')
                val = 1280;
            end
        end
        
        
        
        function val = get.reference_electrode_ids(obj)
        %ID of the electrodes which are used as reference channels
        
            if strcmp(obj.type, '3A')
                % from phase 3A manual, subtract 1 as in that manual they
                % label from 1 but we will index from 0
                val = [ 37,  76, 113, 152, 189, 228, 265, 304, 341, 380, ...
                       421, 460, 497, 536, 573, 612, 649, 688, 725, 764, ...
                       805, 844, 881, 920, 957] - 1;  
            elseif strcmp(obj.type, '24')
                val = [127, 511, 895, 1279];
            end
        end
        
        
        
        function val = get.tip_length_um(obj)
        %length of the probe tip below lowest electrode
            if strcmp(obj.type, '3A')
                val = 200;
            elseif strcmp(obj.type, '24')
                val = 175;  % double check this value
            end
        end
        
        
        
        function val = get.electrode_y_separation_um(obj)
        %separation of electrodes vertically (um)
            if strcmp(obj.type, '3A')
                val = 20;
            elseif strcmp(obj.type, '24')
                val = 15;
            end
        end
        
        
        
        function val = get.electrode_x_separation_um(obj)
        %separation of electrodes on a shank horizontally (um)    
            if strcmp(obj.type, '3A')
                val = 32;
            elseif strcmp(obj.type, '24')
                val = 32;
            end
        end
        
        
        
        function val = get.shank_x_separation_um(obj)
        %separation of shanks in um
            if strcmp(obj.type, '3A')
                val = nan;
            elseif strcmp(obj.type, '24')
                val = 250;
            end
        end
        
        
        
        function val = electrode_from_tip_um(obj, electrode_id)
        %%electrode_from_tip_um Compute distance of electrode from probe tip
        %
        %    DISTANCE = electrode_from_tip_um(ELECTRODE_ID)
        %    computes the distance in um of an electrode or vector of electrode IDs
        %    (ELECTRODE_ID) from the probe tip. Electrode IDs are zero indexed.
            
            row = floor(electrode_id(:)/2); % 0 indexed
            val = obj.tip_length_um + row * obj.electrode_y_separation_um;
        end
        
        
        
        function val = electrode_from_left_um(obj, electrode_id, shank_id)
        %%electrode_from_left_um Compute distance of electrode from leftmost electrode
        %
        %
        %    DISTANCE = electrode_from_left_um(ELECTRODE_ID, SHANK_ID)
        %    computes the distance in um of an electrode or vector of electrode IDs
        %    (ELECTRODE_ID) on shanks with ID SHANK_ID from 
        %    the left-most electrode on *the probe* (i.e. shank ID 0). 
        %    Electrode IDs are zero indexed.
    
            row = floor(electrode_id(:)/2);
            col = mod(electrode_id(:), 2);
            val = col * obj.electrode_x_separation_um();
            
            if strcmp(obj.type, '3A')
                % sites are staggered on 3A so if even row shift by 16um
                idx = mod(row, 2) == 0;
                val(idx) = val(idx) + 16;
            elseif strcmp(obj.type, '24')
                % offset by the shank
                val = shank_id * obj.shank_x_separation_um + val;
            end
        end
    end
end
