classdef SparseNoiseSession < RC2Session
% SparseNoiseSession Class for handling details of a
% sparse_noise session.
%
%   SparseNoiseSession Properties:
%
%   SparseNoiseSession Methods:
%
%   Currently this is a placeholder for handling sessions of pure visual
%   stimuli.
%
%   See also: RC2Session   

    properties
    end
    
    
    methods
        
        function obj = SparseNoiseSession(session)
        % SparseNoiseSession
        %
        %   SparseNoiseSession(SESSION) takes a session structure,
        %   SESSION, saved in the formatted data file and creates an object
        %   to handle the data of a 'sparse_noise' session.
        
            obj = obj@RC2Session(session);
        end
    end
end
