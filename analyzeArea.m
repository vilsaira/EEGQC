%{
EEGQC v5: automated tool for EEG quality analysis.

Copyright (C) 2013, 2015, 2016, Viljami Sairanen (viljami.sairanen@gmail.com)

This software is published under the BSD 2-Clause License
and is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%}

classdef analyzeArea
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        startTime = 0;
        duration = 20;
        window = 2;
        active = true;
    end
    
    methods
        
        %function obj = analyzeArea(startTime, duration, window, active)
        function obj = analyzeArea(startTime)
            if (nargin == 1)
            %switch nargin
            %    case 1
                    obj.startTime = startTime;
            %    case 2
            %        obj.startTime = startTime;
            %        obj.duration = duration;
            %    case 3
            %        obj.startTime = startTime;
            %        obj.duration = duration;
            %        obj.window = window;
            %    case 4
            %        obj.startTime = startTime;
            %        obj.duration = duration;
            %        obj.window = window;                
            %end
            end
        end
        
        function obj = changeStartTime(obj, change)
            obj.startTime = obj.startTime + change;
        end
        
        function obj = changeDuration(obj, change)
            obj.duration = obj.duration + change;
        end
        
        function obj = changeWindow(obj, change)
            obj.window = obj.window + change;
        end

    end
    
end

