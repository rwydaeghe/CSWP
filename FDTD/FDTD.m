classdef cylinder
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        a = 1;
    end
    
    methods
        function obj = cyl(inputArg1,inputArg2)
            %UNTITLED Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = createCyl(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.a + inputArg;
        end
    end
end

