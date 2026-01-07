classdef asic
    properties
        acc acc % accelerometer obj
        gyr gyr % gyroscope     obj
        tmp tmp % temperature sensor  obj
        str str % stress sensors obj
    end
    
    methods
        function obj = asic(nstr, pmode)
            % obj.acc = acc();
            obj.gyr = gyr(pmode);
            obj.tmp = tmp();
            obj.str = repmat(str(), nstr, 1);
        end
    end
end
