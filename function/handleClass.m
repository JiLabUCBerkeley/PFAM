classdef handleClass < handle
    properties
        Shifts
        scatterPoints
        PosRef
        SH_im
        SH_im_bw
        Param
    end 
    methods
        function h = handleClass( Shifts, PosRef, scatterPoints, SH_im, Param)
          h.Shifts = Shifts;
          h.scatterPoints = scatterPoints;
          h.PosRef = PosRef;
          h.SH_im = SH_im;
          h.Param = Param;
        end
    end
end