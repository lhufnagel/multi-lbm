classdef celltype
   enumeration
      Regular, 
      NorthMovSolid, 
      EastSolid, 
      SouthSolid, 
      WestSolid,
      EastPeriodic, 
      WestPeriodic,
      % Additional celltypes for multiphase LBM
      %Fluid1,
      %Fluid2,
      %Interface
   end
end
