#include <octave/oct.h>

DEFUN_DLD( setNewFlagField,args, , "set new flag field")
{
   int           ly = args(0).int_value();
   int           lx = args(1).int_value();
   //int         init = args(2).int_value();
  
   Matrix      fill = args(4).matrix_value();
   Matrix    liquid = args(5).matrix_value();
   Matrix interface = args(6).matrix_value();   
   Matrix       gas = args(7).matrix_value();  
 
   Matrix    liquid_tmp(args(2).matrix_value());
   Matrix       gas_tmp(args(2).matrix_value());
   Matrix interface_tmp(args(2).matrix_value());
   
   for(int i = 1; i < ly+1; i++){
      for( int j=1;j< lx+1; j++){
      if ( fill(i,j) <= 0.0 )
         gas(i,j) = 1.0;
      else if ( (fill(i,j)>0.0) && (fill(i,j)<1.0) )
         interface(i,j) = 1.0;
      else if (fill(i,j) == 1.0 )
         liquid(i,j) = 1.0;
      }
   }
   liquid    = liquid_tmp;
   gas       = gas_tmp;
   interface = interface_tmp;

   return octave_value_list();

}
