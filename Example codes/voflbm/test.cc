#include <octave/oct.h>


DEFUN_DLD (test, args, nargout,
           "test")
{  
  //int number1 = args(0).int_value();
  //int number2 = args(1).int_value();
  Matrix matrix_1 = args(0).matrix_value();
  int lx        = args(1).int_value();
  int ly   	= args(2).int_value();
  int num       = args(3).int_value();
  
  int result = matrix_1( 0 % ly, 0 % lx);


  return octave_value(result);
}
