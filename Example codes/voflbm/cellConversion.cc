#include <octave/oct.h>

// TODO: check boundary because of periodicity!!!! -> change cell conversions!!!

DEFUN_DLD( cellConversion,args,, "do necessary cell conversions")
{
   int    ly     	= args(0).int_value();
   int    lx     	= args(1).int_value();
   Matrix mass      	= args(2).matrix_value();
   Matrix fill      	= args(3).matrix_value();
   Matrix liquid    	= args(4).matrix_value();
   Matrix interface	= args(5).matrix_value();   
   Matrix gas       	= args(6).matrix_value();  
   Matrix conv_I_L  	= args(7).matrix_value();  
   Matrix conv_I_G  	= args(8).matrix_value();
   Matrix fValues       = args(9).matrix_value();
   int    b             = args(10).int_value();
   RowVector weight     = args(11).row_vector_value();
   RowVector cx         = args(12).row_vector_value();
   RowVector cy         = args(13).row_vector_value();  
   
   octave_value_list retval;

   const double tol = 0.01; 

   //---------------------------------------------------------------------//
   // 4 possible cell conversions   
   for(int i = 1; i < ly+1; i++){
      for( int j=1;j< lx+1; j++){
 
      int i_left  = ((i-2) % ly) + 1; // i-1  
      int i_right = (    i % ly) + 1; // i+1

      int j_left  = ((j-2) % lx) + 1; // j-1
      int j_right = ( j    % lx) + 1; // j+1

      // I -> L
      // TODO: conv_I_L = zeros(ly,lx) -> otherwise old values!!!!
      if ( (fill(i,j)>= 1.0 + tol) && (interface(i,j)== 1) )  	// convert cell to liquid cell 
      {
         interface(i,j) = 0;
         liquid(i,j)    = 1;       	   			// counter for gas cells in neighborhood -> how we can distribute mass
         int counter    = 0;
	 
              if( gas(i_right,j)==1)    	 			// search for gas cells! -> convert them to interface cells
              {
                    gas(i_right,j) 		= 0;   	 		// generate new interface cell!
                    interface(i_right,j) 	= 1;
		    conv_I_L(i_right,j) 	= 1;	 		// inter-Flag: remember new interface flags to know which interface flags get mass!
                    counter 			= counter +1;
              }
              else if( gas(i, j_left) == 1)
              {
                    gas(i,j_left) 		= 0;
                    interface(i,j_left) 	= 1;
	            conv_I_L(i,j_left)	 	= 1;
                    counter 			= counter+1;
              }
              else if( gas(i_left,j) ==1 )
              {
		    gas(i_left,j)	    	= 0;
		    interface(i_left,j)		= 1;
		    conv_I_L(i_left,j)		= 1;
                    counter         	        = counter +1;
	      }
              else if( gas(i,j_right) == 1)
              {
    		    gas(i,j_right)	     	= 0;
		    interface(i,j_right) 	= 1;
                    conv_I_L(i,j_right)  	= 1;
                    counter          		= counter + 1;
	      }
	      else if( gas(i_right,j_left) == 1 )
              {
		    gas(i_right,j_left)      	= 0;
	            interface(i_right,j_left) 	= 1;
		    conv_I_L(i_right,j_left)  	= 1;
		    counter            		= counter + 1;
	      } 
              else if( gas(i_left,j_left) == 1)
	      {
		    gas(i_left,j_left)       	= 0;
		    interface(i_left,j_left) 	= 1;
		    conv_I_L(i_left,j_left)  	= 1;
                    counter            		= counter + 1;
	      }
              else if( gas(i_left,j_right) == 1)
              {
		    gas(i_left,j_right)       	= 0;
		    interface(i_left,j_right) 	= 1;
		    conv_I_L(i_left,j_right)  	= 1;
		    counter            		= counter + 1;
	      }
              else if( gas(i_right,j_right) == 1)
              {
		    gas(i_right,j_right)       	= 0;
		    interface(i_right,j_right) 	= 1;
		    conv_I_L(i_right,j_right)  	= 1;
		    counter            		= counter + 1;
	      } 
          if(counter != 0)
          {
             double masspNewInter = (fill(i,j)-1.0)/counter;
             // distribute mass to new interface cells
             // do we have to run through all cells again????
             if(      conv_I_L(i_right,j      ) == 1)  mass(i_right,j	   ) = masspNewInter;
             else if( conv_I_L(i      ,j_left ) == 1)  mass(i      ,j_left ) = masspNewInter;
             else if( conv_I_L(i_left ,j      ) == 1)  mass(i_left ,j	   ) = masspNewInter;
             else if( conv_I_L(i      ,j_right) == 1)  mass(i      ,j_right) = masspNewInter;
             else if( conv_I_L(i_right,j_left ) == 1)  mass(i_right,j_left ) = masspNewInter;
             else if( conv_I_L(i_left ,j_left ) == 1)  mass(i_left ,j_left ) = masspNewInter;
             else if( conv_I_L(i_left ,j_right) == 1)  mass(i_left ,j_right) = masspNewInter;
             else if( conv_I_L(i_right,j_right) == 1)  mass(i_right,j_right) = masspNewInter;
          }
      } // end of I -> L
      //------------------------------------------------------------------//
      //------------------------------------------------------------------//
      // begin conversion: I -> G
      else if ( (fill(i,j) <= 0.0-tol) && (interface(i,j)=1.0) )
      {
           interface(i,j)	  = 0;  		// convert I -> G
           gas(i,j)	       	  = 1;
           int counter 	  	  = 0;			// counter for neighboring fluid cells, which have to be converted to interface cells

           if( liquid( i_right,j ) == 1  ) 		// find all liquid neighbor cells
           {
		liquid(   i_right,j)   	= 0;
		interface(i_right,j)	= 1;
                conv_I_G( i_right,j) 	= 1;
                counter 		= counter + 1;
	   }
	   else if (liquid(i,j_left) == 1)
	   {
		liquid(   i,j_left) 	= 0;
		interface(i,j_left)	= 1;
                conv_I_G( i,j_left) 	= 1;
  		counter 		= counter +1;
	   }
           else if (liquid(i_left,j) == 1)
           {
		liquid(	  i_left,j) 	= 0;
		interface(i_left,j)	= 1;
                conv_I_G( i_left,j) 	= 1;
		counter 		= counter + 1; 
           }
           else if (liquid(i,j_right) == 1)
	   {
		liquid(   i,j_right) 	= 0;
		interface(i,j_right)	= 1;
 		conv_I_G( i,j_right) 	= 1;
		counter 		= counter + 1; 
	   }
	   else if (liquid(i_right,j_left) == 1)		
           {
		liquid(   i_right,j_left)	= 0;
	        interface(i_right,j_left)	= 1;
		conv_I_G( i_right,j_left) 	= 1;
		counter 			= counter + 1;
	   }
           else if (liquid(i_left,j_left) == 1)
           {
		liquid(   i_left,j_left)	= 0;
		interface(i_left,j_left)	= 1;
		conv_I_G( i_left,j_left)	= 1;
		counter 			= counter + 1;
	   }
	   else if (liquid(i_left,j_right) == 1)
	   {
		liquid(   i_left,j_right)	= 0;
		interface(i_left,j_right)	= 1;
		conv_I_G( i_left,j_right) 	= 1;
		counter				= counter +1; 
	   }
	   else if (liquid(i_right,j_right) == 1)
	   {
		liquid(	  i_right,j_right)	= 0;
		interface(i_right,j_right)	= 1;
		conv_I_G( i_right,j_right)	= 1;
		counter				= counter + 1;
	   }
          if(counter != 0)
          {
             double masspNewInter = fill(i,j)/counter;
             // distribute mass to new interface cells
             // do we have to run through all cells again????
             if(      conv_I_G(i_right,j      ) == 1)  mass(i_right,j      ) = masspNewInter;
             else if( conv_I_G(i      ,j_left ) == 1)  mass(i      ,j_left ) = masspNewInter;
             else if( conv_I_G(i_left ,j      ) == 1)  mass(i_left ,j      ) = masspNewInter;
             else if( conv_I_G(i      ,j_right) == 1)  mass(i      ,j_right) = masspNewInter;
             else if( conv_I_G(i_right,j_left ) == 1)  mass(i_right,j_left ) = masspNewInter;
             else if( conv_I_G(i_left ,j_left ) == 1)  mass(i_left ,j_left ) = masspNewInter;
             else if( conv_I_G(i_left ,j_right) == 1)  mass(i_left ,j_right) = masspNewInter;
             else if( conv_I_G(i_right,j_right) == 1)  mass(i_right,j_right) = masspNewInter;
          }

      } // end of I -> G:
      //-----------------------------------------------------------------//
      //------------------------------------------------------------------//
      //  begin conversion: L -> I
      else if ( (fill(i,j)>0.0) && (fill(i,j)<1.0 - tol) && (liquid(i,j)== 1) )
      {
         interface(i,j) = 1;
         liquid(   i,j) = 0;
         // TODO: distribute mass!!!
          // to do someting????
      }// end of L -> I
      //------------------------------------------------------------------//
      //------------------------------------------------------------------//
      // begin conversion: G -> I
      else if( ( fill(i,j) >= 0 + tol) && ( gas(i,j) == 1) )
      {
	interface(i,j) = 1;
	gas(i,j)       = 0;
        
	double rho 	   = 0.0;
        double ux          = 0.0;
        double uy          = 0.0; 
        int    counter_rho = 0;

        // determine ux and uy!
        if(      liquid(i_right,j      ) == 1)
	{
		 // rho
		 double rho_cell = 0.0;
                 for(int k=1; k<=b; k++)
		 {
			rho_cell = rho_cell + fValues(i_right,j,k);
		 }
		 rho = rho + rho_cell;
                 counter_rho = counter_rho+1; 
		 // velocity
		 double ux_cell = ( ( fValues(i_right,j,2) + fValues(i_right,j,6) + fValues(i_right,j,9) ) - (fValues(i_right,j,4) + fValues(i_right,j,7) + fValues(i_right,j,8) ) )*(1.0/rho_cell);
                 double uy_cell = ( ( fValues(i_right,j,3) + fValues(i_right,j,6) + fValues(i_right,j,7) ) - (fValues(i_right,j,5) + fValues(i_right,j,8) + fValues(i_right,j,9) ) )*(1.0/rho_cell);  
			     
		 ux = ux + ux_cell;
		 uy = uy + uy_cell;
	}
        else if( liquid(i      ,j_left ) == 1) 
        {
		 // rho
		 double rho_cell = 0.0;
                 for(int k=1; k<=b; k++)
		 {
			rho_cell = rho_cell + fValues(i_right,j,k);
		 }
		 rho = rho + rho_cell;
                 counter_rho = counter_rho+1; 
		 // velocity
		 double ux_cell = ( ( fValues(i_right,j,2) + fValues(i_right,j,6) + fValues(i_right,j,9) ) - (fValues(i_right,j,4) + fValues(i_right,j,7) + fValues(i_right,j,8) ) )*(1.0/rho_cell);
                 double uy_cell = ( ( fValues(i_right,j,3) + fValues(i_right,j,6) + fValues(i_right,j,7) ) - (fValues(i_right,j,5) + fValues(i_right,j,8) + fValues(i_right,j,9) ) )*(1.0/rho_cell);  
			     
		 ux = ux + ux_cell;
		 uy = uy + uy_cell;
        }
        else if( liquid(i_left ,j      ) == 1)
	{
		 // rho
		 double rho_cell = 0.0;
                 for(int k=1; k<=b; k++)
		 {
			rho_cell = rho_cell + fValues(i_right,j,k);
		 }
		 rho = rho + rho_cell;
                 counter_rho = counter_rho+1; 
		 // velocity
		 double ux_cell = ( ( fValues(i_right,j,2) + fValues(i_right,j,6) + fValues(i_right,j,9) ) - (fValues(i_right,j,4) + fValues(i_right,j,7) + fValues(i_right,j,8) ) )*(1.0/rho_cell);
                 double uy_cell = ( ( fValues(i_right,j,3) + fValues(i_right,j,6) + fValues(i_right,j,7) ) - (fValues(i_right,j,5) + fValues(i_right,j,8) + fValues(i_right,j,9) ) )*(1.0/rho_cell);  
			     
		 ux = ux + ux_cell;
		 uy = uy + uy_cell;

	} 
        else if( liquid(i      ,j_right) == 1)  
	{
		 // rho
		 double rho_cell = 0.0;
                 for(int k=1; k<=b; k++)
		 {
			rho_cell = rho_cell + fValues(i_right,j,k);
		 }
		 rho = rho + rho_cell;
                 counter_rho = counter_rho+1; 
		 // velocity
		 double ux_cell = ( ( fValues(i_right,j,2) + fValues(i_right,j,6) + fValues(i_right,j,9) ) - (fValues(i_right,j,4) + fValues(i_right,j,7) + fValues(i_right,j,8) ) )*(1.0/rho_cell);
                 double uy_cell = ( ( fValues(i_right,j,3) + fValues(i_right,j,6) + fValues(i_right,j,7) ) - (fValues(i_right,j,5) + fValues(i_right,j,8) + fValues(i_right,j,9) ) )*(1.0/rho_cell);  
			     
		 ux = ux + ux_cell;
		 uy = uy + uy_cell;
	}
        else if( liquid(i_right,j_left ) == 1)  
	{
		 // rho
		 double rho_cell = 0.0;
                 for(int k=1; k<=b; k++)
		 {
			rho_cell = rho_cell + fValues(i_right,j,k);
		 }
		 rho = rho + rho_cell;
                 counter_rho = counter_rho+1; 
		 // velocity
		 double ux_cell = ( ( fValues(i_right,j,2) + fValues(i_right,j,6) + fValues(i_right,j,9) ) - (fValues(i_right,j,4) + fValues(i_right,j,7) + fValues(i_right,j,8) ) )*(1.0/rho_cell);
                 double uy_cell = ( ( fValues(i_right,j,3) + fValues(i_right,j,6) + fValues(i_right,j,7) ) - (fValues(i_right,j,5) + fValues(i_right,j,8) + fValues(i_right,j,9) ) )*(1.0/rho_cell);  
			     
		 ux = ux + ux_cell;
		 uy = uy + uy_cell;
	}
        else if( liquid(i_left ,j_left ) == 1)  
        {
		 // rho
		 double rho_cell = 0.0;
                 for(int k=1; k<=b; k++)
		 {
			rho_cell = rho_cell + fValues(i_right,j,k);
		 }
		 rho = rho + rho_cell;
                 counter_rho = counter_rho+1; 
		 // velocity
		 double ux_cell = ( ( fValues(i_right,j,2) + fValues(i_right,j,6) + fValues(i_right,j,9) ) - (fValues(i_right,j,4) + fValues(i_right,j,7) + fValues(i_right,j,8) ) )*(1.0/rho_cell);
                 double uy_cell = ( ( fValues(i_right,j,3) + fValues(i_right,j,6) + fValues(i_right,j,7) ) - (fValues(i_right,j,5) + fValues(i_right,j,8) + fValues(i_right,j,9) ) )*(1.0/rho_cell);  
			     
		 ux = ux + ux_cell;
		 uy = uy + uy_cell;
	}
	else if( liquid(i_left ,j_right) == 1)
	{
		 // rho
		 double rho_cell = 0.0;
                 for(int k=1; k<=b; k++)
		 {
			rho_cell = rho_cell + fValues(i_right,j,k);
		 }
		 rho = rho + rho_cell;
                 counter_rho = counter_rho+1; 
		 // velocity
		 double ux_cell = ( ( fValues(i_right,j,2) + fValues(i_right,j,6) + fValues(i_right,j,9) ) - (fValues(i_right,j,4) + fValues(i_right,j,7) + fValues(i_right,j,8) ) )*(1.0/rho_cell);
                 double uy_cell = ( ( fValues(i_right,j,3) + fValues(i_right,j,6) + fValues(i_right,j,7) ) - (fValues(i_right,j,5) + fValues(i_right,j,8) + fValues(i_right,j,9) ) )*(1.0/rho_cell);  
			     
		 ux = ux + ux_cell;
		 uy = uy + uy_cell;
	}  
        else if( liquid(i_right,j_right) == 1)  
 	{
		 // rho
		 double rho_cell = 0.0;
                 for(int k=1; k<=b; k++)
		 {
			rho_cell = rho_cell + fValues(i_right,j,k);
		 }
		 rho = rho + rho_cell;
                 counter_rho = counter_rho+1; 
		 // velocity
		 double ux_cell = ( ( fValues(i_right,j,2) + fValues(i_right,j,6) + fValues(i_right,j,9) ) - (fValues(i_right,j,4) + fValues(i_right,j,7) + fValues(i_right,j,8) ) )*(1.0/rho_cell);
                 double uy_cell = ( ( fValues(i_right,j,3) + fValues(i_right,j,6) + fValues(i_right,j,7) ) - (fValues(i_right,j,5) + fValues(i_right,j,8) + fValues(i_right,j,9) ) )*(1.0/rho_cell);  
			     
		 ux = ux + ux_cell;
		 uy = uy + uy_cell;
	}
 	// average rho, ux, uy
	rho = rho/counter_rho;
	ux  = ux/counter_rho;
	uy  = uy/counter_rho;

        // set pdf values!!!! use eq. values therefore
        for(int dir=1; dir <= b; dir++)
        {               
             int cu 			= 3.0 * ( cx(dir)*ux + cy(dir)*uy );
             fValues(i,j,dir) 	 	= weight(dir)* ( 1.0 + cu + 0.5*(cu*cu) - 1.5*( ux*ux + uy*uy ));
	}        
     } 
      //------------------------------------------------------------------//
     } // end of for
   } // end of for
   
   retval(0) = mass;
   retval(1) = fill;
   retval(2) = liquid;
   retval(3) = interface;
   retval(4) = gas; 
   retval(5) = fValues;

   return retval;//octave_value_list();

} // end of function
