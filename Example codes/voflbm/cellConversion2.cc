#include <octave/oct.h>

// TODO: check boundary because of periodicity!!!! -> change cell conversions!!!

DEFUN_DLD( cellConversion2,args,, "do necessary cell conversions")
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
   NDArray fValues      = args(9).array_value();
   int    b             = args(10).int_value();
   RowVector weight     = args(11).row_vector_value();
   RowVector cx         = args(12).row_vector_value();
   RowVector cy         = args(13).row_vector_value();  
   
   octave_value_list retval;

   const double tol = 0.01; 

   //---------------------------------------------------------------------//
   // 4 possible cell conversions   
   for(int i = 0; i < ly; i++){
      for( int j=0; j< lx; j++){

      // I -> L
      // TODO: conv_I_L = zeros(ly,lx) -> otherwise old values!!!!
      if ( (fill(i,j)>= 1.0 + tol) && (interface(i,j)== 1) )  	// convert cell to liquid cell 
      {
         interface(i,j) = 0;
         liquid(i,j)    = 1;       	   			// counter for gas cells in neighborhood -> 
								// how we can distribute mass
         int counter    = 0;

         std::cout << "I-> L begin" <<" " << "i: " << i << " j: " << j << std::endl;	
         for(int k=1; k<b; k++)
   	 {

           std::cout << "Search gas neighbors" << std::endl;
           
           /// adjust mod to get positive values
           //int ni = ((i+(int)cy(k)-1) % ly);
   	   //int nj = ((j+(int)cx(k)-1) % lx);
           int ni = ( (i+(int)cy(k) -1 ) % ly);
           int nj = ( (j+(int)cx(k) -1 ) % lx);
           if(ni < 0){
	      ni = ly - (abs( i+(int)cy(k)-1 ) % ly);
           }
	   if(nj < 0){
              nj = lx - (abs( j+(int)cx(i)-1 ) % lx  );
           }
           std::cout << "ni,nj: " << "(" << ni << "," << nj << ")" << std::endl;
		
   		if(gas(ni,nj) == 1)				// search for gas cells! -> convert them to interface cells
   		{
        	    gas(ni,nj)	 		= 0;   	 	// generate new interface cell!
                    interface(ni,nj) 		= 1;
                    std::cout << " generate new interface cells" << std::endl;
   		    conv_I_L(ni,nj)	 	= 1;	 	// inter-Flag: remember new interface flags to know which interface flags get mass!
                    counter 			= counter +1;
		
   		}
               std::cout << "counter: " << counter << std::endl;
   	 } 

         if(counter != 0)
          {
             double masspNewInter = (fill(i,j)-1.0)/counter;
   	     for(int k=1; k<b; k++)
   	     {
   	       int ni = ((i+(int)cy(k)-1) % ly) ;
   	       int nj = ((j+(int)cx(k)-1) % lx) ;
		
                // distribute mass to new interface cells
   		if( conv_I_L(ni,nj) == 1)
   		{
   			mass(ni,nj) = masspNewInter;
   		}

   	     }
          }
      } // end of I -> L
      //------------------------------------------------------------------//
      //------------------------------------------------------------------//
      // begin conversion: I -> G
      else if ( (fill(i,j) <= 0.0-tol) && (interface(i,j)=1) )
      {
           interface(i,j)	  = 0;  		// convert I -> G
           gas(i,j)	       	  = 1;
           int counter 	  	  = 0;			// counter for neighboring fluid cells, which have to be converted to interface cells

         for(int k=1; k<b; k++)
   	 {
   	   int ni = ((i+(int)cy(k)-1) % ly);
   	   int nj = ((j+(int)cx(k)-1) % lx);		

   			if( liquid(ni,nj) == 1  ) 		// find all liquid neighbor cells
           		{
   				liquid(ni,nj)   	= 0;
   				interface(ni,nj)	= 1;
                		conv_I_G(ni,nj) 	= 1;
                		counter 		= counter + 1;
   	   		}
   	 }
	   
         if(counter != 0)
         {
             double masspNewInter = fill(i,j)/counter;
             // distribute mass to new interface cells
             for(int k=1; k<b; k++)
   	     {
   	       int ni = ((i+(int)cy(k)-1) % ly) ;
   	       int nj = ((j+(int)cx(k)-1) % lx) ;
		
                 if( conv_I_G(ni,nj) == 1) 
   		 { 
   			mass(ni,nj) = mass(ni,nj)+masspNewInter; // changed..............
   		 }

   	     }	
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
      else if( ( fill(i,j) > 0.0 + tol) && ( gas(i,j) == 1) )
      {
   	interface(i,j) = 1;
   	gas(i,j)       = 0;
        
   	double rho 	   = 0.0;
        double ux          = 0.0;
        double uy          = 0.0; 
        int    counter_rho = 0;

        // determine ux and uy!
        for(int k=1; k<b; k++)
        {
   	  int ni = ((i+(int)cy(k)-1) % ly) ;
   	  int nj = ((j+(int)cx(k)-1) % lx) ;

	
        	if(liquid(ni,nj) == 1)
   		{
   			 // rho
   		 	double rho_cell = 0.0;
                 	for(int k=1; k<b; k++)
   		 	{
   				rho_cell = rho_cell + fValues(ni,nj,k);
   		 	}
   		 	rho = rho + rho_cell;
                 	counter_rho = counter_rho+1; 
   		 	// velocity
   		 	double ux_cell = ( ( fValues(ni,nj,1) + fValues(ni,nj,5) + fValues(ni,nj,8) ) - (fValues(ni,nj,3) + fValues(ni,nj,6) + fValues(ni,nj,7) ) )*(1.0/rho_cell);
                 	double uy_cell = ( ( fValues(ni,nj,2) + fValues(ni,nj,5) + fValues(ni,nj,6) ) - (fValues(ni,nj,4) + fValues(ni,nj,7) + fValues(ni,nj,8) ) )*(1.0/rho_cell);  
			     
   		 	ux = ux + ux_cell;
   		 	uy = uy + uy_cell;
   		}
   	}

   	// average rho, ux, uy
   	rho = rho/counter_rho;
   	ux  = ux/counter_rho;
   	uy  = uy/counter_rho;

        // set pdf values!!!! use eq. values therefore
        for(int dir=1; dir < b; dir++)
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
