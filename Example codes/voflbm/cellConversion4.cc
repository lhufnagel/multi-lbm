#include <octave/oct.h>

// TODO: check boundary because of periodicity!!!! -> change cell conversions!!!

DEFUN_DLD( cellConversion4,args,, "do necessary cell conversions")
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
   // possible cell conversions   
   for(int i = 0; i < ly; i++){
      for( int j=0; j< lx; j++){

      // I -> L
      // TODO: conv_I_L = zeros(ly,lx) -> otherwise old values!!!!
      if ( (fill(i,j)>= 1.0 + tol) && (interface(i,j)== 1) )  	// convert cell to liquid cell 
      {
         std::cout << " I -> L begin!!! at (i= " << i << "," << "j= " << j << ")" << std::endl; 
         interface(i,j) = 0;
         liquid(i,j)    = 1;       	   			// counter for gas cells in neighborhood ->  how we can distribute mass
         int counter    = 0;
       
         for(int k=1; k<b; k++)
   	 {
            std::cout << "Search gas neighbors: k = " << k  << std::endl;         
           int ni = ( (i+(int)cy(k) ) % ly);
           int nj = ( (j+(int)cx(k) ) % lx);
           if(ni < 0){
	      ni = ly - (abs( i+(int)cy(k) ) % ly);
           }
	   if(nj < 0){
	      nj = lx - (abs( j+(int)cx(k) ) % lx);
           }
           std::cout << "Gas neighbors k= " << k << "( ni= " << ni << ", nj= "  << nj << ")" << std::endl;
		
   	   if( (gas(ni,nj) == 1) || (conv_I_L(ni,nj) ==1) )	// search for gas cells! -> convert them to interface cells
   	   {
             interface(ni,nj) 		= 1;            	// generate new interface cells
             conv_I_L(ni,nj)	 	= 1;	 		// inter-Flag: remember new interface flags to know which interface flags get mass!
             counter 			= counter +1;	
	   }
           std::cout << "Counter of gas neighbors: " << counter << std::endl;
   	 } // end of neighbors 

         if(counter != 0)     
	 {
             double masspNewInter = (mass(i,j)-1.0)/counter;
             mass(i,j) = 1.0;
	     
   	     for(int k=1; k<b; k++)
   	     {
   	        int ni = ((i+(int)cy(k)) % ly) ;
   	        int nj = ((j+(int)cx(k)) % lx) ;
                if(ni < 0){
                	ni = ly - (abs( i+(int)cy(k) ) % ly);
               	}
                if(nj < 0){
                	nj = lx - (abs( j+(int)cx(k) ) % lx);
               	}
               	// distribute mass to new interface cells
   	       	if( conv_I_L(ni,nj) == 1)
   	       	{
			mass(ni,nj) = mass(ni,nj) + masspNewInter;
			gas(ni,nj)  = 0;
                        std::cout << "Distribute masspNewInter: " << masspNewInter << " to ( ni= " << ni << "nj=  " << nj << ")" << std::endl;
			std::cout << " mass("<< ni <<" , " << nj << ") = " << mass(ni,nj) << std::endl; 		   
   		}
   	      }
         }
         /////////////////////////////////////////////////////////////////////////
  	 std::cout << "G -> I begin as a consequence of I -> L !!! " << std::endl;
         std::cout << "start initialize process for fValues    !!! " << std::endl;
         
         double rho         = 0.0;
         double ux          = 0.0;
         double uy          = 0.0;
         int    counter_rho = 0;

         // determine ux and uy of liquid neighbors
         for(int k=1; k<b; k++)
         {
          	int ni = ((i+(int)cy(k)) % ly) ;
          	int nj = ((j+(int)cx(k)) % lx) ;
          	if(ni < 0){
              		ni = ly - (abs( i+(int)cy(k) ) % ly);
          	}
          	if(nj < 0){
              		nj = lx - (abs( j+(int)cx(k) ) % lx);
          	}	

          	if(liquid(ni,nj) == 1)//average values of liquid neighbor cells!
          	{
			std::cout << "liquid neighbor: ni= " << ni << ", nj =" << nj << std::endl;
			// rho
                 	double rho_cell = 0.0;
                 	for(int k=0; k< b; k++)
                 	{
                        	rho_cell = rho_cell + fValues(ni,nj,k);
				//std::cout << ni << ", " << nj << "rho_cell " << rho_cell << std::endl;
                 	}
                 	rho = rho + rho_cell;
			//std::cout << "rho = " << rho << std::endl;
                 	counter_rho = counter_rho+1;
                 	// velocity
                 	double ux_cell = ( ( fValues(ni,nj,1) + fValues(ni,nj,5) + fValues(ni,nj,8) ) - (fValues(ni,nj,3) + fValues(ni,nj,6) + fValues(ni,nj,7) ) )*(1.0/rho_cell);
                 	double uy_cell = ( ( fValues(ni,nj,2) + fValues(ni,nj,5) + fValues(ni,nj,6) ) - (fValues(ni,nj,4) + fValues(ni,nj,7) + fValues(ni,nj,8) ) )*(1.0/rho_cell);

                 	ux = ux + ux_cell;
                 	uy = uy + uy_cell;
			//std::cout << " ux = " << ux << std::endl;
			//std::cout << " uy = " << uy << std::endl;
          	}
          }
          // average rho, ux, uy
          // std::cout << " rho = " << rho << "; counter_rho= " << counter_rho << std::endl;
          // std::cout << " ux = "  << ux  << std::endl;
          // std::cout << " uy = "  << uy  << std::endl; 
          rho = rho/counter_rho;
          ux  = ux/counter_rho;
          uy  = uy/counter_rho;

	  std::cout << " avg rho = " << rho << std::endl;
          std::cout << " avg ux  = " << ux  << std::endl;
          std::cout << " avg uy  = " << uy  << std::endl;

          // set pdf values!!!! use eq. values therefore
          for(int k=1; k<b; k++)
          {
          	// std::cout << "Search ex gas neighbors  << std::endl;         
           	int ni = ( (i+(int)cy(k) ) % ly);
           	int nj = ( (j+(int)cx(k) ) % lx);
           	if(ni < 0){
              		ni = ly - (abs( i+(int)cy(k) ) % ly);
           	}
           	if(nj < 0){
              		nj = lx - (abs( j+(int)cx(k) ) % lx);
           	}
           	
           	if( conv_I_L(ni,nj) ==1 )     // search for ex gas cells! -> initialize them!!!
                {
			//std::cout << "conv_I_L at (ni = " << ni << ", nj= " << nj << ")"  << std::endl;
 			double sum_fValues = 0.0; 
	         	for(int dir = 0; dir < b; dir++)
          		{
                                //std::cout << "dir = " << dir << ", cy(dir)*uy = " << cy(dir)*uy*3.0 << std::endl;
				//std::cout << "dir = " << dir << ", cx(dir)*ux = " << cx(dir)*ux*3.0 << std::endl; 
             			double cu               = 3.0 * ( cx(dir)*ux + cy(dir)*uy );
				//std::cout << "cu = " << cu << std::endl;
             			fValues(ni,nj,dir)   = weight(dir)* ( 1.0 + cu + 0.5*(cu*cu) - 1.5*( ux*ux + uy*uy ));
				sum_fValues          = sum_fValues + fValues(ni,nj,dir);
				// std::cout << "fValues(i = " << i << ", j " << j << ", dir " << dir << ") -> fValues = " << fValues << std::endl;  
          		}
			std::cout << "check sum_fValues = " << sum_fValues << std::endl;
		}
         } 
	 std::cout << " end of conversion G -> I !!! " << std::endl;
         std::cout << " end of conversion I -> L !!! " << std::endl;
      } // end of if for I -> L
      //------------------------------------------------------------------//
      //------------------------------------------------------------------//
      // begin conversion: I -> G
      else if ( (fill(i,j) <= 0.0-tol) && (interface(i,j)=1) )
      {
	std::cout << "I -> G begin!!!, i:= " << i << "j:= " << j << std::endl;
        interface(i,j)	  = 0;  		// convert I -> G
        gas(i,j)	  = 1;
        int counter 	  = 0;			// counter for neighboring fluid cells, which have to be converted to interface cells

        for(int k=1; k<b; k++)
   	{
   		int ni = ((i+(int)cy(k)) % ly);
   	   	int nj = ((j+(int)cx(k)) % lx);		

          	if(ni < 0){
              		ni = ly - (abs( i+(int)cy(k) ) % ly);
          	}
          	if(nj < 0){
              		nj = lx - (abs( j+(int)cx(k) ) % lx);
          	}
	  	//std::cout << "Liquid neighbors: ni: =" << ni << "nj:= " << nj << std::endl;
   	  	if( (liquid(ni,nj) == 1) || (conv_I_G(ni,nj)==1) ) 		// find all liquid neighbor cells
          	{
   	   		//liquid(ni,nj)   	= 0;
   			interface(ni,nj)	= 1;
         		conv_I_G(ni,nj) 	= 1;
           		counter 		= counter + 1;
                }
   	}
	std::cout << "Counter of liquid neighbors = " << counter << std::endl;         
        double mass_new = mass(i,j);
        fill(i,j) = 0.0;
        mass(i,j) = 0.0;

	if(counter != 0)
        {
        	double masspNewInter = mass_new/counter;
             	//std::cout << "massNewInter= " << masspNewInter << std::endl;
             	// distribute mass to new interface cells
             	for(int k=1; k<b; k++)
   	     	{
   	       		int ni = ((i+(int)cy(k)) % ly) ;
   	       		int nj = ((j+(int)cx(k)) % lx) ;
		        if(ni < 0){
                		ni = ly - (abs( i+(int)cy(k) ) % ly);
               		}
               		if(nj < 0){
                  		nj = lx - (abs( j+(int)cx(k) ) % lx);
               		}
               		std::cout << "Liquid neighbors for mass: ni= " << ni << "nj= " << nj << std::endl;
               		if( conv_I_G(ni,nj) == 1) 
   	       		{	 
   	   	  		mass(ni,nj) = mass(ni,nj)+masspNewInter; 
		  		liquid(ni,nj) = 0;
	          		std::cout << "mass(ni,nj)= " << mass(ni,nj)  << std::endl;
   	       		}
   	     	}	
        }
        std::cout << "end of conversion I -> G!" << std::endl;;
      } // end of I -> G:
      //-----------------------------------------------------------------//
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

