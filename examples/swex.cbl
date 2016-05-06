Problem Description
    title = "SWEX Buoy Farm Surface Mooring: Baseline case"
    type  = surface

Analysis Parameters
    static-relaxation = 0.01
    static-iterations = 20000
    static-tolerance  = 0.001

    relax-adapt-up    = 1.02
    relax-adapt-down  = 1.1

    static-outer-iterations = 100
    static-initial-guess    = shooting
    static-solution         = shooting
    static-outer-relaxation = 0.98
    static-outer-tolerance  = 0.01
    static-solution   = catenary

    duration           = 210.0
    ramp-time	       = 10.0
    time-step          = 0.05    /* 0.1 */
    dynamic-relaxation = 1.0
    dynamic-iterations = 15
    dynamic-tolerance  = 1e-6

Environment
    input-type       = random
    x-current	     = (0.0, 0.534) (0.32, 0.534) (0.57, 0.388) (0.82, 0.313) 
			(1.07, 0.296) (1.32, 0.257) (1.57, 0.212) (1.82, 0.194)
		 	(2.07, 0.192) (2.32, 0.192) (2.57, 0.189) (2.82, 0.189) 			(3.07, 0.187) (3.32, 0.193) (3.57, 0.193) (3.82, 0.191) 			(4.07, 0.192) (4.32, 0.187) (4.57, 0.190) (4.82, 0.191) 			(5.07, 0.189) (5.32, 0.193) (5.57, 0.191) (5.82, 0.191) 			(6.07, 0.190) (6.32, 0.190) (6.57, 0.190) (6.82, 0.185) 			(7.07, 0.185) (7.32, 0.184) (7.57, 0.184) (7.82, 0.185) 			(8.07, 0.188) (8.32, 0.185) (8.57, 0.184) (8.82, 0.186) 			(9.07, 0.186) (9.32, 0.189) (9.57, 0.185) (9.82, 0.186) 			(10.07, 0.182) (42.60, 0.182)
    x-wind	     = 7.17
    rho		     = 1027
    gravity	     = 9.81
    depth            = 43.31
    bottom-stiffness = 100.0
    bottom-damping   = 1.0
    forcing-method   = velocity
    velocity-file    = "swex_6dec98.vel"
    
Buoys
   ssar		type = axisymmetric
 		diameters = (0.0, 0.235) (1.39, 0.235) (1.4, 1.27) (2.14, 1.27)
		m = 227
		Cdn = 0.5
		Cdt = 0.5
                Cdw = 1.3

Anchors
   clump	color = red

Materials	
chain_1_2       EA = 6.44e+07    EI = 10.0          GJ = 0.1
                 m = 3.7335        		   wet = 31.8455
                 d = 0.0495      Cdt = 0.01        Cdn = 0.5
				 Cmt = 1.1        Cmn = 2.0
                 comment = "1/2 Trawler Chain"

	/*
	 * revised AxPack Cmt, Cmn from 1.5, 2.0 
	 * after coeff3.dat
	 */

AxPack		EA = 8.0e7	 EI = 30000.0       GJ = 0.1
		 m = 10.02   			   wet = 70.82
		 d = 0.0762	Cdt = 0.065   	   Cdn = 0.8
				 Cmt = 1.1         Cmn = 2.0
		 comment = "Length = 0.76 m, with SRS at both ends"

Connectors
shackle		d = 0.01	/* just something to release moments */
		wet = 0.01
		m = 0.01

Layout
   terminal = {
      anchor = clump
   }
   segment = {
       length = 45
       material = chain_1_2
       nodes = (91, 1.0)
   }
   connector = shackle
   segment = {
       length = 0.76
       material = AxPack
       nodes = (3, 1.0)
   }
   connector = shackle
   segment = {
       length = 3.5
       material = chain_1_2
       nodes = (9, 1.0)
   }
   connector = shackle
   segment = {
       length = 0.76
       material = AxPack
       nodes = (3, 1.0)
   }
   connector = shackle
   segment = {
       length = 7.0
       material = chain_1_2
       nodes = (15, 1.0)
   }
   connector = shackle
   segment = {
       length = 0.76
       material = AxPack
       nodes = (3, 1.0)
   }
   connector = shackle
   segment = {
       length = 22.5 + 0.5	/* 0.5 for srs-swivel-srs at top */
       material = chain_1_2
       nodes = (24, 1.0)
   } 
   terminal = {
      buoy = ssar
   }

End
