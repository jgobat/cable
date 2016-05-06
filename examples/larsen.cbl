Problem Description
    title = "Lazy Wave Validation Problem"
    type = general

Analysis Parameters
    static-iterations = 100000
    static-tolerance      = 0.01
    static-relaxation = 0.01
    current-steps	   = 4

    static-outer-relaxation = 200
    static-outer-iterations = 100
    static-outer-tolerance = 0.01

    dynamic-iterations	 = 20
    dynamic-tolerance    = 1e-6
    dynamic-relaxation   = 1.0
    duration		 = 120.0
    ramp-time	         = 20.0
    time-step		 = 0.01

Environment
    forcing-method = velocity
    input-type     = regular
    x-current	   = (0.0, 1.0) (172.5, 1.0) (355, 0.0)
    depth	   = 355
    z-input        = (6.0, 12.0, 0.0)
    x-input        = (4, 12.0, 1.57)
    rho		   = 1027
    gravity	   = 9.81
    bottom-stiffness = 1e6
    
Buoys
   platform	

Anchors
   clump

Materials	
   segmentA	EA = 10000e3      EI = 6.57e3        GJ = 1.0e3
                 m = 89.0         
		 d = 0.2154       Cdt = 0.05         Cdn = 1.0
   segmentB	EA = 10000e3      EI = 100.0e3       GJ = 1.0e3
                 m = 352          
		 d = 0.855        Cdt = 0.05         Cdn = 1.0


Layout
   terminal = {
      anchor = clump
   }
   segment = {
       length = 200
       material = segmentA
       nodes = (201, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 2.0
       material = segmentA
       nodes = (4, 1.0)
   }
   segment = {
       length = 1.0
       material = segmentB
       nodes = (3, 1.0)
   }
   segment = {
       length = 360
       material = segmentA
       nodes = (361, 1.0)
   }
   terminal = {
       buoy = platform
       z = 375
       x = -350
       z-force = 1.7035e5
       x-force = -3.5773e4
   }

End
