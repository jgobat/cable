Problem Description
    title = "Breaking Tow Cable Demonstration"
    type  = towing

Analysis Parameters
    duration           = 400.0
    time-step          = 0.5
    relaxation         = 1.0
    static-iterations  = 1000
    dynamic-iterations = 50
    tolerance          = 1e-6
    dynamic-integration = temporal

Environment
    rho		   = 1000
    gravity	   = 10.0
    
Buoys
   tow_ship       type = ship 
   rov    	  type = sphere   h = 0.5
                  d = 1.0	  m = 500
		  buoyancy = 4500
                  Cdn = 2.0	  Cdt = 2.0

Materials	
tow_cable  d=0.0131   EA=1.5e6    EI=10.0
           m=0.3      wet=2.0
           Cdt=0.005  Cdn=1.0     GJ=1.0

Layout
   terminal = {
      buoy = rov
   }
   segment = {
       length = 1000
       material = tow_cable
       nodes = (200, 1.0)
   }
   terminal = {
      buoy = tow_ship
      x-speed = 0.5
      release-time = 100
   }

End
