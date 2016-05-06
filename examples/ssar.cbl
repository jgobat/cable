Problem Description
    title = "SSAR Baseline Drifter w/Extra Sinker Weights - Sept 95 Sea Trials"
    type = drifter

Analysis Parameters
    duration       = 118.0
    time-step      = 0.1
    dynamic-relaxation = 1.0
    dynamic-tolerance  = 1e-6
    dynamic-iterations = 20

    static-relaxation  = 1.0
    static-iterations  = 1000
    static-tolerance   = 1e-6
    static-outer-iterations = 25
    static-outer-tolerance  = 0.001
    static-solution   = catenary

Environment
    forcing-method = velocity
    input-type	   = random
    x-current      = (0.0, 0.2) (100.0, 0.1) (500.0, 0.0) (800.0, 0.0)
    y-current      = (0.0, -0.1) (100.0, -0.1) (500.0, 0.0) (800.0, 0.0)
    rho		   = 1027
    gravity	   = 9.81
    depth          = 800
    velocity-file  = "ssar.vel"
      
Buoys
   snubber        type = cylinder h  = 0.66
		  d = 1.07	  m = 227
		  Cdn = 1.0       Cdt = 1.2
		  Cdw = 0.6
   sinker	  type = cylinder h = 0.9
                  d = 0.144	  m = 321
                  Cdn = 1.0	  Cdt = 1.0
		  buoyancy = 480.38

Connectors
   ssf		  d = 0.76	  m = 50
                  Cdn = 0.8	  Cdt = 0.8
                  wet = -2000 
   shackle	  d = 0.08	  m = 2
		  Cdn = 1.0       Cdt = 1.0
   e_pack	  d = 0.237  	  m = 148
                  Cdn = 1.0      Cdt = 1.0
                  wet = 450

Materials	
wire_rope d=0.008    EA=4.2e6    m=0.193     Cdt=0.004   Cdn=1.5  EI=2.0  GJ=1.0
array     d=0.015    EA=1.34e6   m=0.243     Cdt=0.05    Cdn=1.5  EI=10.0 GJ=1.0
EM_cable  d=0.0131   EA=1.5e6    m=0.320     Cdt=0.008   Cdn=1.5  EI=10.0 GJ=1.0
hose      d=0.104    EA=3.4e4    m=9.61      Cdt=0.03    Cdn=1.0  EI=100 GJ=1.0

Anchors
clump

Layout
   terminal = {
      buoy = sinker
   }
   segment = {
       length = 10
       material = wire_rope
       nodes = (5, 1)
   }
   connector = shackle
   segment = {
       length = 60
       material = array
       nodes = (20, 1)
   }
   connector = e_pack
   segment = {
       length = 500
       material = EM_cable
       nodes = (100, 1)
   }
   connector = ssf
   segment = {
       length = 27
       material = hose
       nodes = (27, 1)
   }
   terminal = {
      buoy = snubber
   }

End
