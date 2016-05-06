Problem Description
    title = "Horizontal Array Test Problem"
    type = horizontal

Analysis Parameters
    static-tolerance        = 1e-4
    static-relaxation       = 0.01
    static-iterations       = 100000

    relax-adapt-up = 1.0
    relax-adapt-down = 1.0

    static-outer-tolerance  = 0.01
    static-outer-iterations = 400
    static-outer-relaxation = 15
    duration                = 100.0
    time-step               = 0.05
    dynamic-tolerance       = 1e-6
    dynamic-iterations      = 500
    dynamic-relaxation	    = 0.5
    ramp-time		    = 0.0
    dynamic-integration     = spatial

Environment
    forcing-method  = morison
    input-type      = random
    x-wave	    = (1.0, 10, 0.0)
    y-wave	    = (0.4, 7.5, 0.75)
    rho		= 1027
    gravity	= 9.81
    x-current   = 0.1
    y-current   = 0.05
    depth       = 45
   
Buoys
   sinker	  d = 6*0.0254
		  m = 20
		  Cdn = 1.0
		  type = sphere
 
Anchors
   clump	

Connectors
   sphere_48in    d = 48*0.0254
                  wet = -1590*4.448 - 6000
                  Cdn = 0.8
                  m = 690/2.2
   shackle        m = 4.45/2.205
                  wet = 4.45*4.448*.8667
                  Cdn = 0.5
                  d = 0.168


Materials
   wire_10mm        EA = 4.7e6           EI = 20          GJ = 5
                     m = 0.320           am = 0.08       wet = 2.64
                     d = 0.01           Cdt = 0.01       Cdn = 1.5


Layout
   terminal = {
      anchor = clump
   }
   segment = {
       length = 50
       material = wire_10mm 
       nodes = (10, 1.0)
   }
   connector = sphere_48in
   segment = {
       length = 20
       material = wire_10mm
       nodes = (10, 1.0)
   }
   connector = shackle
   branch = {
      segment = { length = 40
         material = wire_10mm
         nodes = (20, 1.0) }
      terminal = {
 	 node = 50
      }
   }
   segment = { 
      length = 20
      material = wire_10mm
      nodes = (10, 1.0) 	
   }
   connector = shackle
   segment = { 
      length = 20 
      material = wire_10mm
      nodes = (10, 1.0) 
   }
   connector = shackle
   branch = {
      segment = { 
	 length = 40
         material = wire_10mm
         nodes = (20, 1.0) 
      }
      terminal = {
 	 node = 90
      }
   }
   segment = {
      length = 20
      material = wire_10mm
      nodes = (10, 1.0)
   }
   connector = shackle
   segment = {
      length = 20
      material = wire_10mm
      nodes = (10, 1.0)
   }
   connector = sphere_48in
   segment = {
       length = 50.0
       material = wire_10mm
       nodes = (10, 1.0)
   }
   terminal = {
      anchor = clump
      x = 160
      z = 0.0
   }

End 
