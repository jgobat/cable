Problem Description
    title = "Horizontal Array Test Problem"
    type = horizontal

Analysis Parameters
    duration                = 90.0
    time-step               = 0.1
    dynamic-iterations      = 2000
    dynamic-relaxation	    = 0.1
    dynamic-tolerance       = 1e-6

    static-relaxation       = 0.5
    static-iterations       = 30000
    static-tolerance        = 1e-8
    static-outer-iterations = 500
    static-outer-relaxation = 8
    static-outer-tolerance  = 0.001
    current-steps           = 0

Environment
    forcing-method  = morison
    input-type      = random
    x-wave	    = (4.5, 13, 0.0)
    rho		= 1027
    gravity	= 9.81
    x-current   = 0.2
    depth       = 44
   
Buoys
   sinker	  d = 6*0.0254
		  m = 5
		  Cdn = 1.0
		  buoyancy = 10
 
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
      x = -30
      y = -30
      z = 0.0
   }
   segment = {
       length = 50
       material = wire_10mm 
       nodes = (10, 1.0)
   }
   connector = sphere_48in
   branch = {
      segment = {
         length = 50
         material = wire_10mm
         nodes = (11, 1.0)
      }
      terminal = {
 	 anchor = clump
         x = 30.0
	 z = 0.0
         y = -30
      }
   }
   branch = {
      segment = {
         length = 50
         material = wire_10mm
         nodes = (10, 1.0)
      }
      terminal = {
 	 anchor = clump
         x = -30.0
	 z = 0.0
         y = 30
      }
   }
   segment = {
       length = 50.0
       material = wire_10mm
       nodes = (10, 1.0)
   }
   terminal = {
      anchor = clump
      x = 30
      y = 30
      z = 0.0
   }

End 
