Problem Description
    title = "Horizontal Array Test Problem"
    type = horizontal

Analysis Parameters
    tolerance               = 1e-6

    static-relaxation       = 1.0
    static-iterations       = 30000
    static-outer-iterations = 2000
    static-outer-relaxation = 20
    static-outer-tolerance  = 0.001
    current-steps           = 0
    duration                = 100.0
    time-step               = 0.01
    dynamic-iterations      = 20
    dynamic-relaxation	    = 1.0
    ramp-time		    = 10.0

Environment
    forcing-method  = morison
    input-type      = random
    x-wave	    = (2.0, 10, 0.0)
    y-wave          = (0.4, 6.5, 0.75)
    rho		= 1027
    gravity	= 9.81
    x-current   = 0.2
    y-current   = 0.05
    depth       = 70
   
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
       nodes = (5, 1.0)
   }
   connector = sphere_48in
   segment = {
       length = 20
       material = wire_10mm
       nodes = (5, 1.0)
   }
   connector = shackle
   branch = {
      segment = {
         length = 20
         material = wire_10mm
         nodes = (5, 1.0)
      }
      connector = shackle
      segment = {
         length = 10
         material = wire_10mm
         nodes = (5, 1.0)
      }
      terminal = {
         buoy = sinker
      }
   }
   segment = {
      length = 20
      material = wire_10mm
      nodes = (5, 1.0)
   }
   connector = shackle
   branch = {
      segment = {
         length = 20
         material = wire_10mm
         nodes = (5, 1.0)
      }
      terminal = {
 	 buoy = sinker
      }
   }
   segment = {
      length = 20
      material = wire_10mm
      nodes = (5, 1.0)
   }
   connector = sphere_48in
   segment = {
       length = 50.0
       material = wire_10mm
       nodes = (5, 1.0)
   }
   terminal = {
      anchor = clump
      x = 120
      z = 0.0
   }

End 
