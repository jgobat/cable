Problem Description
    title = "Surface Bi-Moor"
#ifdef PHASE1
    type = positioned
# elif defined (PHASE2)
    type = general
# else
    type = surface
# endif

Analysis Parameters
    duration                = 40.0
    time-step               = 0.01
    dynamic-iterations      = 1000
    dynamic-relaxation	    = 1.0
    dynamic-tolerance       = 1e-4

    relax-stall-limit       = 75

    static-relaxation       = 0.01
    static-iterations       = 100000
    static-tolerance        = 0.001

    static-outer-iterations = 500
    static-outer-relaxation = 0.95
    static-outer-tolerance  = 0.01

Environment
    forcing-method  = velocity
    input-type      = random
    x-wave	    = (0.0, 10, 0.0)
    rho		= 1027
    gravity	= 9.81
# ifdef PHASE3
    x-current   = 0.2
# endif
    depth       = 44
    bottom-stiffness = 100
    bottom-damping = 1.0
 
Buoys
   ssar         type = axisymmetric
                diameters = (0.0, 0.235) (1.39, 0.235) (1.4, 1.27) (2.14, 1.27)
                m = 227
                Cdn = 0.5
                Cdt = 0.5
                Cdw = 1.3

Anchors
   clump	

Connectors
   shackle        m = 0.001
                  wet = 0.001
                  Cdn = 0.5
		  Cdt = 0.5
                  d = 0.001


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
       nodes = (51, 1.0)
   }
   connector = shackle
# if defined (PHASE2) || defined (PHASE3)
   branch = {
      segment = {
         length = 50
         material = wire_10mm
         nodes = (51, 1.0)
      }
      terminal = {
# ifdef PHASE2
         x-force = 20.7444
         z-force = -.27526
# endif
         anchor = clump
      }
   }
# endif
   segment = {
       length = 0.1
       material = wire_10mm 
       nodes = (3, 1.0)
   }
   terminal = {
      buoy = ssar
# ifdef PHASE1
      z = 43
      x = 20
# elif defined (PHASE2)
      z-force = 2.0*132.50095
      x-force = 0.000001
# endif
   }

End 
