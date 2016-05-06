Problem Description
    title = "Tethered Mine Example Problem"
    type = deployment

Analysis Parameters
    duration                = 1000.0
    time-step               = 0.2
    dynamic-iterations      = 20
    dynamic-relaxation      = 1.0
    dynamic-tolerance       = 1e-8

    static-relaxation       = 0.2
    static-iterations       = 1000
    static-tolerance        = 1e-6
    static-initial-guess    = shooting
    shooting-iterations     = 100

    static-outer-relaxation = 0.98
    static-outer-iterations = 500
    static-outer-tolerance  = 0.01
    static-solution   = relaxation

Environment
    rho            = 1027
    gravity        = 9.81
    depth          = 1000
    bottom-stiffness = 250    
 
Buoys
   mine         type = sphere
                d = 1.5           m = 1200
                Cdn = 0.3         Cdt = 0.3 
                
Anchors
   clump        color = red
                wet = 2000*4.448 
                m = 2000/0.86/2.2
                d = 0.81
                Cdn = 1.0
                Cdt = 1.0


Materials       
   wire_10mm   EA = 4.7e6           EI = 20.0         GJ = 5
                m = 0.320           am = 0.08       wet = 2.64
                d = 0.01           Cdt = 0.01       Cdn = 1.5
                type = linear



Layout
   terminal = {
      anchor = clump
      x-speed = 2.0*0.514
   }
   segment = {
       length = 1050.0
       material = wire_10mm
       nodes = (201, 1.0)
   }
   terminal = {
       buoy = mine
   }

End
