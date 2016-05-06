Problem Description
    title = "Shallow Water S-tether Mooring"
    type = surface

Analysis Parameters
    static-outer-iterations = 1000
    static-iterations       = 5000
    static-relaxation       = 0.1
    static-tolerance        = 0.01
    static-outer-tolerance  = 0.01
    static-outer-relaxation = 0.5
    static-solution   = catenary
    duration                = 400.0
    time-step               = 0.05
    dynamic-tolerance       = 1e-6
    dynamic-relaxation      = 1.0
    dynamic-iterations      = 20
    ramp-time               = 50.0

Environment
    rho            = 1025
    gravity        = 9.81
    x-current      = 0.12
    depth          = 100
    
Buoys
   float        type = sphere
                m = 1
                d = 1.0
                Cdn = 1.0

Anchors
   clump        color = red

Connectors
   Panther	d = 0.36
		wet = -200
		Cdn = 0.95
		m = 0.32

Materials       
   Umbilical    EA = 3.4e7           EI = 20.0         GJ = 25
                 m = 1.383           
                 d = 0.0221          Cdt = 0.01       Cdn = 1.2

Layout
   terminal = {
       anchor = clump
   }
   segment = {
       length = 150
       material = Umbilical
       nodes = (150, 1.0)
       attachments = Panther : (147, 144, 141, 138, 135, 132, 72, 69, 66, 63, 60, 57, 54, 51)
   }
   terminal = {
       buoy = float
   }

End

