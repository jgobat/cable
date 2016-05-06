Problem Description
    title = "Lowering an Anchor"
    type  = towing

Analysis Parameters
    duration           = 300
    time-step          = 0.05
    dynamic-relaxation = 1.0
    dynamic-iterations = 20
               
    static-relaxation  = 1.0
    static-iterations  = 1000

    tolerance          = 1e-12

Environment
    rho            = 1025       /* no waves - just study tow dynamics */
    gravity        = 9.81
    depth	   = 40.0
    bottom-stiffness = 10000.0
    bottom-friction  = 0.5
    bottom-elevation = 2.0*sin(pi*x/50)
    
Buoys
   towship        type = ship
   sled           type = sphere
                  d = 0.736
                  m = 1135
                  buoyancy = 4448
                  Cdn = 0.77

Materials       
cable   d = 0.0173
        m = 1.12
        wet = 8.89
        am = 0.24
        Cdn = 1.5
        Cdt = 0.01
        EA = 1.2717e7
        EI = 237.88
        GJ = 10.0

Layout
   terminal = {
      buoy = sled
   }
   segment = {
       length   = 10
       nodes    = (11, 1.0)
       material = cable
   }
   terminal = {
      buoy     = towship
      pay-rate = 0.6 
      x-speed  = 0.5
   }

End
