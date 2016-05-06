Problem Description
    title = "Towing the Biomapper Vehicle"
    type  = towing

Analysis Parameters
    duration           = 100
    time-step          = 0.5
    dynamic-relaxation = 1.0
    dynamic-iterations = 20
               
    static-relaxation  = 1.0
    static-iterations  = 1000

    tolerance          = 1e-6
    static-solution   = catenary
    dynamic-integration = temporal

Environment
    rho            = 1025       /* no waves - just study tow dynamics */
    gravity        = 9.81
    x-current	   = -0.05 
    
Buoys
   towship        type = ship
   sled           type = sphere
                  d = 0.736
                  m = 1135
                  buoyancy = 448
                  Cdn = 0.1

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
       length   = 400
       nodes    = (401, 1.0)
       material = cable
       top-spool-pay = 1.0
       top-spool-length = 200
       top-spool-nodes = (201, 1.0)
       bottom-spool-pay = 1.0
       bottom-spool-length = 200
       bottom-spool-nodes = (201, 1.0)

   }
   terminal = {
      buoy     = towship
      x-speed  = 0
   }

End
