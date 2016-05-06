Problem Description
    title = "Towing the Biomapper Vehicle in a Circle"
    type  = towing

Analysis Parameters
    duration           = 1200
    time-step          = 1.0
    dynamic-relaxation = 1.0
    dynamic-iterations = 20
               
    static-relaxation  = 0.3
    static-iterations  = 400

    tolerance          = 1e-6
    static-solution   = catenary

Environment
    rho            = 1025       /* no waves - just study tow dynamics */
    gravity        = 9.81
    
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
        GJ = 100

Layout
   terminal = {
      buoy = sled
   }
   segment = {
       length   = 1000
       nodes    = (100, 1.0)
       material = cable
   }
   terminal = {
      buoy     = towship
      x-speed  = t < 200 || t > 800 ? 0.514*6 : 0.514*6*cos(2*3.14159/600*(t - 200))
      y-speed  = t < 200 || t > 800 ? 0.0     : 0.514*6*sin(2*3.14159/600*(t - 200))
   }

End
