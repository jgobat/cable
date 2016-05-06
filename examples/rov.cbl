Problem Description
    title = "ROV Maneuvering Example"
    type  = towing

Analysis Parameters
    static-outer-iterations = 1000
    static-iterations       = 5000
    static-relaxation       = 0.1
    static-tolerance        = 0.001
    static-outer-tolerance  = 0.01

    duration                = 1200.0
    time-step               = 0.5
    dynamic-tolerance       = 1e-6
    dynamic-relaxation      = 1.0
    dynamic-iterations      = 20
    ramp-time               = 30.0

Environment
    rho            = 1025       /* no waves - just study tow dynamics */
    gravity        = 9.81
    x-current      = -0.5
    input-type     = regular
    forcing-method = velocity
    z-input         = (1.0, 8.0, 0.0)   
   
Buoys
   towship        type = ship

   rov            type = sphere
                  d = 1.37
                  h = 1.80	     
                  m = 4000.
                  buoyancy = 4000*9.81 - 200
                  Cdn = 1.0	/* horz drag, Cd=1.0 pi*d^2/4 = 1.47 */ 
		  Cdt = 1.67	/* vert drag, Cd=1.0 d*h = 2.466     */

Connectors
   depressor      d = 1.24
                  m = 3300.
                  wet = 19620.
                  Cdn = 1.0

Materials       
cable   d = 0.0173
        m = 1.12
        wet = 8.89
        am = 0.24
        Cdn = 2.2
        Cdt = 0.01
        EA = 1.2717e7
        EI = 237.88
        GJ = 10.0

tether  d = 0.0173
        m = 0.241
        wet = 1e-6
        am = 0.241
        Cdn = 1.3
        Cdt = 0.01
        EA = 1.2717e7
        EI = 0.4
        GJ = 0.4

Layout
   terminal = {
      buoy = rov
      x-thrust = 100*4.448
   }
   segment = {
       length   = 100
       nodes    = (50, 1.0)
       material = tether
   }
   connector = depressor
   segment = {
       length   = 1265
       nodes    = (400, 1.0)
       material = cable
   }
   terminal = {
      buoy     = towship
      x-speed  = 0.0
   }

End
