Problem Description
    title = "AOSN Labrador Sea Mooring"   
    type = surface

Analysis Parameters
    static-outer-iterations = 1000
    static-outer-relaxation = 0.75
    static-iterations       = 5000
    static-relaxation       = 0.1
    static-tolerance        = 1e-6
    static-outer-tolerance  = 0.01

    duration                = 200.0
    time-step               = 0.2
    dynamic-tolerance       = 1e-6
    dynamic-relaxation      = 1.0
    dynamic-iterations      = 20
    ramp-time               = 50.0

Environment
    rho            = 1025
    gravity        = 9.81
    x-current      = (0.0, 0.1) (100, 0.1) (1000, 0.05) (3500, 0.05)
    depth          = 3500
    input-type     = random
    forcing-method = wave-follower
    x-wave         = (2.7, 10.0, 0.0)   /* sea state 6 */
# if 0
    x-current-modulation = t < 7200 ? cos(2.0*3.14159/14400*t) : -1
# endif
    bottom-stiffness = 0.0
    bottom-damping   = 0.0
    
Buoys
   float        type = sphere
                m = 65
                d = 1.22
                Cdn = 1.0


Anchors
   clump        color = red

Connectors
   shackle      d = 0.375*0.0254
                wet = 8.0
                Cdn = 0.1
                m = 1.0
   ssf          d = 1.628
                m = 931.3
                wet = -13350
                Cdn = 2.068

Materials       
   wire         EA = 4.4e6           EI = 500          GJ = 25
                 m = 0.160           am = 0.05        wet = 1.15
                 d = 0.0063          Cdt = 0.01       Cdn = 1.5

   array        EA = 4.4e6           EI = 500         GJ = 25
                 m = 1.484           am = 0.631       wet = 8.514
                 d = 0.028           Cdt = 0.01       Cdn = 2.0

   chain        EA = 5.5e7           EI = 0.01        GJ = 0.1
                 m = 4.4             am = 0.584      wet = 38.1
                 d = 0.0265         Cdt = 0.01       Cdn = 1.0

   pole         EA = 5.5e7           EI = 0.01        GJ = 0.1
                 m = 4.4             am = 0.584      wet = 38.1
                 d = 0.0265         Cdt = 0.01       Cdn = 1.0

   b_tether     EA = 1.0e7           EI = 800          GJ = 40
                 m = 4.160           am = 2.08        wet = -4.45
                 d = 0.05            Cdt = 0.01       Cdn = 1.5

   w_tether     EA = 1.0e7           EI = 800          GJ = 40
                 m = 4.160           am = 2.08        wet = 2.67
                 d = 0.05            Cdt = 0.01       Cdn = 1.5

   snubber      EA = 3.13e4          EI = 1000        GJ = 500
                 m = 14.07           am = 9.74        wet = 42.46
                 d = 0.11            Cdt = 0.01       Cdn = 1.5

Layout
   terminal = {
       anchor = clump
   }
   segment = {
       length = 2900
       material = wire
       nodes = (200, 1.0)
   }
   connector = shackle
   segment = {                  /* this is the pinger cage/spare pole   */
       length = 4
       material = pole
       nodes = (5, 1.0)
   }
   connector = shackle
   segment = {
       length = 100
       material = wire
       nodes = (20, 1.0)
   }
   connector = shackle
   segment = {                  /* this is the docking station plus some  */
       length = 6.5             /* stuff fudged above and below it        */
       material = pole
       nodes = (5, 1.0)
   }
   connector = shackle  
   segment = {                  /* this is the EM cable - what is it?   */
       length = 400
       material = array
       nodes = (50, 1.0)
   }
   connector = shackle
   segment = {                          /* 3 m of potted chain below SSF */
       length = 3
       material = chain
       nodes = (5, 1.0)
   }
   connector = ssf
   segment = {                          /* 3 m of potted chain above SSF */
       length = 3
       material = chain
       nodes = (5, 1.0)
   }
   connector = shackle
   segment = {
       length = 75
       material =  b_tether
       nodes = (75, 1.0)
   }
   connector = shackle
   segment = {
       length = 75
       material = w_tether
       nodes = (75, 1.0)
   }
   connector = shackle
   segment = {
      length = 15
      material = snubber
      nodes = (30, 1.0)
   }
   terminal = {
       buoy = float
   }

End
