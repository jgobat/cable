Problem Description
    title = "Shallow Water Surface Mooring"   
    type  = surface

Analysis Parameters
    static-iterations       = 2000
    static-tolerance        = 1e-4
    static-relaxation       = 0.01
    static-outer-iterations = 100
    static-outer-relaxation = 0.8
    static-outer-tolerance  = 0.01
    static-initial-guess = shooting
    duration                = 120.0
    time-step               = 0.01
    dynamic-relaxation      = 1.0
    dynamic-iterations      = 50
    dynamic-tolerance       = 1e-6
    ramp-time               = 16.0          /* minimize transients         */

Environment
    rho            = 1025
    gravity        = 9.81
    x-current      = 0.5*0.514          
    depth          = 100
    input-type     = random                /* use random input spectrum     */
    forcing-method = wave-follower         
    x-wave         = (48*0.0254, 8.0, 0.0) /* 8 ft sig height, 8 sec period */
    bottom-stiffness = 100.0
    bottom-damping   = 1.0
    bottom-friction  = 0.2
    
Buoys
   float        type = cylinder
                m = 23
                h = 8*0.0254
                d = 36*0.0254
                Cdn = 1.0

Anchors
   clump

Connectors
   shackle      d = 0.375*0.0254
                wet = 4.0
                Cdn = 0.0
                m = 0.5 

Materials       
   chain            EA = 5.5e7           EI = 0.01        GJ = 0.1
                     m = 4.4             am = 0.584      wet = 38.1
                     d = 0.0265         Cdt = 0.01       Cdn = 1.0
   wire             EA = 4.7e6           EI = 20          GJ = 5
                     m = 0.320           am = 0.08       wet = 2.64
                     d = 0.01           Cdt = 0.01       Cdn = 1.5




Layout
   terminal = {
      anchor = clump
   }
   segment = {
       length = 72*0.0254
       material = chain
       nodes = (25, 1.0)
   }
   connector = shackle
   segment = {
       length = 150
       material = wire
       nodes = (150, 1.0)
   }
   connector = shackle
   segment = {
       length = 72*0.0254
       material = chain
       nodes = (25, 1.0)
   }
   terminal = {
      buoy = float
   }

End
