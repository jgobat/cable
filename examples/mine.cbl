Problem Description
    title = "Tethered Mine Example Problem"
    type = subsurface

Analysis Parameters
    duration       = 10.0
    time-step      = 0.1
    relaxation     = 0.1
    max-iterations = 300
    tolerance      = 1e-6
    static-solution = relaxation

Environment
    forcing-method = morison
    input-type     = regular
    x-wave         = (1.0, 8.0, 0.0)
    y-current      = 0.2*exp(-0.055*H)
    x-current      = 0.2*exp(-0.055*H)
    rho            = 1027
    gravity        = 9.81
    depth          = 25
    
Buoys
   mine         type = sphere
                d = 1.5           m = 1611
                Cdn = 0.5         Cdt = 0.5
                
Anchors
   clump        

Materials       
   rope         EA = 5.0e4       EI = 20         GJ = 20
                m = 0.09         am = 0.08        wet = 0.087
                d = 0.01         Cdt = 0.01       Cdn = 1.5


Connectors
    foo d = 0.1
        m = 1
        wet = 1
        Cdn = 1.0
        Cdt = 1.0

Layout
   terminal = {
      anchor = clump
   }
   segment = {
       length = 20.0
       material = rope
       nodes = (21, 1.0)
   }
/*   connector = foo 
   segment = {
       length = 10.0
       material = rope
       nodes = (21, 1.0)
   }
*/
   terminal = {
       buoy = mine
   }

End
