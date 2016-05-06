Problem Description
    title = "Seatex WaveScan Buoy Deployment"
    type  = surface

Analysis Parameters
    duration                = 150.0
    time-step               = 0.05

    dynamic-relaxation      = 1.0
    dynamic-iterations      = 50
    dynamic-tolerance       = 1e-8
    ramp-time               = 10.0

    static-tolerance        = 1e-5
    static-relaxation       = 0.001
    static-iterations       = 1000000
   
    static-solution = shooting

    static-outer-iterations = 100
    static-outer-relaxation = 0.95

    static-outer-tolerance  = 0.01

Environment
    rho            = 1025
    gravity        = 9.81
    x-current      = 0.5
    depth          = 41

Buoys
   wavescan     type = axisymmetric
                m = 870
                h = 4.0
                diameters = (0.0, 0.94) (0.39, 0.94)
                            (1.18, 2.7) (1.97, 0.94)
                Cdn = 1.0

Anchors
   clump

Connectors
   panther	d = 0.25
		m = 3.13
		wet = -39.6
		Cdn = 1.0
		
   shackle      d = 0.375*0.0254
                wet = 4.0
                Cdn = 0.0
                m = 0.5

   sphere_48in  d = 1.2307
                wet = -1480*4.448
                Cdn = 0.6
                m = 366

Materials
   chain_3_8    EA = 3.5e+07      EI = 1e-4         GJ = 0.1
                 m = 2.365                         wet = 20.1715
                 d = 0.0364      Cdt = 0.05        Cdn = 0.55

Layout
   terminal = {
      anchor = clump
   }
   segment = {
       length = 11 
       material = chain_3_8
       nodes = (12, 1.0)
   }
   segment = {
       length = 40
       material = chain_3_8
       nodes = (42, 1.0)
       attachments = panther : (2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
				20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 
				35 36 37 38 39 40 41)
   }
   segment = {
       length = 30 
       material = chain_3_8
       nodes = (62, 1.0)
       attachments = panther : (2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19
			        20 21 22 23 24 25 26 27 28 29 30 31 32 33 34
				35 36 37 38 39 40 41 42 43 44 45 46 47 48 49
				50 51 52 53 54 55 56 57 58 59 60 61)
   }
   terminal = {
      buoy = wavescan
   }

End
