Problem Description
    title = "Horizontal Array - August 1998 Deployment"
    type = horizontal


Analysis Parameters
    static-relaxation       = 0.3
    static-iterations       = 100000
    static-tolerance        = 1e-4

    relax-adapt-up = 1.02
    relax-adapt-down = 1.1

    static-outer-iterations = 500
    static-outer-tolerance  = 0.001
    static-outer-relaxation = 10

    duration                = 200.0
    time-step               = 0.5
    dynamic-iterations      = 50
    dynamic-relaxation	    = 1.0
    dynamic-tolerance       = 1e-8


Environment
    forcing-method  = morison
    input-type      = random
    x-current 	    = 0.3
    x-wave 	    = (1.0, 10.5, 0.0)
    rho		    = 1025
    gravity	    = 9.81
    depth           = 84.5

Anchors
   clump	


Buoys
   float        type = cylinder
                m = 23
                h = 8*0.0254
                d = 38.5*0.0254
                Cdn = 1.0

   tp_terminator  type = cylinder 
		  d = 2.35*0.0254
   		  h = 0.5
  		  m = 2.68
		  Cdn = 1.0
		  Cdt = 1.0


   microcat_terminator  type = cylinder
			d = 5.3*0.0254
			h = 17.65*0.0254
			m = 9.1/2.205
			buoyancy = 3.0*4.448
			Cdn = 1.0
			Cdt = 1.0

    sinker		type = cylinder	
			d = 2*.0254
			h = 15*.0254
			m = 75/2.205
			Cdn = 1.0
			Cdt = 1.0

    sinker_2		type = cylinder
			d = 2*.0254
			h = 15*.0254
			m = 50/2.205
			Cdn = 1.0
			Cdt = 1.0

Connectors
   glass        d = 17*0.0254           /* 17 in glass spec from
Benthos */
                m = 18
                wet = -25.4*9.81
                Cdn = 0.6
                Cdt = 0.6


   sphere17_48in  d = 1.2307
                  wet = -1480*4.448 + 2*4.45*4.448*.8667
                  Cdn = 0.6
                  m = 366 + 2*4.45/2.205


   sphere8_48in   d = 1.2307
                  wet = -1170*4.448 + 2*4.45*4.448*.8667 + 20*4.448
                  Cdn = 0.6
                  m = 518 + 2*4.45/2.205 + 10.0


   tp_pod	  m = 2.68
		  wet = 1.66*4.448
		  Cdn = 1.0
		  d = 2.35*0.0254


   termination_B  m = 4.45/2.205
		  wet = 4.45*4.448*.8667
		  Cdn = 0.5
		  d = 0.168


   acm_3d	  m = 20.65 + 2.31*2 + 6.88 + 3.18*2
		  wet = 94.34 - 39.60*2
		  Cdn = 1.0
		  d = 0.585 + .041*2


   acm_3d_glass	  m = 20.65 + 2.31*2 + 6.88 + 3.18*2 + 18
		  wet = 94.34 - 2*250 
		  Cdn = 1.0
		  d = 0.585 + .041*2


   acm_3d_x3	  m = 20.65 + 2.31*2 + 6.88 + 3.18*2
		  wet = 94.34 - 39.60*3
		  Cdn = 1.0
		  d = 0.585 + .041*2


   acm_3d_x6	  m = 20.65 + 2.31*2 + 6.88 + 3.18*2
		  wet = 94.34 - 39.60*6
		  Cdn = 1.0
		  d = 0.585 + .041*2


   microcat       m = 9.04/2.2
		  wet = 6.1*4.448
  		  d = 5.4*0.0254
		  Cdn = 1.0


   T_string_junc  m = 3.89/2.2
		  wet =1.14*4.448 - 39.6*5.0
		  d = 6.5*0.0254
		  Cdn = 1.0


   T_string_junc_glass  m = 3.89/2.205 + 2*17.68
		        wet = 1.14*4.448 - 2*56*4.448
			d = 17*.0254
			Cdn = 1.0	


   SBE_39	  m = 1.672
		  wet = .901*9.81
		  d = 1.9*0.0254
		  Cdn = 1.0

 			 

   T_P_recorder   m = 2.67
		  wet = 1.66*4.448
		  d = 2.35*.0254
		  Cdn = 1.0


Materials
   wire_10mm        EA = 4.7e6           EI = 20          GJ = 5
                     m = 0.320           am = 0.08       wet = 2.64
                     d = 0.01           Cdt = 0.01       Cdn = 1.5


   chain	    EA = 5.5e7           EI = 0.01        GJ = 0.1
                     m = 4.4             am = 0.584      wet = 38.1
                     d = 0.0265         Cdt = 0.01       Cdn = 0.55


   acoustic_rel_8202 EA = 1.0e8           EI = 1000        GJ = 100
                     m = 45/1.13         am = 18/1.13	  wet = 314/1.13
		     d = 0.130		 Cdt = 0.18	  Cdn = 1.2	


   steamer	    EA = 1.0e8		 EI = 0.01	   GJ = 0.1
		     m = 24.9     	 am = 1.0	  wet = 192.94
		     d = 0.05		Cdt = 0.01        Cdn = 0.55


    wire_3_16       EA = 4.4e6           EI = 5          GJ = 1 
                     m = 0.160           am = 0.05        wet = 1.15
                     d = 0.0063          Cdt = 0.01       Cdn = 1.5


    wire_string       EA = 4.4e6           EI = 500          GJ = 25
                     m = 0.160 + 0.067     am = 0.05        wet = 1.15 + 0.23
                     d = 0.0127          Cdt = 0.01       Cdn = 1.0


    spectra_string  EA = 1.0e7       EI = 10.0	       GJ = 1.0
                     m = 0.20 + 0.067                  wet = -0.20*9.81*0.06 + 0.23
                     d = 0.0127      Cdt = 0.01        Cdn = 1.0


 

Layout

   terminal = {

      anchor = clump

   }

   segment = {

       length = 17 - 4*0.3048

       material = steamer

       nodes = (17, 1.0)

   }

   connector = termination_B

   segment = {

       length = 15

       material = chain

       nodes = (15, 1.0)

   }

   connector = termination_B

   segment = {

       length = 63.5

       material = wire_10mm 

       nodes = (53, 1.0)

   }

   connector = termination_B

   segment = {

       length = 2

       material = chain

       nodes = (10, 1.0)

   }

   connector = sphere8_48in		
/*
   branch = {
      segment = {
         length = 30.0
         nodes = (31, 1.0)
         material = wire_3_16
      }
      terminal = {
         buoy = float
      }
   }
*/
   segment = {

       length = 20.

       material = wire_10mm

       nodes = (21, 1.0)

   } 

   connector = T_string_junc_glass

   branch = {

      segment = {

         length = 26

         nodes = (27, 1.0)

	 attachments = SBE_39 : (1 6 11 16 21)

	 attachments = T_P_recorder : (26)

         material = spectra_string

      }

      terminal = {

         buoy = sinker

      }

   }

   segment = {

      length = 5

      material = wire_10mm

      nodes = (6, 1.0)

   }

   connector = acm_3d

   segment = {

       length = 25.

       material = wire_10mm

       nodes = (22, 1.0)

   }

   connector = T_string_junc_glass

   branch = {

      segment = {

         length = 26

         nodes = (27, 1.0)

	 attachments = SBE_39 : (1 6 11 16 21)

	 attachments = T_P_recorder : (26)

         material = spectra_string

      }

      terminal = {

         buoy = sinker

      }

   } 

   segment = {

       length = 30.

       material = wire_10mm

       nodes = (31, 1.0)

   }

   connector = acm_3d_glass

   branch = {

      segment = {

         length = 26

         material = wire_3_16

         attachments = microcat : (6 11 16 21 26)

         nodes = (27, 1.0)

      }

      terminal = {

         buoy = sinker_2 

      }

   }

   segment = {

       length = 30.

       material = wire_10mm

       nodes = (31, 1.0)

   }

   connector = T_string_junc_glass

   branch = {

      segment = {

         length = 26

         nodes = (27, 1.0)

         attachments = SBE_39 : (1 6 11 16 21)

	 attachments = T_P_recorder : (26)

         material = spectra_string

      }

      terminal = {

         buoy = sinker

      }

   }

   segment = {

       length = 25.

       material = wire_10mm

       nodes = (26, 1.0)

   }

   connector = acm_3d

   segment = {

       length = 5.

       material = wire_10mm

       nodes = (6, 1.0)

   }

   connector = T_string_junc_glass

   branch = {

      segment = {

         length = 26

         nodes = (27, 1.0)

	 attachments = SBE_39 : (1 6 11 16 21)

	 attachments = T_P_recorder : (26)

         material = spectra_string

      }

      terminal = {

         buoy = sinker

      }

   }

   segment = {

      length = 20

      material = wire_10mm

      nodes = (21, 1.0)

   }

   connector = sphere17_48in
/*
   branch = {
      segment = {
         length = 30.0
         nodes = (31, 1.0)
         material = wire_3_16
      }
      terminal = {
         buoy = float
      }
   }
*/
   segment = {

       length = 2

       material = chain

       nodes = (10, 1.0)

   }

   connector = termination_B

   segment = {

       length = 46.5

       material = wire_10mm

       nodes = (50, 1.0)

   }

   connector = termination_B

   segment = {

        length = 15

        material = chain

        nodes = (5, 1.0)

   } 

   connector = termination_B

   segment = {

       length = 1.2

       material = acoustic_rel_8202

       nodes = (5, 1.0)

   }

   connector = termination_B

   segment = {

        length = 5.0

        material = chain

        nodes = (15, 1.0)

   }

   connector = termination_B

   segment = {

       length = 17 

       material = steamer

       nodes = (20, 1.0)

   }

   terminal = {

      anchor = clump

      z = 0.0           x = 290

   }


End
