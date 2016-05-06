Problem Description
	title = "Example of a Taut Surface Mooring"
	type  = surface

Analysis Parameters
	duration          = 100
	time-step         = 0.01
	relaxation        = 1.0
	max-iterations    = 100
	tolerance         = 0.01
    dynamic-tolerance = 1e-6
	ramp-time	  = 10
	static-solution   = catenary

Environment
	forcing-method = wave-follower
	input-type     = regular
	x-wave         = (1.5, 5.0, 0.0)
	x-current      = 0.25
	rho            = 1026
	gravity        = 9.81		
	depth          = 73
	bottom-stiffness = 100.0
	bottom-damping   = 1.0 
    
Buoys
	aass	type = cylinder
		d = 0.165
		m = 27.2
		h = 1.83
		Cdn = 1.0

Anchors
	clump        

Connectors
   shackle	d   = 0.03 /* 0.03 m = 1.25 in */
		m   = 0.05 /* 1/4" shackle = 1.9 ounces (.12 lb) */
		wet = 0.44 /* approx. 90% of dry weight */
		Cdn = 1.0 
		
Materials       
	rope	EA = 8734.69 /* E = 4.59683e8 Pa, A = 60%(pi)R^2 = 1.9e-5 m^2 */
		EI = 0.018   /* Irod = (pi*R^4)/4, EIrope = 50% EIrod         */
		GJ = 0.01    /* Jrod = 1/2pi*R^4 G = ? Not used in 2D         */
		d = 0.00635  /* 0.00635 m = 0.25 in                           */
		m = 0.022    /* 0.022 kg/m = 0.15 lbf/ft                      */
		wet = 0.0028 /* 12.6% weight in air                           */
		Cdn = 1.4
		Cdt = 0.01

	wire	EA = 228666  /* E = 1.9305e11 Pa, A = .47*D^2 = 1.1845e-6 m^2 */
		EI = 0.045   /* Irod = (pi*R^4)/4, EIcable = 75% EIrod        */
		GJ = 0.03    /* Jrod = 1/2pi*R^4 G = ? Not used in 2D         */
		d = 0.001588 /* 0.001588 m = 0.0625 in = 0.005208 ft          */
		m = 0.011161 /* 0.011161 kg/m = 0.0075 lbf/ft                 */
		wet = 0.0975 /* = mg - (.47D^2)*rho*g                         */
		Cdn = 1.4
		Cdt = 0.01
	
Layout
	terminal = {
           anchor = clump
        }
	segment = {
	   length = 5.0	
           material = rope
	   nodes = (81, 0.8) (40, 0.2)   
	}
	connector = shackle
	segment = {
	   length = 70.0
	   material = wire
	   nodes = (141, 0.05) (266, 0.95) 
	}
	terminal = {
	   buoy = aass
        }

End
