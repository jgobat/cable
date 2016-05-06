Problem Description
    title = "Stratus 2 Surface Mooring"
    type  = surface

Analysis Parameters
    static-relaxation = 0.1
    static-iterations = 5000
    static-tolerance  = 0.001

    relax-adapt-up    = 1.02
    relax-adapt-down  = 1.1

    static-outer-iterations = 1000

    static-outer-relaxation = 0.98
    static-outer-tolerance  = 0.01
  

    duration           = 200.0
    ramp-time          = 10.0
    time-step          = 0.1    
    dynamic-relaxation = 1.0
    dynamic-iterations = 15
    dynamic-tolerance  = 1e-6
    dynamic-rho        = -0.5

Environment
    rho              = 1027     
    gravity          = 9.81
    depth            = 4440
    bottom-stiffness = 0      
    bottom-damping   = 0
    x-current        = (0.0, 0.1)(100,.05)(500,.01)(4440.00, 0.0)
    x-wind           = 8
    forcing-method   = velocity
    z-input	     = (1.5,9,0)
    input-type 	     = random


Buoys
    discus_3m	type = axisymmetric	
		m = 1145
		diameters = (0.0, 0.02) (1.99, 0.02) (2.0, 0.91) 
 		  	    (2.3, 2.6) (2.8, 3.0) (2.96, 3.0)
		Cdn = 1.0   Cdt = 1.0
		Sw = 6.0    Cdw = 1.5

Anchors
    clump

Connectors
    pivot	d = 0.01	/* to release moments - use segments mostly */
		m = 0.01
		wet = 0.01
		Cdn = 0.0
		length = 0.0

Materials

chain_3/4in	 EA = 1.3e+08      EI = 1e-4         GJ = 1e-4
                  m = 8.643                         wet = 73.7084
                  d = 0.0681       Cdt = 0.05       Cdn = 0.55
				   Cat = 0.1	    Can = 1.0
		  comment = "3/4 Crosby Proof Coil"

nilspin_7/16in  EA = 6.7304e+06 EI = 1.0        GJ = 1.0
                 m = 0.58316                    wet = 4.1009
                 d = 0.5625     Cdn = 1.5       Cdt = 0.005
                                Can = 1.0       Cat = 0.0
                 comment = "Nilspin jacketed 7/16in wire rope (WRCA)"
nilspin_3/8in   EA = 4.944e+06  EI = 1.0        GJ = 1.0
                 m = 0.446                      wet = 3.109
                 d = 0.5        Cdn = 1.5       Cdt = 0.005
                                Can = 1.0       Cat = 0.0
                 comment = "Nilspin jacketed 3/8in wire rope (WRCA)"
nylon_7/8in     EA = 84071      EI = 0.01       GJ = 0.01
                 m = 0.29                       wet = 0.2846
                 d = 0.02222    Cdn = 1.5       Cdt = 0.005
                                Can = 1.0       Cat = 0.0
                 comment = "7/8 in nylon - Columbian Rope Corportation"
nylon_1in       EA = 100084     EI = 0.01       GJ = 0.01
                 m = 0.3719                     wet = 0.3648
                 d = 0.0254     Cdn = 1.5       Cdt = 0.005
                                Can = 1.0       Cat = 0.0 
                 comment = "1 in nylon - Columbian Rope Corportation"
polypro_1_1/8in EA = 77265      EI = 0.01       GJ = 0.01
                 m = 0.3942                     wet = -0.4931
                 d = 0.02858    Cdn = 1.5       Cdt = 0.005
                                Can = 1.0       Cat = 0.0
                 comment = "1-1/8 in polypro - Columbian Rope Corporation"
balls_1/2in     EA = 6.0e+07       EI = 1e-4         GJ = 1e-4
                 m = 21.712                         wet = -217.5
                 d = 0.0475        Cdn = 2.6        Cdt = 0.7
                                   amn = 25.0       amt = 25.0
                 comment = "17in glass balls on 1/2in trawler, 1 per meter"
trawler_1/2in   EA = 6.0e+07       EI = 1e-4         GJ = 1e-4
                 m = 3.712                          wet = 31.64
                 d = 0.0475        Cdt = 0.05       Cdn = 0.55
                                   Cat = 0.1        Can = 1.0
                 comment = "1/2 in Trawlex - Acco"
nystron_1in     EA = 748784     EI = 0.01       GJ = 0.01
                 m = 0.4700                     wet = 0.7957
                 d = 0.0254     Cdn = 1.5       Cdt = 0.005
                                Can = 1.0       Cat = 0.0
                 comment = "1 in Nystron - The American Group"

release_322     EA = 5.0e7         EI = 10000        GJ = 10000
                 m = 30.6                           wet = 129.8
                 d = 0.15         Cdt = 0.05        Cdn = 1.0
            length = 1.8          Cat = 0.01        Can = 1.0
                 comment = "EG&G 322 acoustic release - made up"




Layout

terminal = {
       anchor = clump
      }
segment = {
	length = 5
	material = trawler_1/2in
	nodes = (5,1.0)
}
connector = pivot

segment = {
	length = 20
	material = nystron_1in
	nodes = (5,1.0)
}
connector = pivot

segment = {
	length = 5
	material = trawler_1/2in
	nodes = (5,1.0)
}
connector = pivot

segment = {
	length = 1.8
	material = release_322
	nodes = (2,1.0)
}
connector = pivot

segment = {
	length = 5
	material = trawler_1/2in
	nodes = (5,1.0)
}
connector = pivot

segment = {
	length = 80
	material = balls_1/2in
	nodes = (10,1.0)
}
connector = pivot

segment = {
	length = 1400
	material = polypro_1_1/8in
	nodes = (100,1.0)
}
connector = pivot

segment = {
	length = 100
	material = nylon_1in
	nodes = (10,1.0)
}
segment = {
	length = 1850
	material = nylon_7/8in
	nodes = (100,1.0)
}
connector = pivot
	
segment = {
	length = 1750
	material = nilspin_3/8in
	nodes = (200,1.0)
}
connector = pivot

segment = {
	length = 210
	material = nilspin_7/16in
	nodes = (20,1.0)
}
connector = pivot
segment = {
 	length = 40
	material = chain_3/4in
	nodes = (5,1.0)
}
connector = 	pivot

terminal = {
	buoy = discus_3m
}

end



