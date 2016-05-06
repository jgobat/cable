Problem Description
    title = "CMO Central Discus Mooring: Baseline case"
    type  = surface

Analysis Parameters
    static-relaxation = 0.1
    static-iterations = 20000
    static-tolerance  = 0.001

    relax-adapt-up    = 1.02
    relax-adapt-down  = 1.1

    static-outer-iterations = 100

    static-outer-relaxation = 0.98
    static-outer-tolerance  = 0.01
    static-initial-guess = shooting
    static-solution   = catenary
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
    depth            = 70.0
    bottom-stiffness = 100      
    bottom-damping   = 20
    x-current        = 0.5
    x-wind           = 17.5
    forcing-method   = wave-follower
    x-wave           = (0.0, 10.0, 0.0)

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
    chain_3_4   EA = 5.0e7	EI = 10.0		GJ = 10.0
		 m = 8.0                      		wet = 68.0
		 d = 0.066	Cdt = 0.006		Cdn = 0.5
				Cmt = 1.1		Cmn = 2.0

    A_connect   EA = 5.0e7	EI = 10.0		GJ = 10.0
		 m = 15.5 + 9.2/0.3 /* MTR */		wet = 132
		 d = 0.1	Cdt = 0.01		Cdn = 0.5
				Cmt = 1.1		Cmn = 2.0

    B_connect   EA = 5.0e7	EI = 10.0		GJ = 10.0
		 m = 12.18				wet = 104.07
		 d = 0.083	Cdt = 0.005		Cdn = 0.5
				Cmt = 1.2		Cmn = 2.0

    seacat	EA = 5.0e7	EI = 1000.0		GJ = 1000.0
		 m = 22.5				wet = 160
		 d = 0.1	Cdt = 0.01 		Cdn = 1.0
				Cmt = 1.09		Cmn = 2.0

    VMCM	EA = 5.0e7	EI = 1000.0		GJ = 1000.0
		 m = 39.1 				wet = 203.86
		 d = 0.30 	Cdt = 0.01 		Cdn = 1.0
				Cmt = 1.2  		Cmn = 1.5

    depressor   EA = 5.0e7	EI = 10000.0		GJ = 10000.0
		 m = 573.16				wet = 4886.2
		 d = 0.305	Cdt = 0.01		Cdn = 1.0
				Cmt = 1.26		Cmn = 2.0

    swivel_5T   EA = 5.0e7	EI = 1000.0		GJ = 1000.0
		 m = 26.48				wet = 225.7
		 d = 0.057	Cdt = 0.01 		Cdn = 1.0
		 		Cmt = 1.18		Cmn = 2.0

    release	EA = 5.0e7	EI = 10000.0		GJ = 10000.0
		 m = 30.6				wet = 129.8
		 d = 0.1	Cdt = 0.01 		Cdn = 1.0
				Cmt = 1.05		Cmn = 2.0

Layout
   terminal = {
      anchor = clump
      release-time = 5
   }
   segment = {
      length = 70.0
      material = chain_3_4
      nodes = (71, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = { 
      length = 0.21
      material = swivel_5T
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 15.07
      material = chain_3_4
      nodes = (16, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = { 
      length = 1.44
      material = release
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 5.0
      material = chain_3_4
      nodes = (6, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = { 
      length = 0.76
      material = depressor
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 5.0
      material = chain_3_4
      nodes = (6, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = { 
      length = 0.76
      material = seacat
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 2.83
      material = chain_3_4
      nodes = (4, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = { 
      length = 2.07
      material = VMCM
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 3.26
      material = chain_3_4
      nodes = (4, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = { 
      length = 0.76
      material = seacat
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 2.83
      material = chain_3_4
      nodes = (4, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = { 
      length = 2.07
      material = VMCM
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 2.26
      material = chain_3_4
      nodes = (4, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = { 
      length = 2.07
      material = VMCM
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 1.03
      material = chain_3_4
      nodes = (3, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = { 
      length = 0.76
      material = seacat
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.33
      material = chain_3_4
      nodes = (3, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = { 
      length = 2.07
      material = VMCM
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.83
      material = chain_3_4
      nodes = (3, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = { 
      length = 0.76
      material = seacat
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.78
      material = chain_3_4
      nodes = (3, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = { 
      length = 2.07
      material = VMCM
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = { 
      length = 0.21
      material = swivel_5T
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.27
      material = B_connect
      nodes = (2, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.35
      material = chain_3_4
      nodes = (3, 1.0)
   }
   connector = pivot
   segment = {
      length = 0.3
      material = A_connect
      nodes = (2, 1.0)
   }
   terminal = {
      buoy = discus_3m
   }

End
