#Density - table of variables of LB Node to stream
#	Pressure Evolution:
AddDensity( name="g[0]", dx= 0, dy= 0, group="g")
AddDensity( name="g[1]", dx= 1, dy= 0, group="g")
AddDensity( name="g[2]", dx= 0, dy= 1, group="g")
AddDensity( name="g[3]", dx=-1, dy= 0, group="g")
AddDensity( name="g[4]", dx= 0, dy=-1, group="g")
AddDensity( name="g[5]", dx= 1, dy= 1, group="g")
AddDensity( name="g[6]", dx=-1, dy= 1, group="g")
AddDensity( name="g[7]", dx=-1, dy=-1, group="g")
AddDensity( name="g[8]", dx= 1, dy=-1, group="g")

#	Phase Field Evolution:
AddDensity( name="h[0]", dx= 0, dy= 0, group="h")
AddDensity( name="h[1]", dx= 1, dy= 0, group="h")
AddDensity( name="h[2]", dx= 0, dy= 1, group="h")
AddDensity( name="h[3]", dx=-1, dy= 0, group="h")
AddDensity( name="h[4]", dx= 0, dy=-1, group="h")
AddDensity( name="h[5]", dx= 1, dy= 1, group="h")
AddDensity( name="h[6]", dx=-1, dy= 1, group="h")
AddDensity( name="h[7]", dx=-1, dy=-1, group="h")
AddDensity( name="h[8]", dx= 1, dy=-1, group="h")

if (Options$thermo){
	AddDensity( name="f[0]", dx= 0, dy= 0, group="f")
	AddDensity( name="f[1]", dx= 1, dy= 0, group="f")
	AddDensity( name="f[2]", dx= 0, dy= 1, group="f")
	AddDensity( name="f[3]", dx=-1, dy= 0, group="f")
	AddDensity( name="f[4]", dx= 0, dy=-1, group="f")
	AddDensity( name="f[5]", dx= 1, dy= 1, group="f")
	AddDensity( name="f[6]", dx=-1, dy= 1, group="f")
	AddDensity( name="f[7]", dx=-1, dy=-1, group="f")
	AddDensity( name="f[8]", dx= 1, dy=-1, group="f")
}

if (Options$Outflow) {
	AddDensity( name=paste("gold",0:8,sep=""), dx=0, dy=0, group="gold")
	AddDensity( name=paste("hold",0:8,sep=""), dx=0, dy=0, group="hold")
}

#	Fields required for solid contact
AddDensity(name="nw_x", dx=0, dy=0, group="nw", comment='phase field normal at the wall in x direction, pointing into fluid')
AddDensity(name="nw_y", dx=0, dy=0, group="nw")

#	Velocity Fields
AddDensity(name="U", dx=0, dy=0, group="Vel")
AddDensity(name="V", dx=0, dy=0, group="Vel")

#	Phase-field stencil for finite differences
AddField('PhaseF',stencil2d=1, group="PF")

#	Additional access required for outflow boundaries
if (Options$Outflow){
	for (d in rows(DensityAll)){
    		AddField( name=d$name,  dx=-d$dx-1, dy=-d$dy )
    		AddField( name=d$name,  dx=-d$dx, dy=-d$dy-1 )
	}
	AddField('U',dx=c(-1,0))
	AddField('V',dx=c(0,-1))
}

if (Options$thermo){
	AddDensity("Temp", dx=0, dy=0, group="Thermal")
	AddDensity("Cond", dx=0, dy=0,group="Thermal")
	AddDensity("SurfaceTension", dx=0, dy=0, group="Thermal")

	# Add as fields so we can copy if not running hydrodynamics
	AddField( name="g[0]",group="g")
	AddField( name="g[1]",group="g")
	AddField( name="g[2]",group="g")
	AddField( name="g[3]",group="g")
	AddField( name="g[4]",group="g")
	AddField( name="g[5]",group="g")
	AddField( name="g[6]",group="g")
	AddField( name="g[7]",group="g")
	AddField( name="g[8]",group="g")
	AddField( name="h[0]",group="h")
	AddField( name="h[1]",group="h")
	AddField( name="h[2]",group="h")
	AddField( name="h[3]",group="h")
	AddField( name="h[4]",group="h")
	AddField( name="h[5]",group="h")
	AddField( name="h[6]",group="h")
	AddField( name="h[7]",group="h")
	AddField( name="h[8]",group="h")
	AddField( name="f[0]",group="f")
	AddField( name="f[1]",group="f")
	AddField( name="f[2]",group="f")
	AddField( name="f[3]",group="f")
	AddField( name="f[4]",group="f")
	AddField( name="f[5]",group="f")
	AddField( name="f[6]",group="f")
	AddField( name="f[7]",group="f")
	AddField( name="f[8]",group="f")
	AddField("Temp",stencil2d=1, group="Thermal")
	AddField("Cond",stencil2d=1, group="Thermal")
	AddField("SurfaceTension",stencil2d=1, group="Thermal")
	AddField(name="U",group="Vel")
	AddField(name="V",group="Vel")

	# Temperature outputs:
	AddQuantity(name="temperature",unit="K")
	AddQuantity(name="surfaceTension",unit="N/m")
	AddQuantity(name="thermalconductivity")

	# Temperature inputs:
	AddSetting("surfPower",	default="1", comment="Use for parabolic representation of surface tension")
	AddSetting("sigma_T",		default="-0.00017",	comment="Derivative describing how surface tension changes with temp unit=[N/m2]")
	AddSetting("sigma_TT",			default="1",comment="Derivative describing how surface tension changes with temp unit=[N/m3]")
	AddSetting("T_ref",			default="1",	comment="Reference temperature at which sigma is set unit=[K]")
	AddSetting("T_init",zonal=T, 	comment="Initial temperature field unit=[K]")
	AddSetting("cp_h",			default="2",	comment="specific heat for heavy phase unit=[J/kg/K]")
	AddSetting("cp_l",			default="1",	comment="specific heat for light phase unit=[J/kg/K]")
	AddSetting("k_h", 			default="0.5",	comment="thermal conductivity for heavy phase unit=[W/m/K]")
	AddSetting("k_l", 			default="0.024",	comment="thermal conductivity for light phase unit=[W/m/K]")
	AddSetting("Beta_h",			default="0.000002",	comment="thermal expansion coefficient for heavy phase unit=[]")
	AddSetting("Beta_l",			default="0.004",	comment="thermal expansion coefficient for light phase unit=[]")
	AddSetting("dT",				comment="Application of vertical temp gradient to speed up initialisation unit=[K]")
	AddSetting("dTx", default="0",	comment="Application of horizontal temp gradient to speed up initialisation unit=[K]")	
	AddSetting("stabiliser", default="1",	comment="If not solving flow field, can adjust temperature timestep")	
	AddGlobal("TempChange")
	
	# Add stage for temperature update
	AddStage("CopyDistributions","TempCopy",save=Fields$group %in% c("g","h","f","Thermal","Vel","nw"))
	AddStage("CopyThermal","ThermalCopy", save=Fields$name %in% c("Temp","Cond","SurfaceTension"), load=DensityAll$name %in% c("Temp","Cond","SurfaceTension"))
	AddStage("TempInit", "Init_temp", save=Fields$group %in% c("Thermal"))
	AddStage("TempIter", "Iter_temp", save=Fields$group %in% c("Thermal"), load=DensityAll$group %in% c("f"))
	AddStage("Init_temp_distro", "temp_init_f", save=Fields$group %in% c("f"))
	AddStage("Iter_temp_distro", "temp_iter_f", save=Fields$group %in% c("f"), load=DensityAll$group %in% c("f"))
}	

#	Stages - processes to run for initialisation and each iteration
if (Options$RT) {
 	AddField('PhaseOld', group="PF")

	# initialisation
	AddStage("PhaseInit" , "Init_phase" 		, save=Fields$group %in% c("PF"))
	AddStage("WallInit"  , "Init_wallNorm"    	, save=Fields$group %in% c("nw"))
	AddStage("BaseInit"  , "Init_distributions" , save=Fields$group %in% c("g","h","Vel"))
	
	# iteration
	AddStage("BaseIter"  , "calcHydroIter" 		, save=Fields$group %in% c("g","h","Vel")		, load=DensityAll$group %in% c("g","h","Vel"))
    AddStage("PhaseIter" , "calcPhaseFIter"		, save=Fields$name  %in% c("PhaseF","PhaseOld") , load=DensityAll$group=="h")
  	AddStage("WallIter"  , "calcWallPhaseIter"  , save=Fields$group %in% c("PF")				, load=DensityAll$group=="nw") 
#} else if (Options$Thermo) {
#	AddStage("", save=Fields$group %in% c("f"))	#Learn Options, Stages & Actions so this can be properly added
} else if (Options$Outflow) {
	# initialisation
	AddStage("PhaseInit" , "Init_phase"			, save=Fields$group %in% c("PF"))
	AddStage("WallInit"  , "Init_wallNorm"		, save=Fields$group %in% c("nw"))
	AddStage("BaseInit"  , "Init_distributions"	, save=Fields$group %in% c("g","h","Vel","gold","hold")) 
	# iteration
	AddStage("BaseIter"  , "calcHydroIter"      , save=Fields$group %in% c("g","h","Vel","nw","gold","hold"), 
												  load=DensityAll$group %in% c("g","h","Vel","nw","gold","hold")) 
	AddStage("PhaseIter" , "calcPhaseFIter"		, save=Fields$group %in% c("PF"), load=DensityAll$group %in% c("g","h","Vel","nw","gold","hold"))
	AddStage("WallIter"  , "calcWallPhaseIter"	, save=Fields$group %in% c("PF"), load=DensityAll$group=="nw")	
} else if (Options$thermo) {
	# initialisation
	AddStage("PhaseInit" , "Init_phase"			, save=Fields$group %in% c("PF"))
	AddStage("WallInit"  , "Init_wallNorm"		, save=Fields$group %in% c("nw"))
	AddStage("BaseInit"  , "Init_distributions" , save=Fields$group %in% c("g","h","f","Vel"))

	# iteration
	AddStage("BaseIter"  , "calcHydroIter"      , save=Fields$group %in% c("g","h","f","Vel","nw") , load=DensityAll$group %in% c("g","h","f","Vel","nw"))  # TODO: is nw needed here?
	AddStage("PhaseIter" , "calcPhaseFIter"		, save=Fields$group %in% c("PF")			   , load=DensityAll$group %in% c("g","h","f","Vel","nw"))
	AddStage("WallIter"  , "calcWallPhaseIter"	, save=Fields$group %in% c("PF")			   , load=DensityAll$group %in% c("nw"))	
} else {
	# initialisation
	AddStage("PhaseInit" , "Init_phase"			, save=Fields$group %in% c("PF"))
	AddStage("WallInit"  , "Init_wallNorm"		, save=Fields$group %in% c("nw"))
	AddStage("BaseInit"  , "Init_distributions" , save=Fields$group %in% c("g","h","Vel"))

	# iteration
	AddStage("BaseIter"  , "calcHydroIter"      , save=Fields$group %in% c("g","h","Vel","nw") , load=DensityAll$group %in% c("g","h","Vel","nw"))  # TODO: is nw needed here?
	AddStage("PhaseIter" , "calcPhaseFIter"		, save=Fields$group %in% c("PF")			   , load=DensityAll$group %in% c("g","h","Vel","nw"))
	AddStage("WallIter"  , "calcWallPhaseIter"	, save=Fields$group %in% c("PF")			   , load=DensityAll$group %in% c("nw"))	
}

if (Options$thermo){
	AddAction("Iteration", c("BaseIter", "PhaseIter","WallIter","TempIter"))
	AddAction("Init"     , c("PhaseInit","WallInit", "WallIter","TempInit","BaseInit"))
} else {
	AddAction("Iteration", c("BaseIter", "PhaseIter","WallIter"))
	AddAction("Init"     , c("PhaseInit","WallInit", "WallIter","BaseInit"))
}
# 	Outputs:
AddQuantity(name="Rho",	  unit="kg/m3")
AddQuantity(name="PhaseField",unit="1")
AddQuantity(name="U",	  unit="m/s",vector=T)
AddQuantity(name="NormalizedPressure",	  unit="Pa")
AddQuantity(name="Pressure",	  unit="Pa")
AddQuantity(name="Normal", unit="1", vector=T)

#	Initialisation States
AddSetting(name="Period", default="0", comment='Number of cells per cos wave')
AddSetting(name="Perturbation", default="0", comment='Size of wave perturbation, Perturbation Period')
AddSetting(name="MidPoint", default="0", comment='height of RTI centerline')

AddSetting(name="Radius" , default="0", comment='Radius of diffuse interface circle')
AddSetting(name="CenterX", default="0", comment='Circle center x-coord')
AddSetting(name="CenterY", default="0", comment='Circle Center y-coord')
AddSetting(name="BubbleType", default="1", comment='Drop/bubble')

#	Inputs: For phasefield evolution
AddSetting(name="IntWidth", default=4,    comment='Anti-diffusivity coeff')#Added for F_surftens
AddSetting(name="Density_h", comment='High density fluid')
AddSetting(name="Density_l", comment='Low  density fluid')
AddSetting(name="PhaseField_h", default=1, comment='PhaseField in high density fluid')
AddSetting(name="PhaseField_l", default=0, comment='PhaseField in low density fluid')
AddSetting(name="PhaseField_init", 	   comment='Initial/Inflow PhaseField distribution', zonal=T)
AddSetting(name="W", default=4,    comment='Anti-diffusivity coeff (phase interfacial thickness) ')
AddSetting(name="omega_phi", comment='one over relaxation time (phase field)')
AddSetting(name="M", omega_phi='1.0/(3*M+0.5)', default=0.02, comment='Mobility')
AddSetting(name="sigma", 		default=0.076 ,  comment='surface tension')
AddSetting(name="radAngle",default='1.570796', comment='Contact angle in radians, can use units -> 90d where d=2pi/360', zonal=T)

# 	Inputs: Fluid Properties
AddSetting(name="tau_l", comment='relaxation time (low density fluid)')
AddSetting(name="tau_h", comment='relaxation time (high density fluid)')
AddSetting(name="dynamic_viscosity_l", default=0.1,comment='dynamic viscosity light')
AddSetting(name="dynamic_viscosity_h", default=0.1,comment='dynamic viscosity heavy')
AddSetting(name="Viscosity_l", tau_l='(3*Viscosity_l)', default=0.16666666, comment='kinematic viscosity')
AddSetting(name="Viscosity_h", tau_h='(3*Viscosity_h)', default=0.16666666, comment='kinematic viscosity')

AddSetting(name="omega_bulk", comment='inverse of bulk relaxation time', default=1.0)
AddSetting(name="bulk_visc", omega_bulk='1.0/(3*bulk_visc+0.5)',  comment='bulk viscosity')
#	Inputs: Flow Properties
AddSetting(name="VelocityX", default=0.0, comment='inlet/outlet/init velocity', zonal=T)
AddSetting(name="VelocityY", default=0.0, comment='inlet/outlet/init velocity', zonal=T)
AddSetting(name="Pressure" , default=0.0, comment='inlet/outlet/init density', zonal=T)
AddSetting(name="GravitationX", default=0.0, comment='applied (rho)*GravitationX', zonal=T)
AddSetting(name="GravitationY", default=0.0, comment='applied (rho)*GravitationY', zonal=T)
AddSetting(name="BuoyancyX", default=0.0, comment='applied (rho-rho_h)*BuoyancyX')
AddSetting(name="BuoyancyY", default=0.0, comment='applied (rho-rho_h)*BuoyancyY')
AddSetting(name="fixedIterator", default=2.0, comment='fixed iterator for velocity calculation')
#	Globals - table of global integrals that can be monitored and optimized
AddGlobal(name="PressureLoss", comment='pressure loss', unit="1mPa")
AddGlobal(name="OutletFlux", comment='pressure loss', unit="1m2/s")
AddGlobal(name="InletFlux", comment='pressure loss', unit="1m2/s")
AddGlobal(name="TotalDensity", comment='Mass conservation check', unit="1kg/m3")
AddNodeType(name="SpikeTrack", group="ADDITIONALS")
AddNodeType(name="BubbleTrack", group="ADDITIONALS")
AddGlobal(name="RTIBubble", comment='Bubble Tracker', op="MAX")
AddGlobal(name="RTISpike",  comment='Spike Tracker', op="MAX")
AddGlobal(name="NMovingWallForce", comment='force exerted on the N Moving Wall')
AddGlobal(name="NMovingWallPower", comment='implented: Vx* incoming momentum (precollision)')

if (Options$debug){
	AddGlobal(name="MomentumX", comment='Total momentum in the domain', unit="")
	AddGlobal(name="MomentumY", comment='Total momentum in the domain', unit="")
	AddGlobal(name="MomentumX_afterCol", comment='Total momentum in the domain', unit="")
	AddGlobal(name="MomentumY_afterCol", comment='Total momentum in the domain', unit="")

	AddGlobal(name="F_pressureX", comment='Pressure force X', unit="")
	AddGlobal(name="F_pressureY", comment='Pressure force Y', unit="")
	AddGlobal(name="F_bodyX", comment='Body force X', unit="")
	AddGlobal(name="F_bodyY", comment='Body force Y', unit="")
	AddGlobal(name="F_surf_tensionX", comment='Surface tension force X', unit="")
	AddGlobal(name="F_surf_tensionY", comment='Surface tension force Y', unit="")
	AddGlobal(name="F_muX", comment='Viscous tension force X', unit="")
	AddGlobal(name="F_muY", comment='Viscous tension force Y', unit="")
	AddGlobal(name="F_total_hydroX", comment='Total hydrodynamic force X', unit="")
	AddGlobal(name="F_total_hydroY", comment='Total hydrodynamic force Y', unit="")

	AddGlobal(name="F_phiX", comment='Forcing term for interface tracking X', unit="")
	AddGlobal(name="F_phiY", comment='Forcing term for interface tracking Y', unit="")
}
#	Node things
if (Options$CM){
	AddNodeType(name="CM", group="COLLISION")  # Central Moments collision
}
AddNodeType(name="Smoothing", group="ADDITIONALS")  #  To smooth phase field interface during initialization.

#	Boundary things
AddNodeType(name="MovingWall_N", group="BOUNDARY")
AddNodeType(name="MovingWall_S", group="BOUNDARY")
AddNodeType(name="NVelocity", group="BOUNDARY")
AddNodeType(name="WVelocity", group="BOUNDARY")

AddNodeType(name="Body", group="BODY")  # To measure force exerted on the body.

AddGlobal(name="FDrag", comment='Force exerted on body in X-direction', unit="N")
AddGlobal(name="FLift", comment='Force exerted on body in Y-direction', unit="N")
AddGlobal(name="FTotal", comment='Force exerted on body in X+Y -direction', unit="N")

if (Options$Outflow) {
	AddNodeType(name="Convective_E", group="BOUNDARY")
	AddNodeType(name="Convective_N", group="BOUNDARY")
	AddNodeType(name="Neumann_E", group="BOUNDARY")
}
AddNodeType(name="Solid", group="BOUNDARY")
AddNodeType(name="Wall", group="BOUNDARY")
AddNodeType(name="BGK", group="COLLISION")
AddNodeType(name="MRT", group="COLLISION")
