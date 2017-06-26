# Density - table of variables of LB Node to stream
#  name - variable name to stream
#  dx,dy,dz - direction of streaming
#  comment - additional comment
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
AddDensity(name="U", dx=0, dy=0, group="g")
AddDensity(name="V", dx=0, dy=0, group="g")
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

AddField("Solidflag",stencil2d=1,group="nw")
AddField("PhaseF",   stencil2d=1,group="nw")

# Stages - processes to run for initialisation and each iteration
AddStage("PhaseInit"    , 
          save=Fields$group=="nw")
AddStage("BaseInit"     , "Init", 
	  save=Fields$group=="g" | Fields$group=="h" )

AddStage("calcPhaseF",
	  save=Fields$name=="PhaseF",
	  load=DensityAll$group=="h")
AddStage("calcWallPhaseF", 
	  save=Fields$name=="PhaseF")

AddStage("BaseIter"     , "Run" , 
	  save=Fields$group=="g" | Fields$group=="h" , 
	  load=DensityAll$group=="g" | DensityAll$group=="h")

AddAction("Iteration", c("BaseIter", "calcPhaseF","calcWallPhaseF"))
AddAction("Init"     , c("PhaseInit","calcWallPhaseF","BaseInit"))


# 	Outputs:
AddQuantity(name="Rho",	  unit="kg/m3")
AddQuantity(name="PhaseField",unit="1")
AddQuantity(name="U",	  unit="m/s",vector=T)
AddQuantity(name="P",	  unit="Pa")
#	Cases:
########SPHERE
AddSetting(name="Radius", default="0", comment='Use to create air bubble')
AddSetting(name="CenterX", default="0", comment='Location of air bubble')
AddSetting(name="CenterY", default="0", comment='Location of air bubble')
########RTI
AddSetting(name="Period", default="0", comment='Number of cells per cos wave')
AddSetting(name="Perturbation", default="0", comment='Size of wave perturbation, Perturbation*Period')
AddSetting(name="MidPoint", default="0", comment='height of RTI centerline')

# 	Three phase contact
AddSetting(name="theta", default="90", comment='Contact angle')

#	Inputs: For phasefield evolution
AddSetting(name="Density_h", comment='High density')
AddSetting(name="Density_l", comment='Low  density')
AddSetting(name="PhaseField_h", default=1, comment='PhaseField in Liquid')
AddSetting(name="PhaseField_l", default=0, comment='PhaseField gas')
AddSetting(name="PhaseField", 	   comment='Initial PhaseField distribution', zonal=T)
AddSetting(name="W", default=4,    comment='Anti-diffusivity coeff')
AddSetting(name="omega_phi", comment='one over relaxation time (phase field)')
AddSetting(name="M", omega_phi='1.0/(3*M+0.5)', default=0.02, comment='Mobility')
AddSetting(name="sigma", 		   comment='surface tension')
# 	Inputs: Fluid Properties
AddSetting(name="tau_l", comment='relaxation time (low density fluid)')
AddSetting(name="tau_h", comment='relaxation time (high density fluid)')
AddSetting(name="Viscosity_l", tau_l='(3*Viscosity_l)+0.5', default=0.16666666, comment='kinematic viscosity')
AddSetting(name="Viscosity_h", tau_h='(3*Viscosity_h)+0.5', default=0.16666666, comment='kinematic viscosity')

#Viscosity commented to alter to Guo Force Scheme
#AddSetting(name="Viscosity_l", tau_l='(3*Viscosity_l)', default=0.16666666, comment='kinematic viscosity')
#AddSetting(name="Viscosity_h", tau_h='(3*Viscosity_h)', default=0.16666666, comment='kinematic viscosity')
#AddSetting(name="S0", default=1.0, comment='Relaxation Param')
#AddSetting(name="S1", default=1.0, comment='Relaxation Param')
#AddSetting(name="S2", default=1.0, comment='Relaxation Param')
#AddSetting(name="S3", default=1.0, comment='Relaxation Param')
#AddSetting(name="S4", default=1.0, comment='Relaxation Param')
#AddSetting(name="S5", default=1.0, comment='Relaxation Param')
#AddSetting(name="S6", default=1.0, comment='Relaxation Param')
#	Inputs: Flow Properties
AddSetting(name="VelocityX", default=0.0, comment='inlet/outlet/init velocity', zonal=T)
AddSetting(name="VelocityY", default=0.0, comment='inlet/outlet/init velocity', zonal=T)
AddSetting(name="Pressure" , default=0.0, comment='inlet/outlet/init density', zonal=T)
AddSetting(name="GravitationX", default=0.0, comment='applied (rho)*GravitationX')
AddSetting(name="GravitationY", default=0.0, comment='applied (rho)*GravitationY')
AddSetting(name="BuoyancyX", default=0.0, comment='applied (rho-rho_h)*BuoyancyX')
AddSetting(name="BuoyancyY", default=0.0, comment='applied (rho-rho_h)*BuoyancyY')
# Globals - table of global integrals that can be monitored and optimized
AddGlobal(name="PressureLoss", comment='pressure loss', unit="1mPa")
AddGlobal(name="OutletFlux", comment='pressure loss', unit="1m2/s")
AddGlobal(name="InletFlux", comment='pressure loss', unit="1m2/s")
AddGlobal(name="TotalDensity", comment='Mass conservation check', unit="1kg/m3")

# For Contact angles
#AddNodeType(name="SouthWall", group="BOUNDARY")
