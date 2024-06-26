# Setting permissive access policy.
#  * This skips checks of fields being overwritten or read prematurely.
#  * Otherwise the model compilation was failing.
#  * This should be removed if the issue is fixed
SetOptions(permissive.access=TRUE)  ### WARNING


# Fluid Density Populations
AddDensity( name="f[0]", dx= 0, dy= 0, group="f")
AddDensity( name="f[1]", dx= 1, dy= 0, group="f")
AddDensity( name="f[2]", dx= 0, dy= 1, group="f")
AddDensity( name="f[3]", dx=-1, dy= 0, group="f")
AddDensity( name="f[4]", dx= 0, dy=-1, group="f")
AddDensity( name="f[5]", dx= 1, dy= 1, group="f")
AddDensity( name="f[6]", dx=-1, dy= 1, group="f")
AddDensity( name="f[7]", dx=-1, dy=-1, group="f")
AddDensity( name="f[8]", dx= 1, dy=-1, group="f")
# Pseudopotential field
AddField("psi", stencil2d=1)
# Stages and Actions
AddStage("BaseIteration", "Run",  save=Fields$group=="f", load=DensityAll$group=="f")
AddStage("calcPsi"    , save="psi", load=DensityAll$group == "f" )
AddStage("BaseInit"     , "Init", save=Fields$group=="f", load=DensityAll$group=="f")
AddAction("Iteration", c("BaseIteration","calcPsi"))
AddAction("Init"     , c("BaseInit",     "calcPsi"))
# Output Values
AddQuantity( name="Rho",unit="kg/m3"         )
AddQuantity( name="U"  ,unit="m/s"  ,vector=T)
AddQuantity( name="F"  ,unit="N"    ,vector=T)
AddQuantity( name="P"  ,unit="Pa"            )
AddQuantity( name="Psi",unit="1")
# Model Specific Parameters
AddSetting(name="G", 	 default=-1.0,    comment='interaction strength' )
AddSetting(name="T", 	 default=0.0585,  comment='effective temperature')
AddSetting(name="alpha", default=0.25, 	  comment='CS EoS parameter'     )
AddSetting(name="R"    , default=0.25,    comment='CS EoS parameter'     )
AddSetting(name="beta" , default=1 , 	  comment='CS EoS parameter'     )
AddSetting(name="kappa", default=0 , 	  comment='surface tension parameter' )
AddSetting(name="eps_0", default=2 , 	  comment='mechanical stability coef' )
AddSetting(name="betaforcing",default=1.0,comment='beta forcing scheme')
# Flow Properties:
AddSetting(name="omega", S7='1-omega', comment='one over relaxation time')
AddSetting(name="tempomega",default=1, comment='omega seems to get overwritten in preamble??')
AddSetting(name="nu", omega='1.0/(3*nu + 0.5)', default=0.16666666, comment='viscosity')
AddSetting(name="Velocity", default=0, comment='inlet/outlet/init velocity', zonal=T)
AddSetting(name="VelocityY", default=0, comment='init velocity in y dirn',   zonal=T)
AddSetting(name="Density",  default=1, comment='inlet/outlet/init density'  , zonal=T)
AddSetting(name="GravitationY", comment='Gravitation in the direction of y')
AddSetting(name="GravitationX", comment='Gravitation in the direction of x')
# Relaxation Properties
AddSetting(name="S0", default="0"          ,comment='MRT Sx')
AddSetting(name="S1", default="0"		   ,comment='MRT Sx')
AddSetting(name="S2", default="0"		   ,comment='MRT Sx')
AddSetting(name="S3", default="-.333333333",comment='MRT Sx')
AddSetting(name="S4", default="0"		   ,comment='MRT Sx')
AddSetting(name="S5", default="0"		   ,comment='MRT Sx')
AddSetting(name="S6", default="0"		   ,comment='MRT Sx')
AddSetting(name="S7", default="1.-omega"   ,comment='MRT Sx')
AddSetting(name="S8", default="1.-omega"   ,comment='MRT Sx')
# Globals - table of global integrals that can be monitored and optimized
AddGlobal(name="PressureLoss", comment='pressure loss', unit="1mPa")
AddGlobal(name="OutletFlux"  , comment='pressure loss', unit="1m2/s")
AddGlobal(name="InletFlux"   , comment='pressure loss', unit="1m2/s")
# Additional NodeTypes
AddNodeType(name="BottomSymmetry", group="BOUNDARY")
AddNodeType(name="TopSymmetry", group="BOUNDARY")
AddNodeType(name="RightSymmetry", group="BOUNDARY")
AddNodeType(name="EPressure", group="BOUNDARY")
AddNodeType(name="EVelocity", group="BOUNDARY")
AddNodeType(name="Solid", group="BOUNDARY")
AddNodeType(name="Wall", group="BOUNDARY")
AddNodeType(name="WPressure", group="BOUNDARY")
AddNodeType(name="WVelocity", group="BOUNDARY")
AddNodeType(name="MRT", group="COLLISION")
