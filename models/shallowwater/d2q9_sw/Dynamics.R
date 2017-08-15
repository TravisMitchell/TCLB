## Model for d2q9 SRT BKG-LBM
#     Density - performs streaming operation for us
#	

# Add the particle distribution functions as model Densities:
AddDensity( name="f[0]", dx= 0, dy= 0, group="f")
AddDensity( name="f[1]", dx= 1, dy= 0, group="f")
AddDensity( name="f[2]", dx= 0, dy= 1, group="f")
AddDensity( name="f[3]", dx=-1, dy= 0, group="f")
AddDensity( name="f[4]", dx= 0, dy=-1, group="f")
AddDensity( name="f[5]", dx= 1, dy= 1, group="f")
AddDensity( name="f[6]", dx=-1, dy= 1, group="f")
AddDensity( name="f[7]", dx=-1, dy=-1, group="f")
AddDensity( name="f[8]", dx= 1, dy=-1, group="f")

AddField( name="DensityF", stencil2d=1)
AddField( name="PressureF",stencil2d=1)

AddStage("BaseIter", "Run"  , save=Fields$group=="f", load=DensityAll$group=="f")

AddStage("calcDensity","calcDensityF",save=Fields$name=="DensityF",load=DensityAll$group=="f")
AddStage("calcPressure","calcPressureF",save=Fields$name=="DensityF" | Fields$name=="PressureF")

AddStage("Base", "Init", save=Fields$group=="f" 
			    | Fields$name=="DensityF"
			    | Fields$name=="PressureF" )

AddAction("Iteration", c("BaseIter", "calcDensity","calcPressure"))
AddAction("Init",      c("Base","calcPressure"))

# Add the quantities we wish to be exported
#    These quantities must be defined by a function in Dynamics.c
AddQuantity( name="U",unit="m/s", vector=TRUE )
AddQuantity( name="Rho",unit="kg/m3" )
AddQuantity( name="P", unit="Pa")

# Add the settings which describes system constants defined in a .xml file
AddSetting( name="omega", comment='inverse of relaxation time')
AddSetting( name="viscosity", omega='1.0/(3*viscosity+0.5)', default=0.16666666, comment='viscosity')
AddSetting( name="Velocity",default=0, comment='inlet/outlet/init velocity', zonal=TRUE)
AddSetting( name="Velocity_x",default=0, comment='inlet/outlet/init velocity in x', zonal=TRUE )
AddSetting( name="Velocity_y",default=0, comment='inlet/outlet/init velocity in y', zonal=TRUE )
AddSetting( name="Gravitation",default=0, comment='body/external acceleration', zonal=TRUE)
AddSetting( name="Density",default=1, comment='Density', zonal=TRUE)

# D2q9_sw specific parameters
AddSetting( name="Substrate", default=1, comment='On/off substrate', zonal=TRUE)
AddSetting( name="SlipLength", default=1, comment='Slip length', zonal=TRUE)
AddSetting( name="hc", default=0.05, comment='Critical height')
AddSetting( name="hstar", default=0.25, comment='Related to disjoining pressure')
AddSetting( name="n", default=9, comment='Exponents for disjoining pressure')
AddSetting( name="m", default=3, comment='Exponents for disjoining pressure')
AddSetting( name="ContactAngle", default=90, comment='Contact angle', zonal=TRUE)
AddSetting( name="SurfaceTension", default=0.01, comment='Surface Tension', zonal=TRUE)

# Initialisation
AddSetting( name="Rad", default="0")
AddSetting( name="cX",  default="0")
AddSetting( name="cY",  default="0")
AddSetting( name="dHigh",default="1")
AddSetting( name="dLow" ,default="0.05")

AddGlobal( name="Mass" )
