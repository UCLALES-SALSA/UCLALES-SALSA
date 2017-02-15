def get_ts_variable(fname,var_name,target_time):
	# Function for extracting scalar data from NetCDF file based on the given time value.
	# 
	#import numpy
	import netCDF4 as netcdf
	# Open the target NetCDF file
	ncid = netcdf.Dataset(fname,'r+')
	# Variable
	var = ncid.variables[var_name]
	# Time
	times = ncid.variables['time']
	#
	# Find the correct time value
	i=0
	out=-999
	for tt in times:
		if tt==target_time:
			out=var[i]
			break
		i=i+1
	# Close file
	ncid.close()
	return out

#
fname='/cygdrive/k/LES/UCLALES-SALSA/bin/ascos_30ccn_2D.ts.nc'
target_time=3600	# Time (s)

var_name='cfrac'	# Cloud fraction
cf=get_ts_variable(fname,var_name,target_time)
print "Cloud fraction=",cf

var_name='zc'	#  Cloud top height
dt=600				# Time step
zc1=get_ts_variable(fname,var_name,target_time)
zc2=get_ts_variable(fname,var_name,target_time+dt)
dzdt=(zc2-zc1)/600
print "Cloud top height tendency=",(zc2-zc1)/600,"m/s"
# Entrainment rate
div=1.5e-6	# Divergency (from NAMELIST)
E=dzdt-div*zc1
print "Entrainment rate",E,"m/s"
