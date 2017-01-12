#!/usr/bin/python
#
# Combine complete (name.########.nc) or profile (name.ps.########.nc) or temporal (name.ts.########.nc) statistics files from parallel runs. 
#
# Execution:
#	python combine.py rf01.ps
#	python combine.py rf01.ts
#	python combine.py rf01
#
import sys
import os
#
# Input files
# ==========
# Input file name
if len(sys.argv)<2: raise Exception('At least file name required!')
tmp=sys.argv[1]
infile=tmp.strip('\r')	# Mainly a windows problem
#
# At least one file must exists
if not os.path.isfile('%s.%04u%04u.nc' % (infile, 0, 0)):
	print '%s.%04u%04u.nc' % (infile, 0, 0)
	raise Exception( 'Data not found (expecting file %s.%04u%04u.nc)!' % (infile, 0, 0) )
#
# Limits for the data file indices
imax,jmax=0,0
while imax<1000 and os.path.isfile( '%s.%04u%04u.nc' % (infile, imax, 0) ): imax+=1
while jmax<1000 and os.path.isfile( '%s.%04u%04u.nc' % (infile, 0, jmax) ): jmax+=1
#
# More than one file required
if imax*jmax<=1: raise Exception('Not enough files: '+str(imax)+' x '+str(jmax)+'!')
#
# Test that everything is there
for i in range(imax):
	for j in range(jmax):
		if not os.path.isfile( '%s.%04u%04u.nc' % (infile, i, j) ): raise Exception('File %s.%04u%04u.nc not found!' % (infile, i, j)  )
#		
#
def reduce_ts_ps(infile,imax,jmax):
	# Statistics files
	#
	import numpy
	import netCDF4 as netcdf
	#from Scientific.IO import NetCDF as netcdf
	from shutil import copy2
	#
	# Different methods to combine data - edit these lists when new model outputs have been added
	#
	# From reduce.ncl
	# ------------------
	#maxnms = ['cfl','maxdiv','wmax','rlmax','bflxmx','bflxrmx','precip_m']
	#minnms = ['bflxmn','bflxrmn']
	#sumnms = ['wr_cs1','wr_cs2','wv_cs1','wv_cs2','wt_cs1','wt_cs2' \
	#	    ,'rt_cs1','rt_cs2','rl_cs1','rl_cs2','tv_cs1','tv_cs2' \
	#	    ,'tl_cs1','tl_cs2', 'w_cs1', 'w_cs2','cnt_cs1','cnt_cs2']
	#
	# Updated definitions
	# ----------------------
	# a) Averages - when the ouput should be the average over model domain
	#	Note: this is the 'default' method
	# b) Maximums - when the output should be the maximum over all model grid cells
	maxnms = ['zc','cfl','maxdiv','wmax','lmax']
	# c) Minimum - when the output should be the minimum over all model grid cells
	minnms = ['zb']
	# d) Sums - when the output is the sum over model grid cells
	#	Note: the output really depends on the size of the number of grid cells!
	sumnms = ['nrcnt','nccnt','nrain','pfrac','cnt_cs1','w_cs1','tl_cs1','tv_cs1','rt_cs1','rl_cs1','wt_cs1','wv_cs1','wr_cs1',\
	'cnt_cs2','w_cs2','tl_cs2','tv_cs2','rt_cs2','rl_cs2','wt_cs2','wv_cs2','wr_cs2']
	# e) Weighted averages (weights must be specified below) - needed for conditional averages
	#	Note: weight ~ conditinal coverage
	weighted = []
	# f) Variances (weights must be specified below) 
	#	Note 1: combining variances requires information about averages
	# 	Note 2: some units refer to standard deviation = sqrt(variance)?
	#	Note 3: some of the reported variances are not actually variances, but averages of squares
	#	           Example 1: lwp_var or 'Liquid-water path variance (kg/m^2)' is actually the average of LWP^2 (kg^2/m^4)
	variances = ['l_2','q_2','t_2','rflx2','sflx2',]
	# g) Third moments
	thirdmom = ['l_3','q_3']
	# h) Static - these should have the same values
	static = ['dn0','u0','v0','fsttm','lsttm','nsmp']
	#
	# *.ts.nc
	# ======
	# nrain         	{'units': '#/l', 'longname': 'Conditionally sampled rain number mixing ratio'}	<= Sum of rain drop concentrations over rainy grid cells (not avearge)
	# nrcnt         	{'units': '#', 'longname': 'Rain cell counts'}										<= Number of rainy grid cells
	# pfrac			{'units': '-', 'longname': 'Precipitation fraction'}									<= Number of precipitating grid cells (-999,0 means zero)
	# lwp_var		{'units': 'kg/m^2', 'longname': 'Liquid-water path variance'}						<= Actually just average of LWP^2
	# CCN         		{'units': '#', 'longname': 'Cloud condensation nuclei'}								<= Actually concentration (#/cc)
	#
	# *.ps.nc
	# ======
	# evap			{'units': 's^-1', 'longname': 'Net evap  of rain-water'}
	# frc_prc		{'units': '-', 'longname': 'Conditionally sampled rain fraction'}
	# frc_ran        	{'units': '-', 'longname': 'Rain water fraction'}
	# hst_srf         	{'units': '-', 'longname': 'Histogram of surface rain rates'}
	# fsttm         	{'units': 'kg/m^3', 'longname': 'First sample time'}
	# lsttm         	{'units': 'kg/m^3', 'longname': 'Basic state density'}							<= the last time point (s)
	# nsmp         		{'units': 'kg/m^3', 'longname': 'Basic state density'}							<= number of time steps for the average
	# prc_prc         	{'units': 'm/s', 'longname': 'Conditionally sampled rain rate'}
	# rflx2       		{'units': 'W/m^2', 'longname': 'Variance of total radiative flux'}			<= it is variance
	# sflx2         	{'units': 'W/m^2', 'longname': 'Variance of shortwave radiative flux'}		<= it is variance
	# l_2        		{'units': '-', 'longname': 'Variance of liquid'}									<= it is variance
	# l_3         		{'units': '-', 'longname': 'Third moment of liquid'}								<= it is the third central moment
	# q_2         		{'units': '-', 'longname': 'Variance of total water'}							<= it is variance
	# q_3         		{'units': '-', 'longname': 'Third moment of total water'}						<= it is the third central moment
	# t_2         		{'units': 'K^2', 'longname': 'Variance of theta'}								<= it is variance
	# t_3         		{'units': 'K^3', 'longname': 'Third moment of theta'}							<= average of (t-avg(v))^3, which must be an error
	# u_2         		{'units': 'm^2/s^2', 'longname': 'Variance of u wind'}							<= no, it is average of u^2
	# v_2         		{'units': 'm^2/s^2', 'longname': 'Variance of v wind'}							<= no, it is average of v^2
	# w_2         		{'units': 'm^2/s^2', 'longname': 'Variance of w wind'}							<= no, it is average of w^2
	# w_3         		{'units': 'm^3/s^3', 'longname': 'Third moment of w wind'}					<= it is the third raw moment of w (i.e. average of w^3)
	# cs1         		{'units': '-', 'longname': 'Conditionally sampled fraction of flow'}			<= average number of flag values
	# cnt_cs1			{'units': '#', 'longname': 'Sum of I_cs1'}										<= number of flag values
	# w_cs1         	{'units': 'm', 'longname': 'Conditional average of w over cs1'}				<= sum of flagged values
	# ...
	# wr_cs1         	{'units': 'W/m^2', 'longname': 'Covariance of wr_t flux and cs1'}			<= sum of flagged w*r values
	# wt_cs1         	{'units': 'W/m^2', 'longname': 'Covariance of wtheta_l flux and cs1'}		<= sum of flagged w*thetal_l values
	# wv_cs1         	{'units': 'W/m^2', 'longname': 'Covariance of wtheta_v flux and cs1'}		<= sum of flagged w*thetal_v values
	# ...
	#
	nfiles=imax*jmax
	#
	# Calculations
	# ============
	# 1) Generate target
	src='%s.%04u%04u.nc' % (infile, 0, 0)
	dst='%s.nc' % (infile)
	copy2(src, dst)
	#
	# Open the target NetCDF file
	ncid = netcdf.Dataset(dst,'r+')
	#ncid = netcdf.NetCDFFile(dst,'r+')
	#
	# Can use variables:
	#	NPTS = (nxp-4)*(nyp-4)
	#	NZ = nzp
	#
	# Variable list - without dimensions
	var_list=ncid.variables.keys()	# Variables
	dim_list=ncid.dimensions.keys()	# Dimensions
	#
	print ' '
	print 'Generating file '+dst+' from '+str(nfiles)+' files'
	print ' '
	print '%-9s  %3s  %s' % ('name','op.','longname')
	print '---------------------------------------------'
	for name in var_list:
		# Skip dimension variables
		if name in dim_list: continue
		#
		# What operation to apply
		# => Avg, min or max
		if name in maxnms:
			# Maximum
			op='max'
		elif name in minnms:
			# Minimum
			op='min'
		elif name in sumnms:
			# Sum
			op='sum'
		elif name in weighted:
			# Weighted averages
			op='wgh'
			if name=='zb' or name=='zc':
				# Cloud base or top heights: average weighted by cloud fraction
				wname='cfrac'
			else:
				print 'Unknown method for '+name+': default (avg) used!'
				op='avg'
		elif name in variances:
			# Average of variances: need both variance and average
			# 	Note: this is for true variances and not for sqaures
			op='var'
			if name=='lwp_var':	# Liquid water path
				aname='lwp_bar'
			elif name=='u_2':		# u wind
				aname='u'
			elif name=='v_2':		# v wind
				aname='v'
			elif name=='w_2':		# w wind
				aname='w'
			elif name=='t_2':		# temperature
				aname='t'
			elif name=='l_2':		# liquid water
				aname='l'
			elif name=='q_2':		# humidity
				aname='q'
			elif name=='rflx2':	# Total radiative flux
				aname='rflx'
			elif name=='sflx2':	# Shortwave radiative flux
				aname='sflx'
			else:
				print 'Unknown method for '+name+': default (avg) used!'
				op='avg'
		elif name in thirdmom:
			# Third momenets: need the second (variance) and first (average) moments
			# 	Note: this is for true third moments and not for cubes
			op='thr'
			if name=='l_3':
				vname='l_2'
				aname='l'
			elif name=='q_3':
				vname='q_2'
				aname='q'
			else:
				print 'Unknown method for '+name+': default (avg) used!'
				op='avg'
		elif name in static:
			# Should have constant values
			op='con'
		elif name[-3:]=='_ic' or name[-4:]=='_int':
			# Cloud region: these numbers are actually averages of column averages, 
			# so cloud fraction is a proper weight. This should be replaced the overall 
			# average (then shoul use the number of cloudy grid cells)!
			op='wgh'
			wname='cfrac'
			#if 'nccnt' in var_list:
			#	# This is update from December 2016
			#	wname='nccnt'
			#else:
			#	# Use cloud fraction (better than nothing!)
			#	wname='cfrac'
		elif 0>1: #name[-3:]=='_oc':
			# Outside cloud region: should be weighted by the number of clear grid cells
			#if 'nccnt' in var_list:
			#	# This is update from December 2016
			#	op='wgh_neg'
			#	wname='nccnt'
			#else:
			op='avg'	# Simple average, because most grid cells are anyhow outside clouds
		else:
			# Average (default)
			op='avg'
		#
		# Object
		obj = ncid.variables[name]
		#val = obj.getValue()
		val = obj[:] #.getValue()
		#
		# Calculations for the given variable
		tmp=obj.__dict__	# Dictionary of information
		print '%-9s  %3s   %s (%s)' % (name,op,tmp['longname'],tmp['units'])
		#
		info_prints=0
		#
		# All files
		for i in range(imax):
			for j in range(jmax):
				src='%s.%04u%04u.nc' % (infile, i, j)
				#ncid_src = netcdf.NetCDFFile(src,'r')
				ncid_src = netcdf.Dataset(src,'r')
				# Object
				obj_src = ncid_src.variables[name]
				val_src = obj_src[:] #.getValue()
				#
				# Are there flagged values?
				if numpy.any(val_src==-999):
					# Flagged values found
					if name=='pfrac':
						# Precipitation fraction (actually the number of precipitating grid cells): -999 = 0
						val_src[val_src==-999]=0
						if info_prints==0: print '	Flagged values (-999) set to zero!'
						info_prints+=1
					else:
						if info_prints==0: print '	Flagged values (-999) found!'
						info_prints+=1
					#
					if (i+j==0): val=val_src
				#
				if op=='max':
					# Maximum (i=j=0 is already the value for val)
					if i+j>0: val=numpy.maximum(val,val_src)
				elif op=='min':
					# Minimum
					if i+j>0:	val=numpy.minimum(val,val_src)
				elif op=='sum':
					# Sum
					if i+j>0:	val+=val_src
				elif op=='wgh':
					# Weighted average
					obj_src = ncid_src.variables[wname]	# Weight
					val_weight = obj_src[:] #.getValue()
					if i+j==0:
						val=val_src*val_weight
						weight=val_weight
					else:
						val+=val_src*val_weight
						weight+=val_weight
				elif op=='var':
					# Variance
					obj_src = ncid_src.variables[aname]	# Average
					val_avg = obj_src[:] #.getValue()
					if i+j==0:
						val=(val_src+val_avg**2)/float(nfiles)	# Variance
						avg=val_avg/float(nfiles)						# Average
					else:
						val+=(val_src+val_avg**2)/float(nfiles)
						avg+=val_avg/float(nfiles)
				elif op=='thr':
					# Third central moment
					obj_src = ncid_src.variables[aname]	# Average
					val_avg = obj_src[:] #.getValue()
					obj_src = ncid_src.variables[vname]	# Variance
					val_var = obj_src[:] #.getValue()
					if i+j==0:
						val=(val_src+3.0*val_avg*val_var+val_avg**3)/float(nfiles)
						var=(val_var+val_avg**2)/float(nfiles)	# Variance
						avg=val_avg/float(nfiles)						# Average
					else:
						val+=(val_src+3.0*val_avg*val_var+val_avg**3)/float(nfiles)
						var+=(val_var+val_avg**2)/float(nfiles)
						avg+=val_avg/float(nfiles)
				elif op=='con':
					# Should be constant (allow 1e-10 absolute difference)
					if i+j==0:
						val=val_src
					elif numpy.amax( numpy.absolute(val-val_src) )>1e-10:
						print 'Values of '+name+' change (expecting constant values)!'
				else:
					# Average
					if i+j==0:
						val=val_src/float(nfiles)
					else:
						val+=val_src/float(nfiles)
				# Close current file
				ncid_src.close()
		#
		if op=='wgh':
			weight[weight==0]=1e-20	# Avoid dividing by zero
			val/=weight
		if op=='var': val-=avg**2
		if op=='thr': val-=3.0*avg*var+avg**3
		#
		# Save updated data
		#obj.assignValue(val)
		obj[:]=val
		#ncid.flush()
		ncid.sync()
	#	
	# Close file
	ncid.close()
#
#
def reduce_full(infile,imax,jmax):
	# Full data files
	#	-Shapes cannot be modified so a new NetCDF file and variables must be created
	import numpy
	import netCDF4 as netcdf
	#from Scientific.IO import NetCDF as netcdf
	from shutil import copy2
	#
	nfiles=imax*jmax
	#
	#
	def find_sequence(base,seq):
		for start in range(len(base)):
			if abs(base[start]-seq[0])<1e-5:
				# Old: if base[start]==seq[0]:
				for end in range(len(seq)):
					if abs(base[start+end]-seq[end])>1e-5:
						# Old: base[start+end]!=seq[end]
						end=-1
						break
				if end>=0: return start
		raise Exception('Range not found!')
	#
	#
	def set_values(out,inp,pointer):
		if len(pointer)==1:
			# Vector
			out[pointer:pointer+len(inp)]=inp
		elif len(pointer)==2:
			# 2D
			out[pointer[0]:pointer[0]+inp.shape[0],pointer[1]:pointer[1]+inp.shape[1]]=inp
		elif len(pointer)==3:
			# 3D
			out[pointer[0]:pointer[0]+inp.shape[0],pointer[1]:pointer[1]+inp.shape[1],pointer[2]:pointer[2]+inp.shape[2]]=inp
		elif len(pointer)==4:
			# 4D
			out[pointer[0]:pointer[0]+inp.shape[0],pointer[1]:pointer[1]+inp.shape[1],pointer[2]:pointer[2]+inp.shape[2],pointer[3]:pointer[3]+inp.shape[3]]=inp
		elif len(pointer)==5:
			# 5D
			out[pointer[0]:pointer[0]+inp.shape[0],pointer[1]:pointer[1]+inp.shape[1],pointer[2]:pointer[2]+inp.shape[2],pointer[3]:pointer[3]+inp.shape[3],pointer[4]:pointer[4]+inp.shape[4]]=inp
		else:
			print len(pointer)
			raise Exception('Unknown dimensions')
	#
	#
	# Calculations
	# ============
	# Existing NetCDF file
	src='%s.%04u%04u.nc' % (infile, 0, 0)
	#ncid_src = netcdf.NetCDFFile(src,'r')
	ncid_src = netcdf.Dataset(src,'r')
	# Variable and dimension list
	var_list=ncid_src.variables.keys()		# Variables
	dim_list=ncid_src.dimensions.keys()	# Dimensions
	#
	# Output NetCDF file
	dst='%s.nc' % (infile)
	ncid_dst = netcdf.Dataset(dst,'w',format='NETCDF3_CLASSIC')
	#ncid_dst = netcdf.NetCDFFile(dst,'w')
	# Copy global attributes
	print 'Copying global attributes...'
	for att in ncid_src.ncattrs():
		ncid_dst.setncattr(att,ncid_src.getncattr(att))
	#
	#skip=['close', 'createDimension', 'createVariable', 'flush', 'sync','__class__'] # These are always included, but should not be saved!
	#for att in dir(ncid_src):
	#	if att in skip: continue
	#	setattr(ncid_dst,att,getattr(ncid_src,att))
	#
	print 'Done'
	ncid_src.close()
	#
	# a) Dimensions
	# Generate a map based on dimensions
	indices=numpy.empty([len(dim_list),imax,jmax], dtype=int)
	sizes=numpy.zeros(len(dim_list))
	k=0
	print 'Creating dimensions...'
	for name in dim_list:
		# All files
		for i in range(imax):
			for j in range(jmax):
				# Source
				src='%s.%04u%04u.nc' % (infile, i, j)
				ncid_src = netcdf.Dataset(src,'r')
				#ncid_src = netcdf.NetCDFFile(src,'r')
				# Source data
				obj_src = ncid_src.variables[name]
				val_src = obj_src[:] #.getValue()
				#
				if i==0 and j==0:
					# Initialize output
					val=val_src
					indices[k,i,j]=0
					#
					# Information about the dimension
					old_len=len(val)			# Original length (in file 0000 0000)
					var_info=obj_src.__dict__	# Attributes
					typ=obj_src.dtype
					#typ=obj_src.typecode() 	# Type code
				elif numpy.amin(val_src)>numpy.amax(val):
					# Append after previous data
					indices[k,i,j]=len(val)
					val=numpy.append(val,val_src)
				elif numpy.amax(val_src)<numpy.amin(val):
					# Append before => error!
					raise Exception('Monotonic order expected!')
				else:
					# There should be a matching sequence
					start=find_sequence(val,val_src)
					indices[k,i,j]=start
				ncid_src.close()
		# New size of the dimension
		sizes[k]=len(val)
		#
		# Save updated dimensions
		# a) Create dimension
		ncid_dst.createDimension(name,len(val))
		# b) Save value as that of a variable
		id=ncid_dst.createVariable(name,typ,(name,))
		#id.assignValue(val)
		id[:]=val
		# Copy attributes
		#for att in ncid_src.ncattrs(): id.setncattr(att,ncid_src.getncattr(att))
		#print var_info.keys()
		#for att in var_info.keys(): id.setncattr(att,var_info[att])	#setattr(id,att,var_info[att])
		for att in var_info.keys(): setattr(id,att,var_info[att])
		#
		print '  %-8s %4u => %-4u  %s' % (name, old_len, len(val), var_info)
		k+=1
	print 'Done'
	#
	# b) Variables
	print 'Creating varibles...'
	for name in var_list:
		# Skip dimension variables
		if name in dim_list: continue
		#
		# There are variables that are independent of processor (e.g. u0, v0, w0,...)
		# => Should be constants
		IsConst=False
		#
		# All files
		for i in range(imax):
			for j in range(jmax):
				# Source
				src='%s.%04u%04u.nc' % (infile, i, j)
				ncid_src = netcdf.Dataset(src,'r')
				#ncid_src = netcdf.NetCDFFile(src,'r')
				obj_src = ncid_src.variables[name]
				val_src = obj_src[:] #.getValue()	
				#
				if i==0 and j==0:
					# Generate new dimension by reshaping single processor data
					#
					# Information about the variable
					var_info=obj_src.__dict__	# Attributes
					typ=obj_src.dtype
					#typ=obj_src.typecode() 	# Type code
					dims=obj_src.dimensions	# Dimensions
					#
					# Index to current dimensions
					dim_ind=[]
					for dim in dims: dim_ind.append(dim_list.index(dim))
					#
					# Update output size
					val=numpy.copy(val_src)
					val.resize(sizes[dim_ind])
					#
					# Current starting index should be unique
					#ind_ref=indices[dim_ind,0,0]	# Starting index
				elif (indices[dim_ind,0,0]==indices[dim_ind,i,j]).all():
					# The data will be overwritten
					# => OK when the data is aways identical
					if (i+j==1 or IsConst) and (val_src==val).all():
						# Identical data => can be constant
						IsConst=True
						continue
					# Data will be overwritten => problem
					raise Exception('Overlapping data for '+name+'!')
				#
				# Set the data
				ind=indices[dim_ind,i,j]	# Starting index
				set_values(val,val_src,ind)
		#
		# Create variable
		id=ncid_dst.createVariable(name,typ,dims)
		#id.assignValue(val)
		id[:]=val
		# Copy attributes
		for att in var_info.keys(): setattr(id,att,var_info[att])
		#
		# Save data
		#ncid_dst.flush()
		ncid_dst.sync()
		#
		# Print info
		if IsConst:
			# Constant data
			print '  *%-8s %30s  %s' % (name,dims,var_info)
		else:
			print '  %-8s  %30s  %s' % (name,dims,var_info)
	#
	print '*Identical data in each netCDF file'
	print 'Done'
	#	
	# Close file
	ncid_dst.close()
#
#
def nc_PrintFileContent(filename,filename_out=''):
	# Print information on all variables from a file
	import netCDF4 as netcdf
	#from Scientific.IO import NetCDF as netcdf
	#
	# Save results to a file
	if len(filename_out)>0:
		# Save output to a file - the old file is overwritten!
		out_file=open(filename_out,'w')
	#
	# Open
	ncid = netcdf.Dataset(filename,'r')
	#ncid = netcdf.NetCDFFile(filename,'r')
	#
	# Dimension variables
	allDimNames = ncid.dimensions.keys()
	ndims=len(allDimNames)
	#
	# List of all variables
	vars=ncid.variables
	#
	# Header line
	if len(filename_out)==0:
		# Print to screen
		print '%3s %10s %10s %5s   %5s' % ('No.','dimens.','name','size','attributes')
	else:
		print >>out_file, '%3s %10s %10s %5s   %5s' % ('No.','dimens.','name','size','attributes')
	#
	# 1) Print dimension variables
	i=0				# Counter
	for name in allDimNames:
		# Object
		obj = ncid.variables[name]
		#
		# Variable information: no dimension id, but size
		dstr='-';
		dimlen=obj.shape
		dimlen_str=str(dimlen[0])
		#
		# Attribute information
		attrs=str( obj.__dict__ )
		#
		if len(filename_out)==0:
			# Print to screen
			print '%3u %10s %10s %5s   %5s' % (i,dstr,name,dimlen_str,attrs)
		else:
			# Print to file
			print >>out_file, '%3u %10s %10s %5s   %5s' % (i,dstr,name,dimlen_str,attrs)
		i+=1
	#
	# 2) Print data
	for name in sorted(vars.iterkeys()):
		# Skip dimension variables
		if name in allDimNames: continue
		#
		# Object
		obj = ncid.variables[name]
		#
		# Variable information: no size, but dimension id list
		dimlen_str='';
		# All dimensions
		dims=obj.dimensions
		dstr=''
		for dim in dims:
			dstr+=str( allDimNames.index(dim)  )+','
		dstr=dstr[0:len(dstr)-1]
		#
		# Attribute information
		attrs=str( obj.__dict__ )
		#
		if len(filename_out)==0:
			# Print to screen
			print '%3u %10s %10s %5s   %5s' % (i,dstr,name,dimlen_str,attrs)
		else:
			# Print to file
			print >>out_file, '%3u %10s %10s %5s   %5s' % (i,dstr,name,dimlen_str,attrs)
		i+=1
	#
	# Close file
	ncid.close()
	#
	# Close output
	if len(filename_out): out_file.close()
#
#
if '.ts' in infile or '.ps' in infile:
	# Time or profile statistics
	reduce_ts_ps(infile,imax,jmax)
else:
	# Full 4D data
	reduce_full(infile,imax,jmax)
#
# Save information about file content?
# nc_PrintFileContent(infile+'.nc',infile+'.nc.txt')
#
# End of function combine.py
