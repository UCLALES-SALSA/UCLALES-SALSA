#!/usr/bin/python
#
# Functions for post-processing UCLALES-SALSA outputs from parallel runs.
#	1a) Combine 4D outputs (name.########.nc)
# 	1b) Combine 3D column outputs (name.cs.########.nc)
#	2) Average profile (name.ps.########.nc) or temporal (name.ts.########.nc) statistics
#
# Execution:
#	One or more files
#		python combine.py rf01
#		python combine.py rf01.ts
#		python combine.py rf01.ts rf01.ps rf01.cs
#	One file and list if variables to be included in the outputs
#		python combine.py rf01 l q t
#
# Last update: 13.4.2017 Tomi Raatikainen
#
# Required modules
#		Python 2.7.10
#
# Edit this tag when outputs are changed
version='combine_1.0.0'
#
import sys
import os
import time
#
# Time when data is post-processed
date=time.asctime()
#
def reduce_ts_ps(infile,imax,jmax,var_list):
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
	maxnms = ['zc','cfl','maxdiv','wmax','lmax','imax','smax']
	# c) Minimum - when the output should be the minimum over all model grid cells
	minnms = ['zb']
	# d) Sums - when the output is the sum over model grid cells
	#	Note: the output really depends on the number of grid cells!
	sumnms = ['nrcnt','nccnt','cnt_cs1','cnt_cs2','nicnt','nscnt']
	# e) Weighted averages (weights must be specified below) - needed for conditional averages
	#	Note: weight ~ conditinal coverage
	weighted = []
	# f) Variances (weights must be specified below)
	#	Note 1: combining variances requires information about averages
	# 	Note 2: some units refer to standard deviation = sqrt(variance)?
	#	Note 3: some of the reported variances are not actually variances, but averages of squares
	variances = ['l_2','q_2','t_2','rflx2','sflx2',]
	# g) Third moments
	thirdmom = ['l_3','q_3']
	# h) Static - these should have the same values
	static = ['dn0','u0','v0','fsttm','lsttm','nsmp']
	#
	nfiles=imax*jmax
	#
	# Calculations
	# ============
	# 1) Generate target
	src='%s.%04u%04u.nc' % (infile, 0, 0)
	dst='%s.nc' % (infile)
	if len(var_list)==0:	# Variables not filtered
		# Just copy the first file
		copy2(src, dst)
		#
		# Add information about post processing? Not possible for NetCDF3!
	else:
		# Generate new file with same global attributes and dimensions, but with selected variables
		ncid_src = netcdf.Dataset(src,'r')
		dim_list=ncid_src.dimensions.keys()	# Dimensions
		fmt=ncid_src.data_model	# Typically NETCDF3_CLASSIC, but this cannot handle 2+ Gb files!
		# Output NetCDF file
		dst='%s.nc' % (infile)
		ncid = netcdf.Dataset(dst,'w',format=fmt)
		#
		# Copy global attributes
		for att in ncid_src.ncattrs():
			ncid.setncattr(att,ncid_src.getncattr(att))
		#
		# Add information about post processing
		ncid.setncattr('PP_version',version)
		ncid.setncattr('PP_date',date)
		#
		# Copy dimensions and relevant variables
		for name in ncid_src.variables.keys():
			# Source data
			obj_src = ncid_src.variables[name]
			val = obj_src[:]	# Values
			typ=obj_src.dtype	# Type
			#
			if name in dim_list: 	# Create dimension
				# Target
				ncid.createDimension(name,len(val))
				# Save value as that of a variable
				id=ncid.createVariable(name,typ,(name,))
				id[:]=val
			elif name in var_list:	# Create variable
				dims=obj_src.dimensions	# Dimensions
				id=ncid.createVariable(name,typ,dims)
				id[:]=val
			else:
				# Skip
				continue
			#
			# Copy attributes
			var_info=obj_src.__dict__	# Attributes
			for att in var_info.keys(): setattr(id,att,var_info[att])
		#
		# Save data
		ncid.sync()
		#
		ncid_src.close()
		ncid.close()
	#
	# Open the target NetCDF file
	ncid = netcdf.Dataset(dst,'r+')
	#
	#
	# Can use variables:
	#	NPTS = (nxp-4)*(nyp-4)
	#	NZ = nzp
	#
	# Information about file version: 1.0 is the original one and later versions are greater than that
	if 'IO_version' in ncid.ncattrs():
		io=getattr(ncid,'IO_version')
	else:
		# IO before version v1.0.5
		io=1.0
		sumnms+=['nrain','pfrac']
		sumnms+=['w_cs1','tl_cs1','tv_cs1','rt_cs1','rl_cs1','wt_cs1','wv_cs1','wr_cs1',\
			'w_cs2','tl_cs2','tv_cs2','rt_cs2','rl_cs2','wt_cs2','wv_cs2','wr_cs2']
	#
	# Variable and dimension lists
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
			# No weighted averages?
			print 'Unknown method for '+name+': default (avg) used!'
			op='avg'
		elif name in variances:
			# Average of variances: need both variance and average
			# 	Note: this is for true variances, but not for sqaures
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
			# 	Note: this is for true third moments, but not for cubes
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
		else:
			# Average (default)
			op='avg'
		#
		# Object
		obj = ncid.variables[name]
		val = obj[:]
		#
		# Calculations for the given variable
		tmp=obj.__dict__	# Dictionary of information
		print '%-9s  %3s   %s (%s)' % (name,op,tmp['longname'],tmp['units'])
		#
		info_prints=0
		#
		set_missing=-999
		#
		# All files
		for i in range(imax):
			for j in range(jmax):
				src='%s.%04u%04u.nc' % (infile, i, j)
				ncid_src = netcdf.Dataset(src,'r')
				# Object
				obj_src = ncid_src.variables[name]
				val_src = obj_src[:]
				#
				# Are there flagged values?
				if numpy.any(val_src==-999):
					# Flagged values found
					if name=='pfrac':
						# Precipitation fraction: missing values can be set to zero
						val_src[val_src==-999]=0
						if info_prints==0: print '	Flagged values (-999) set to zero!'
						info_prints+=1
					elif name=='zb':
						# Cloud base height: set to a large number (and convert back to -999 when all values are unavailable)
						set_missing=99999	# Will be changed to -999 if any left
						val_src[val_src==-999]=set_missing
						if info_prints==0: print '	Flagged values (-999) found!'
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
					val_weight = obj_src[:]
					if i+j==0:
						val=val_src*val_weight
						weight=val_weight
					else:
						val+=val_src*val_weight
						weight+=val_weight
				elif op=='var':
					# Variance
					obj_src = ncid_src.variables[aname]	# Average
					val_avg = obj_src[:]
					if i+j==0:
						val=(val_src+val_avg**2)/float(nfiles)	# Variance
						avg=val_avg/float(nfiles)						# Average
					else:
						val+=(val_src+val_avg**2)/float(nfiles)
						avg+=val_avg/float(nfiles)
				elif op=='thr':
					# Third central moment
					obj_src = ncid_src.variables[aname]	# Average
					val_avg = obj_src[:]
					obj_src = ncid_src.variables[vname]	# Variance
					val_var = obj_src[:]
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
		# May need to convert some values back to -999
		if numpy.any(val==set_missing):
			val[val==set_missing]=-999
		#
		# Save updated data
		obj[:]=val
		ncid.sync()
	#
	# Close file
	ncid.close()
#
#
def reduce_full(infile,imax,jmax,var_list):
	# Full data files
	#	-Shapes cannot be modified so a new NetCDF file and variables must be created
	import numpy
	import netCDF4 as netcdf
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
			out[pointer[0]:pointer[0]+len(inp)]=inp
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
	ncid_src = netcdf.Dataset(src,'r')
	# Variable and dimension lists
	if len(var_list)==0: var_list=ncid_src.variables.keys()		# Variables
	dim_list=ncid_src.dimensions.keys()	# Dimensions
	fmt=ncid_src.data_model	# netCDF data model version (typically NETCDF3_CLASSIC)
	# **********************************************************
	# Note: NETCDF3_CLASSIC does not handle 2+ Gb files!
	b = os.path.getsize(src)	# File size in bytes
	if nfiles*b>2e9:
		print 'Warning: file size exceeds 2 GB (total approx. %.2f GB) - changing to NetCDF4' % (nfiles*b/1073741824.0)
		fmt='NETCDF4'
	# **********************************************************
	#
	# Output NetCDF file
	dst='%s.nc' % (infile)
	ncid_dst = netcdf.Dataset(dst,'w',format=fmt)
	# Copy global attributes
	print 'Copying global attributes...'
	for att in ncid_src.ncattrs():
		ncid_dst.setncattr(att,ncid_src.getncattr(att))
		print '  %-15s  %s' %(att,ncid_src.getncattr(att))
	# Add information about post processing
	ncid_dst.setncattr('PP_version',version)
	ncid_dst.setncattr('PP_date',date)
	#
	print 'Done'
	ncid_src.close()
	#
	# a) Dimensions
	# Generate a map based on dimensions
	indices=numpy.empty([len(dim_list),imax,jmax], dtype=int)
	sizes=numpy.zeros(len(dim_list),dtype=int)
	k=0
	print 'Creating dimensions...'
	for name in dim_list:
		# All files
		for i in range(imax):
			for j in range(jmax):
				# Source
				src='%s.%04u%04u.nc' % (infile, i, j)
				ncid_src = netcdf.Dataset(src,'r')
				# Source data
				obj_src = ncid_src.variables[name]
				val_src = obj_src[:]
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
				del val_src
				ncid_src.close()
		# New size of the dimension
		sizes[k]=len(val)
		#
		# Save updated dimensions
		# a) Create dimension
		ncid_dst.createDimension(name,len(val))
		# b) Save value as that of a variable
		id=ncid_dst.createVariable(name,typ,(name,))
		id[:]=val
		# Copy attributes
		for att in var_info.keys(): setattr(id,att,var_info[att])
		#
		#print '  %-8s %4u => %-4u  %s' % (name, old_len, len(val), var_info)
		print '  %-8s %4u => %-4u  %s (%s)' % (name,old_len, len(val),var_info['longname'],var_info['units'])

		k+=1
	print 'Done'
	#
	# b) Variables
	print 'Creating varibles...'
	for name in var_list:
		# Skip dimension variables
		if name in dim_list: continue
		#if name in ignore_list: continue
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
				obj_src = ncid_src.variables[name]
				val_src = obj_src[:] #.getValue()
				#
				if i==0 and j==0:
					# Generate new dimension by reshaping single processor data
					#
					# Information about the variable
					var_info=obj_src.__dict__	# Attributes
					typ=obj_src.dtype	# Type code
					dims=obj_src.dimensions	# Dimensions
					#
					# Index to current dimensions
					dim_ind=[]
					for dim in dims: dim_ind.append(dim_list.index(dim))
					#
					# Update output size
					val=numpy.copy(val_src)
					val.resize(sizes[dim_ind])
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
				# Clean
				del val_src
		#
		# Create variable
		#print 'creating',name,typ,dims,val.shape
		id=ncid_dst.createVariable(name,typ,dims)
		id[:]=val
		# Copy attributes
		for att in var_info.keys(): setattr(id,att,var_info[att])
		#
		# Save data
		ncid_dst.sync()
		#
		# Print info
		if IsConst:
			# Constant data
			#print '  *%-8s %30s  %s' % (name,dims,var_info)
			print '  *%-8s %30s  %s (%s)' % (name,dims,var_info['longname'],var_info['units'])
		else:
			#print '  %-8s  %30s  %s' % (name,dims,var_info)
			print '  %-8s  %30s  %s (%s)' % (name,dims,var_info['longname'],var_info['units'])
	#
	print '*Identical data in each netCDF file'
	print 'Done'
	#	
	# Close file
	ncid_dst.close()
#
#
def examine_fname(infile):
	# Examine if the given file name pattern represents a valid output file
	import os
	#
	# At least one file must exists
	if not os.path.isfile('%s.%04u%04u.nc' % (infile, 0, 0)): return 0,0
	#
	# Limits for the data file indices
	imax,jmax=0,0
	while imax<1000 and os.path.isfile( '%s.%04u%04u.nc' % (infile, imax, 0) ): imax+=1
	while jmax<1000 and os.path.isfile( '%s.%04u%04u.nc' % (infile, 0, jmax) ): jmax+=1
	#
	# More than one file required
	if imax*jmax<=1: return 0,0
	#
	# Test that everything is there
	for i in range(imax):
		for j in range(jmax):
			if not os.path.isfile( '%s.%04u%04u.nc' % (infile, i, j) ): return 0,0
	#
	# Valid file found
	return imax,jmax
#
#
def examine_var_name(file,var):
	# Examine if variable can be found from given NetCDF file (given file name pattern)
	import netCDF4 as netcdf
	# File name
	fname=('%s.%04u%04u.nc' % (file, 0, 0))
	# OPen NetCDF file
	ncid = netcdf.Dataset(fname,'r')
	# Compare name with variable list
	if var in ncid.variables.keys():
		# Found
		ncid.close()
		return True
	# Not found
	ncid.close()
	return False
#
#
# Input files
# ==========
# Input file name must be given
if len(sys.argv)<2: raise Exception('At least file name required!')
#
# Input options:
#	a) One or more file name patterns
#	b) Single file name pattern and a list of data to be included
names=[]		# List of file names
variables=[]	# List of variable (empty=all)
i=1
for tmp in sys.argv:
	# The first argument is not an input
	if tmp==sys.argv[0]: continue
	#
	infile=tmp.strip('\r')	# Mainly a windows problem
	#
	# The first argument must be a file name. The second argument determines the following arguments after that.
	# If the second argument is file name, variable name list is empty.
	if i==1 or (i>2 and len(variables)==0):
		# Expect file name
		imax,jmax=examine_fname(infile)
		if imax==0:
			raise Exception( 'File %s.%04u%04u.nc not found!' % (infile, 0, 0) )
		# It is a valid file
		names.append(infile)
	elif i==2:
		# The second argument can be either file or variable name
		imax,jmax=examine_fname(infile)
		if imax>0:
			# It is a valid file
			names.append(infile)
		elif examine_var_name(names[0],infile):
			# It is a valid variable
			variables.append(infile)
		else:
			raise Exception( 'Data not found (file or variable name %s)!' % (infile) )
	else:
		# Expecting variable
		if not examine_var_name(names[0],infile):
			raise Exception( 'Variable %s not found from file %s.%04u%04u.nc)!' % (infile, names[0], 0, 0) )
		variables.append(infile)
	i+=1
#
#
for infile in names:
	# Determine the true number of files
	imax,jmax=examine_fname(infile)
	if '.ts' in infile or '.ps' in infile:
		# Time or profile statistics
		reduce_ts_ps(infile,imax,jmax,variables)
	else:
		# Full 4D data
		reduce_full(infile,imax,jmax,variables)
#
# End of function combine.py