


def examine_fname(infile):
    # Examine if the given file name pattern represents a valid output file
    #
    # At least one file must exists
    if not os.path.isfile('%s.%04u%04u.nc' % (infile, 0, 0)):
        print("Input files not found: ",'%s.%04u%04u.nc' % (infile, 0, 0))
        return 0,0
    #
    # Limits for the data file indices
    imax,jmax=0,0
    while imax<1000 and os.path.isfile( '%s.%04u%04u.nc' % (infile, imax, 0) ): imax+=1
    while jmax<1000 and os.path.isfile( '%s.%04u%04u.nc' % (infile, 0, jmax) ): jmax+=1
    
    return imax,jmax

def create_output_template(inputPath,inputPrefix,outputPath,outputPrefix,listOfVars,nFX,nFY,useBinned):
    
    print "Create output template with experiment prefix",inputPrefix
    
    ifile = nc.NcInput(inputPath+'/'+inputPrefix+".00000000.nc",useBinned=useBinned,varNames=listOfVars)

    if listOfVars == []: listOfVars = ifile.data.keys()

    ofile = nc.NcNew(outputPath+'/'+outputPrefix+".nc",ifile.dimensions,ifile.data,nFX,nFY)
    
    return listOfVars

def concatenateDomain(inputPath,inputPrefix,outputPath,outputPrefix,listOfVars):

    allfiles = os.listdir(inputPath)
    infile_list = [ff for ff in allfiles if inputPrefix in ff.split(".") and \
                                            "nc" in ff.split(".")        and \
                                            not 'ps' in ff.split(".")    and \
                                            not 'ts' in ff.split(".")]    
    
    print"Beging sequence for subdomain concatenation"
    print ""
    print"Input files:"
    print infile_list
    print""

    ntot = len(infile_list)
    i = 1
    
    print"Processing",ntot,"files and",len(listOfVars),"variables"

    ofile = nc.NcExtend(outputPath+'/'+outputPrefix+".nc")
    for infile in infile_list:
        print "Processing input from",infile, str(i)+"/"+str(ntot)
        print ""
        ifile = nc.NcInput(inputPath+'/'+infile,varNames=listOfVars)
        ofile.processData(ifile.dimensions,ifile.data)
        ifile.NcClose()
        i+=1
        del ifile
        gc.collect()
        print "Done", infile
        print "======================================="
        print ""
    ofile.NcClose()


def parse_args(arguments):

    ifile = ''
    ipath = '../bin/'
    ofile = ''
    opath = '../bin/'
    lvars = []
    useBinned = True

    try: 
        opts,args = getopt.getopt(arguments,'i:v:o:',["input-path=","output-path=","no-binned"])
    except getopt.GetoptError:
        print("Error options parser")
        sys.exit(2)

    check=False
    for opt,arg in opts:
        if opt == '-i':
            check = True
    if not check:
        print "Config error: Must define at least the input prefix -i"
        sys.exit(2)

    for opt,arg in opts:
        if opt == '-i':
            # Input file prefix
            ifile = arg
            # Use same prefix for outputfile. Overwritten if defined by -o
            ofile = arg+"_postp"
        elif opt == '-o':
            ofile = arg
        elif opt == '-v':
            lvars = arg.split(',')
        elif opt == '--input-path':
            # If other than UCLALES-SALSA/bin
            ipath = arg
        elif opt == '--output-path':
            # If other than UCLALES-SALSA/bin
            opath = arg
        elif opt == "--no-binned":
            # Don't process binned (==5d) variables 
            useBinned=False

    print"CONFIG:"
    print"--------------"
    print"Input path:",ipath,", with file prefix:",ifile
    print"Processed output path:",opath,", with file prefix",ofile
    print"--------------"
    print""

    return ifile,ipath,ofile,opath,lvars,useBinned


import numpy as np
import sys, getopt, os, gc
import classNCFile as nc

argu = sys.argv[1:]

# Parse arguments
inputPrefix,inputPath,outputPrefix,outputPath,listOfVars,useBinned = parse_args(argu)

# Number of files in x and y
nFX,nFY = examine_fname(inputPath+'/'+inputPrefix)

listOfVars = create_output_template(inputPath,inputPrefix,outputPath,outputPrefix,listOfVars,nFX,nFY,useBinned)

concatenateDomain(inputPath,inputPrefix,outputPath,outputPrefix,listOfVars)
