************************************************************************************
! The post-processing fortran routine was obtained from :

! https://github.com/UCLALES-SALSA/UCLALES-SALSA/blob/IceDevelOrg/script/combine.f90


!These are instructions to use this routine with model outputs of the DEV branch
-----------------------------------------------------------------------------------
!

***********************************************************************************
!      To compile the fortran routine to do post-processing in PC

! You can use combine.f90 to include all variables including binned variables
! (i.e. Ncba, Dwcba, etc.) in the final nc.file

! Warning : 
! It is not recommended to include all variables if the simulation ran for a long time.
! The file size can be too large to be easily handle.
--

!To have all 4D and 5D variables run

gfortran -O2 -I/usr/include -o ./pples combine.f90 -lnetcdff

--
! You can use combined_nobinned.f90 to include just 4D variables (e.g. rc(t,x,y,z)).

! To have just 4D variables run

gfortran -O2 -I/usr/include -o ./pples_nobinned combine_nobinned.f90 -lnetcdff

--
! You can use combined_4d_single_var.f90 to include just one 4D variables (e.g. rc(t,x,y,z)).

! To have just 4D variables run

gfortran -O2 -I/usr/include -o ./pples_4d_single_var combine_4d_single_var.f90 -lnetcdff

--
! Binned variables can be processed separately using combined_binned.f90

gfortran -O2 -I/usr/include -o ./pples_binned combine_binned.f90 -lnetcdff


************************************************************************************
!                      To perform post processing in PC
! Run
./pples  <file name prefix> <key1=value1 key2=value2 ...>

! Examples
./pples SPICULE_20210605_RF04b_control_L5 deflate_level=4

--
./pples_nobinned SPICULE_20210605_RF04b_control_L5 deflate_level=4

--
./pples_4d_single_var SPICULE_20210605_RF04b_control_L5 deflate_level=4 variable_name=wwind

--
./pples_binned SPICULE_20210605_RF04b_control_L5 deflate_level=4 variable_name=Ncba

************************************************************************************
!                  To compile the fortran routine in Puhti
!**
! Post-processing UCLALES-SALSA analysis and column outputs.
!    Tomi Raatikainen, March 2020, FMI
!
! Versions
!   20200316    The first version
!
! Compile and run
! ===============
! a) Cygwin
!  Compile
!    gfortran -O2 -I/usr/include -o ./pples combine.f90 -lnetcdff
!  Run
!    ./pples <file name prefix>
! b) Puhti
!  Compile
!    module swap gcc intel-oneapi-compilers-classic
!    module load intel-oneapi-mpi/2021.6.0
!    module load netcdf-fortran/4.5.4
!    ifort -O2 -o ./pples combine.f90 -lnetcdff
!  Run
!    srun --ntasks=1 --time=00:0:10 --partition=<partition> --account=<project> pples <file name prefix> <key1=value1 key2=value2 ...>
!  Examples
!    srun --ntasks=1 --time=00:10:00 --mem=4000 --partition=fmi --account=project_2001823 pples dycoms_l4 deflate_level=4
!
!***********************************************************************************


************************************************************************************
!                    To perform post processing in Puhti
--
! You can build an sbatch file like this


-----------------------------------------------------------------------------------
#!/bin/bash
#SBATCH --job-name=pples
#SBATCH --account=project_xxxxxx
#SBATCH --partition=<partition>
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=24GB
#SBATCH --mail-type=FAIL,END
#SBATCH --mail-user=your_login@xx.xx

module swap gcc intel-oneapi-compilers-classic
module load intel-oneapi-mpi/2021.6.0
module load netcdf-fortran/4.5.4

! Remember to upload the right pples file  obtained after compilation in the
! same folder where your segmented files are

!  Run 
!    srun --ntasks=1 --time=00:0:10 --partition=<partition> --account=<project> pples <file name prefix> <key1=value1 key2=value2 ...>

! Examples
!  srun pples OMAA2019AUG12_lowcape2_phfull_M22_charge_hr11 deflate_level=4
!  srun pples_binned OMAA2019AUG12_lowcape2_phfull_M22_charge_hr11 deflate_level=4 variable_name=Ncba
!  srun pples_binned OMAA2019AUG12_lowcape2_phfull_M22_charge_hr11 deflate_level=4 variable_name=Dwcba
!  srun pples_4d_single_var OMAA2019AUG12_lowcape2_phfull_M22_charge_hr11 deflate_level=4 variable_name=wwind
------------------------------------------------------------------------------------

