#!/bin/csh
#
##BSUB -n 256
##BSUB -J MPAS_fcst
##BSUB -o output.fc
##BSUB -e output.fc
##BSUB -q regular
##BSUB -P NMMM0032
##BSUB -W 60
##BSUB -R "span[ptile=16]"

#PBS -S /bin/csh
#PBS -N upp
#PBS -A P48503002
#PBS -l walltime=90:00
#PBS -q regular
#PBS -o upp.out
#PBS -j oe 
#PBS -k oed
#PBS -l select=1:ncpus=1:mpiprocs=1
#PBS -m n
#PBS -V

##SBATCH -J upp
##SBATCH -o logs/upp.%j
##SBATCH -e logs/upp.%j
##SBATCH -n 1
##SBATCH -N 1
##SBATCH --exclusive
##SBATCH --partition=hera
##SBATCH -t 02:30:00
##SBATCH -A hmtb

#
# This script runs UPP
#
# MH Load modules and set envars for MPAS on Hera (updated for Rocky 8):
#module purge
#
#module load cmake/3.28.1
#module load gnu
#module load intel/2023.2.0
#module load impi/2023.2.0
#module load pnetcdf/1.12.3
#module load szip
#module load hdf5parallel/1.10.5
#module load netcdf-hdf5parallel/4.7.0
#setenv PNETCDF /apps/pnetcdf/1.12.3/intel_2023.2.0-impi
#setenv CMAKE_C_COMPILER mpiicc
#setenv CMAKE_CXX_COMPILER mpiicpc
#setenv CMAKE_Fortran_COMPILER mpiifort
#
#setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/scratch2/BMC/fv3lam/HWT/code/jasper/miniconda3_RL/lib
#
#export PNETCDF=/apps/pnetcdf/1.12.3/intel_2023.2.0-impi
#export CMAKE_C_COMPILER=mpiicc
#export CMAKE_CXX_COMPILER=mpiicpc
#export CMAKE_Fortran_COMPILER=mpiifort
#export WRF_DIR=/scratch2/BMC/fv3lam/ajohns/WRFV3
#export JASPERLIB=/lib64
#export JASPERINC=/usr/include/jasper
#export I_MPI_CC=icc
#export I_MPI_CXX=icpc
#export I_MPI_F77=ifort
#export I_MPI_F90=ifort

#module use /glade/work/wmayfield/dtc_ncar_mpas/ufs-srweather-app/modulefiles
#module load build_derecho_intel

module --force purge
module use /glade/work/wmayfield/dtc_ncar_mpas/UPP/modulefiles
module load derecho_intel 

# -----------------------------------------------------
# Deal with batch stuff and specifying ensemble member
# -----------------------------------------------------
if ( $batch_system == LSF ) then
   set mem = $LSB_JOBINDEX
else if ( $batch_system == PBS ) then
   set mem = $PBS_ARRAY_INDEX #PBS with "-J" flag (PBS pro)
else if ( $batch_system == none ) then
   set mem = $1
else if ( $batch_system == SBATCH ) then
   set mem = $1
   echo "running UPP on driver.csh node for member $mem"
else
   echo "what is your batch sytem? exit"
   exit 2
endif

set start_init = $DATE   # From driver
set DATE = $DATE   # From driver
set diag_output_interval = $diag_output_interval   # From driver
set end_time = `$TOOL_DIR/da_advance_time.exe $DATE $FCST_RANGE`

while ( $DATE <= $end_time)

   # -------------------------------------
   # Get current date into proper format
   # -------------------------------------
   setenv date_file_format `${TOOL_DIR}/da_advance_time.exe $DATE 0 -f ccyy-mm-dd_hh.nn.ss`
   setenv date_file_format_colon `${TOOL_DIR}/da_advance_time.exe $DATE 0 -w`
   echo "Date in format is: ${date_file_format}"
   set vhr = `echo ${date_file_format}`

   # ----------------------------------
   # Make and go to working directory
   # ----------------------------------
   set rundir = ${UPP_OUTPUT_DIR}/$start_init/ens_${mem}/${vhr}
   mkdir -p $rundir
   cd $rundir

   #-----------------------------------------------------------
   # Link necessary input files and code and fill namelist
   #-----------------------------------------------------------
   #ln -sf ${UPP_CODE_DIR}/bin/unipost.exe .
   ln -sf ${UPP_CODE_DIR}/upp.x .
   ln -sf ${SCRIPT_DIR}/upp_files/params_grib2_tbl_new params_grib2_tbl_new
   ln -sf ${SCRIPT_DIR}/upp_files/postxconfig-NT.txt postxconfig-NT.txt
#   ln -sf ${SCRIPT_DIR}/upp_files/ETAMPNEW_DATA nam_micro_lookup.dat
#   ln -sf ${SCRIPT_DIR}/upp_files/ETAMPNEW_DATA.expanded_rain hires_micro_lookup.dat
   ln -sf ${SCRIPT_DIR}/upp_files/*.dat .


   ln -sf ${UPP_CODE_DIR}/../crtm/fix/SpcCoeff/Big_Endian/*Coeff.bin .
   ln -sf ${UPP_CODE_DIR}/../crtm/fix/TauCoeff/ODPS/Big_Endian/*Coeff.bin .
   ln -sf ${UPP_CODE_DIR}/../crtm/fix/AerosolCoeff/Big_Endian/*Coeff.bin .
   ln -sf ${UPP_CODE_DIR}/../crtm/fix/CloudCoeff/Big_Endian/*Coeff.bin .
   ln -sf ${UPP_CODE_DIR}/../crtm/fix/EmisCoeff/IR_Ice/SEcategory/Big_Endian/*Coeff.bin .
   ln -sf ${UPP_CODE_DIR}/../crtm/fix/EmisCoeff/IR_Land/SEcategory/Big_Endian/*Coeff.bin .
   ln -sf ${UPP_CODE_DIR}/../crtm/fix/EmisCoeff/IR_Snow/SEcategory/Big_Endian/*Coeff.bin .
   ln -sf ${UPP_CODE_DIR}/../crtm/fix/EmisCoeff/IR_Water/Big_Endian/*Coeff.bin .
   ln -sf ${UPP_CODE_DIR}/../crtm/fix/EmisCoeff/MW_Water/Big_Endian/*Coeff.bin .

   setenv mpassit_file $MPASSIT_OUTPUT_DIR/$start_init/ens_${mem}/proc.${date_file_format}.nc

   $NAMELIST_TEMPLATE upp

   #----------------------------------------------------
   # Run MPASSIT
   #----------------------------------------------------
   rm -f ./*.log
   rm -f ./core*
   rm -f ./*.err

   #   $run_cmd ./unipost.exe  # Run UPP!
   #mpirun -np 1 ./unipost.exe > run.ouput  # Run UPP! # MH note - this run command works but should be updated to not be hardcoded.
#   mpiexec -n 2 ./upp.x > run.ouput
   ./upp.x > run.ouput

   set fhr = `ls WRFPRS* | cut -c12-13`
   mv WRFPRS.GrbF${fhr}* ../mpas.t00z.prslev.f0${fhr}.conus_3km.grib2

   if ( $status != 0 ) then
      echo "UPP failed. Exit." >> ./FAIL
      exit 6
   endif

   # Done with this foreecast hour; go to next one
   echo "The current date is ${DATE}"
   set DATE = `$TOOL_DIR/da_advance_time.exe ${DATE} $diag_output_interval`
   echo "The new date is ${DATE}"
end # loop over time/initializations

exit 0
