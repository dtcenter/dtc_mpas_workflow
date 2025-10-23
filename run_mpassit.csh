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
#PBS -N mpassit
#PBS -A P48503002
#PBS -l walltime=50:00
#PBS -q regular
#PBS -o mpassit.out
#PBS -j oe 
#PBS -k oed
#PBS -l select=1:ncpus=128:mpiprocs=128:mem=209GB
#PBS -m n
#PBS -V

##SBATCH -J mpassit
##SBATCH -o logs/mpassit.%j
##SBATCH -e logs/mpassit.%j
##SBATCH -n 1200
##SBATCH --exclusive
##SBATCH --partition=hera
##SBATCH -t 02:00:00
##SBATCH -A hmtb

#
# This script runs MPASSIT
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

#setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:/scratch2/BMC/fv3lam/HWT/code/jasper/miniconda3_RL/lib

#Derecho
#ulimit -s unlimited

unlimit

module --force purge
module load ncarenv/23.09    
module load intel-oneapi/2023.2.1    
module load mkl/2023.2.0    
module load cray-mpich/8.1.27    
module load netcdf-mpi/4.9.2    
module load hdf5-mpi/1.12.2    
module load parallel-netcdf/1.12.3    
module load parallelio/2.6.2    
    
module load esmf/8.6.0    
module load cmake/3.26.3 

#module --force purge
#module load ncarenv/23.09 craype/2.7.23 intel/2023.2.1 ncarcompilers/1.0.0 cray-mpich/8.1.27
#module load cmake/3.26.3
#module load esmf/8.6.0
#module load hdf5-mpi/1.12.2
#module load mkl/2023.2.0
#module load netcdf-mpi/4.9.2
#module load parallel-netcdf/1.12.3
#module load parallelio/2.6.2

setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${PNETCDF}/lib
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${PIO}/lib

echo "$LD_LIBRARY_PATH"

module list

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
   echo "running MPASSIT on driver.csh node for member $mem"
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
   set date_file_format =  `${TOOL_DIR}/da_advance_time.exe $DATE 0 -f ccyy-mm-dd_hh.nn.ss`
   echo "Date in format is: ${date_file_format}"
   set vhr = `echo ${date_file_format}`

   # ----------------------------------
   # Make and go to working directory
   # ----------------------------------
   set rundir = $MPASSIT_OUTPUT_DIR/$start_init/ens_${mem}/${vhr}
   mkdir -p $rundir
   cd $rundir

   #-----------------------------------------------------------
   # Link necessary input files and code and fill namelist
   #-----------------------------------------------------------
   ln -sf ${MPASSIT_CODE_DIR}/bin/mpassit .
   ln -sf ${SCRIPT_DIR}/mpassit_files/* .

   setenv grid_file $MPAS_INIT_DIR/$start_init/ens_${mem}/init.nc
   setenv hist_file $EXP_DIR/$start_init/ens_${mem}/history.${date_file_format}.nc
   setenv diag_file $EXP_DIR/$start_init/ens_${mem}/diag.${date_file_format}.nc
   setenv output_file $MPASSIT_OUTPUT_DIR/$start_init/ens_${mem}/proc.${date_file_format}.nc

   $NAMELIST_TEMPLATE mpassit

   #----------------------------------------------------
   # Run MPASSIT
   #----------------------------------------------------
   rm -f ./*.log
   rm -f ./core*
   rm -f ./*.err

   #   $run_cmd ./mpassit namelist.input # Run MPASSIT!
   mpirun -np 128 ./mpassit namelist.input

   if ( $status != 0 ) then
      echo "MPASSIT failed. Exit." >> ./FAIL
      exit 6
   endif

   # Modify dx in MPASSIT output (it is wonky)
   #module load intel/2022.1.2
   #module use /contrib/METplus/modulefiles
   #module load metplus/5.1.0
   #echo "python changedx.py ${output_file}"
   #/scratch1/BMC/dtc/miniconda/miniconda3/envs/metplus_v5.1_py3.10/bin/python changedx.py ${output_file} 


   module use /glade/work/dtcrt/METplus/derecho/components/METplus/installations/modulefiles
   module load metplus/5.1.0
   echo "python changedx.py ${output_file}"
   /glade/work/dtcrt/METplus/derecho/miniconda/miniconda3/envs/metplus_v5.1_py3.10/bin/python changedx.py ${output_file} 

   # Done with this foreecast hour; go to next one
   set DATE = `$TOOL_DIR/da_advance_time.exe ${DATE} $diag_output_interval`
end # loop over time/initializations

exit 0
