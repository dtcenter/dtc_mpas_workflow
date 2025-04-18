#!/bin/csh
#########################################################################
# Script: driver.csh
#
# Purpose: Top-level driver script to run MPAS forecasts
#
# 2024: MH modifying Craig Schwartz's scripts for WPO HWT project (Hera)
#
#########################################################################
#
##BSUB -n 1
##BSUB -J MPAS_driver
##BSUB -o output.driver
##BSUB -e output.driver
##BSUB -q caldera
##BSUB -P P64000510
##BSUB -W 90
##BSUB -R "span[ptile=32]"
#
# Derecho
#PBS -S /bin/csh
#PBS -N driver
#PBS -A P48503002
#PBS -l walltime=60:00
#PBS -q main
#PBS -o ./output_file
#PBS -j oe 
#PBS -k eod 
#PBS -l select=1:ncpus=1:mpiprocs=1
#PBS -m n    
#PBS -V 
#
# Hera uses slurm
##SBATCH -J driver_expt
##SBATCH -o driver_expt.o%j
##SBATCH -N 1
##SBATCH -n 1
##SBATCH -p hera
##SBATCH -t 02:00:00 
##SBATCH -A fv3lam

# Basic system information
setenv batch_system         PBS # Type of submission system. 
setenv num_procs_per_node   12  # Number of processors per node on the machine
setenv run_cmd              "mpirun" # Command to run MPI executables #"mpirun.lsf"

# Total number of processors to use (not the number of nodes)
setenv NUM_PROCS_MPAS_ATM   12
setenv NUM_PROCS_MPAS_INIT  128 # 180
setenv num_mpas_procs_per_node  12  # Number of processors per node you want to use for MPAS forecasts (not intialization) -- not used for slurm

setenv mpas_init_walltime 120      #Run-time (minutes) for MPAS initialization
setenv mpas_fcst_walltime 240      #Run-time (minutes) for MPAS forecasts
setenv mpas_account   "P48503002"  # core-hour accoun
setenv mpas_queue     "main"   # system queue

# Decide which stages to run (run if true; lowercase):
setenv RUN_UNGRIB              false  # (true, false )
setenv RUN_MPAS_INITIALIZE     false
setenv RUN_MPAS_FORECAST       true
setenv RUN_MPASSIT             false
setenv RUN_UPP                 false

#######################################
# Directories pointing to source code #
#######################################

setenv   SCRIPT_DIR           /glade/campaign/ral/jntp/mayfield/dtc_ncar_mpas/workflow_dtc_mpas  # Location of all these .csh scripts

# Path to MPAS initialization code
setenv   MPAS_INIT_CODE_DIR   /glade/campaign/ral/jntp/mayfield/mpas_stoch/merge/MPAS-Model #the stoch_physics code compiled in MPAS 8.2.2 with intel
#setenv   MPAS_INIT_CODE_DIR   /glade/campaign/ral/jntp/mayfield/dtc_ncar_mpas/MPAS-Model #MPAS 8.2.2 with gnu

# Path to MPAS model code (could be different from initialization code)
setenv   MPAS_CODE_DIR        $MPAS_INIT_CODE_DIR

# Path to MPASSIT code
setenv   MPASSIT_CODE_DIR     /scratch2/BMC/fv3lam/HWT/code/MPASSIT

# Path to UPP code
setenv   UPP_CODE_DIR        /scratch2/BMC/fv3lam/HWT/code/UPP_NSSL
setenv   UPP_CODE_DIR        /scratch2/BMC/fv3lam/ajohns/mpas_pp

# Path to WPS; need ungrib.exe for MPAS initialization
setenv   WPS_DIR              /glade/work/wrfhelp/derecho_pre_compiled_code/wpsv4.6.0

setenv   TOOL_DIR             /glade/u/home/schwartz/utils/derecho  # Location of ./da_advance_time.exe, built from WRFDA
setenv   VTABLE_DIR           $SCRIPT_DIR # Location of where ungrib.exe VTables are
setenv   WPS_GEOG_DIR         /glade/u/home/wrfhelp/WPS_GEOG  #Directory with WPS geogrid information, needed for MPAS_INIT

############################################################
# This controls the top-level directory of the experiment  #
############################################################

setenv MESH      conus_15km_ensemble  # The MPAS mesh. Can really name whatever you want
setenv EXPT      test_nonuniform_stoch      # The experiment name that you are running

###########################
# Time/Experiment control #
###########################

setenv start_init    2022050100  # starting and ending forecast initialization times
setenv end_init      2022050200
setenv inc_init      24  #Time (hours) between forecast initializations

setenv FCST_RANGE              12    #Length of MPAS forecasts (hours)
setenv diag_output_interval    1     # Diagnostic file output frequency (hours)
setenv restart_output_interval 24000 # Restart file output frequency (hours)

setenv update_sst  .false.    # Use SST from a different dataset, other than the IC files?

setenv MPAS_REGIONAL .true.   # If .true., we are doing ungrib and mpas_init for a regional run, and running MPAS for a regional run
setenv LBC_FREQ 6          # LBC frequency for a regional run (hours)

#This is the model providing initial conditions for MPAS. Mostly needed to tell the MPAS initialization
#  how many vertical levels to expect in the GRIB files.  See run_mpas_init.csh
setenv  COLD_START_INITIAL_CONDITIONS_MODEL GEFS # (GFS, GEFS, RRFS, HRRR.pressure)
setenv  COLD_START_BOUNDARY_CONDITIONS_MODEL_CTL GFS #atj: new variable                                                                                  
setenv  COLD_START_BOUNDARY_CONDITIONS_MODEL_PERT GEFS #atj: new variable 
#setenv  COLD_START_BOUNDARY_CONDITIONS_MODEL GFS

setenv ENS_SIZE             2 # Ensemble size for the forecasts. Ensemble forecasts for all members are run all at once.
set    ie        =          1  # What ensemble member to start with?  Usually set to 1. Members ${ie} - ${ENS_SIZE} will be run

####################################################
# These control input into the sequence of scripts #
####################################################

setenv      NAMELIST_TEMPLATE        ${SCRIPT_DIR}/namelist_template.csh  #MPAS namelist tempalte. Much is "hard-wired" here, some filled on-the-fly.
setenv      STREAMS_TEMPLATE         ${SCRIPT_DIR}/streams_template.csh # Template for MPAS streams. Much is hard-wired, but some filled on-the-fly.

  #  this directory needs to be there, with the data. For the input directories, the date needs to be last (after the member).
setenv      GRIB_INPUT_DIR_MODEL      /glade/campaign/ral/jntp/mayfield/dtc_ncar_mpas/ic_bc_data# Location of global GRIB files (sub-directories by ensemble member and initialization time)
setenv      GRIB_INPUT_DIR_SST       $GRIB_INPUT_DIR_MODEL  # Where you could put GRIB files for SST

setenv      ungrib_prefx_model   "FILE"  # Probably never need to change, but there might be some need to in the future
setenv      ungrib_prefx_sst     "FILE"  #"SST"  # Usually "SST"...but can set to "FILE" (or $ungrib_prefx_model) if wanting to use SST from GFS

####################################################
# Where input files are stored
####################################################

   # example of directory "by ensemble member and date": /glade/scratch/schwartz/MPAS/ungrib_met/2023052300/ens_3
   #   these will be created


setenv   UNGRIB_OUTPUT_DIR_MODEL      /glade/campaign/ral/jntp/mayfield/dtc_ncar_mpas/expts_dtc_ncar/${EXPT}/ungrib_met  #Top-level WPS/Ungrib output directory for model data (sub-dirs by ensemble member and date)
setenv   UNGRIB_OUTPUT_DIR_SST        /glade/campaign/ral/jntp/mayfield/dtc_ncar_mpas/expts_dtc_ncar/${EXPT}/ungrib_sst  #Top-level WPS/Ungrib output directory for SST data    (sub-dirs by ensemble member and date)
setenv   MPAS_INIT_OUTPUT_DIR_TOP     /glade/campaign/ral/jntp/mayfield/dtc_ncar_mpas/expts_dtc_ncar/${EXPT}/${MESH}/mpas_init   #Top-level directory holding MPAS initialization files (sub-dirs by ens member date)
setenv   EXP_DIR_TOP                  /glade/campaign/ral/jntp/mayfield/dtc_ncar_mpas/expts_dtc_ncar/${EXPT}/${MESH}/mpas_atm  #Directory where MPAS forecasts are run
setenv   MPASSIT_OUTPUT_DIR_TOP       /glade/campaign/ral/jntp/mayfield/dtc_ncar_mpas/expts_dtc_ncar/${EXPT}/${MESH}/mpassit  #Directory where MPAS forecasts are post-processed by MPASSIT
setenv   UPP_OUTPUT_DIR_TOP           /glade/campaign/ral/jntp/mayfield/dtc_ncar_mpas/expts_dtc_ncar/${EXPT}/${MESH}/upp  #Directory where MPAS forecasts are post-processed by UPP

########################
# MPAS mesh settings
########################

# -------------------------------------------------------
# Settings for regional or global forecasts
# -------------------------------------------------------

setenv    MPAS_GRID_INFO_DIR      /glade/campaign/ral/jntp/mayfield/dtc_ncar_mpas/meshes/ # Directory containing MPAS grid files, must be there
setenv    graph_info_prefx        conus_15km.graph.info.part. #x1.${num_mpas_cells}.graph.info.part.     # Should be located in $MPAS_GRID_INFO_DIR, files must be there
setenv    grid_file_netcdf        ${MPAS_GRID_INFO_DIR}/conus_15km.grid.nc #${MPAS_GRID_INFO_DIR}/x1.${num_mpas_cells}.grid.nc  # Needs to be there before doing anything for global...not for regional though b/c you'll start with static.nc
setenv    mpas_static_data_file   ${MPAS_GRID_INFO_DIR}/conus_15km.static.nc      # Will be created during initialization, just needs to be created once.

# --------------------------------------------------------------------------------------------
# Vertical grid dimensions, same for both the ensemble and high-res determinsitic forecasts
# --------------------------------------------------------------------------------------------
setenv    num_mpas_vert_levels   55      # Number of vertical levels ( mass levels )
setenv    num_mpas_soil_levels   4       # Number of soil levels
#setenv    num_soilcat            16     #atj: added to be consistent with NSSL
#setenv    z_top_meters           25878.712      # MH commenting out since not an integer
#setenv    z_top_km               20      # MPAS model top (km)

############################################################################
# MPAS model settings; can also hardcode in $NAMELIST_TEMPLATE if you want
############################################################################
setenv time_step            60.0   # Seconds. Typically should be 4-6*dx; use closer to 4 for cycling DA
setenv radiation_frequency  30     # Minutes. Typically the same as dx (for dx = 15 km, 15 minutes)
setenv config_len_disp      3000. # Meters, diffusion length scale, which should be finest resolution in mesh (not needed in MPASv8.0+)
setenv soundings_file       dum #${SCRIPT_DIR}/sounding_locations.txt   # set to a dummy to disable soundings
setenv physics_suite        "convection_permitting" #"mesoscale_reference" #"convection_permitting_wrf390" 
#setenv deep_conv_param      "cu_ntiedtke" #"cu_grell_freitas" # "tiedtke"

# You can override physics suite individual parameterizations with these
#  Make sure you put these options in the MPAS namelist
#setenv microphysics_param        "mp_thompson"
#setenv lsm_param                 "noah"
#setenv pbl_param                 "bl_mynn"
#setenv longwave_rad_param        "rrtmg_lw"
#setenv shortwavae_rad_param      "rrtmg_sw"
#setenv sfc_layer_param           "sf_mynn"

#########################################################
#
# NOTHING BELOW HERE SHOULD NEED CHANGING (HOPEFULLY...)
#
#########################################################

set ff = ( $NAMELIST_TEMPLATE  $STREAMS_TEMPLATE )
foreach f ( $ff ) 
   if ( ! -e $f ) then
     echo "$f doesn't exist. Abort!"
     exit
   endif
end

#-------------------------

setenv DATE $start_init
while ( $DATE <= $end_init )

   echo "Processing $DATE"

   setenv FCST_RUN_DIR   ${EXP_DIR_TOP}/${DATE} # Where MPAS forecasts are run

   if ( $RUN_UNGRIB == true ) then
      foreach mem ( `seq $ie 1 $ENS_SIZE` )
         ${SCRIPT_DIR}/run_ungrib.csh $mem # send the ensemble member into run_ungrib.csh
      end
   endif



   if ( $RUN_MPAS_INITIALIZE == true ) then
      # Took out the option for creating new static file since it will exist for these experiments
      set queue_opts = "-q economy -n $NUM_PROCS_MPAS_INIT -P $mpas_account -W 30"
      set this_num_procs_per_node = $num_procs_per_node
      set this_num_needed_nodes = `echo "$NUM_PROCS_MPAS_INIT / $this_num_procs_per_node" | bc`
      if ( $this_num_needed_nodes == 0 ) then # Mostly useful for testing in serial mode to avoid lots of log files
	set this_num_needed_nodes = 1
	set this_num_procs_per_node = 1
      endif
      set last_member = $ENS_SIZE
      if ( $batch_system == LSF ) then
	 bsub -J "mpas_init_${DATE}" $queue_opts < ${SCRIPT_DIR}/run_mpas_init.csh
      else if ( $batch_system == PBS ) then
	 set pp_init = `qsub -N "mpas_init_${DATE}" -A "$mpas_account" -q "$mpas_queue" -V \
	                     -l walltime=${mpas_init_walltime}:00 -J ${ie}-${last_member} -S "/bin/csh" \
	                     -l "select=${this_num_needed_nodes}:ncpus=${this_num_procs_per_node}:mpiprocs=${this_num_procs_per_node}" \
			      ${SCRIPT_DIR}/run_mpas_init.csh`
      else if ( $batch_system == SBATCH ) then
        foreach mem ( `seq $ie 1 $ENS_SIZE` )
          sbatch ${SCRIPT_DIR}/run_mpas_init.csh $mem
        end
      else if ( $batch_system == none ) then
         ${SCRIPT_DIR}/run_mpas_init.csh 1 # really just for testing
      endif
   endif

   if ( $RUN_MPAS_FORECAST == true ) then
      if ( $batch_system == LSF ) then
	 set pp_fcst = `bsub -J "run_mpas_${DATE}[${ie}-${ENS_SIZE}]" -n $NUM_PROCS_MPAS_ATM -W $mpas_fcst_walltime -q regular < ${SCRIPT_DIR}/run_mpas.csh`
      else if ( $batch_system == PBS ) then
	 set num_needed_nodes = `echo "$NUM_PROCS_MPAS_ATM / $num_mpas_procs_per_node" | bc`
	     # Use the next line if you DON'T want a dependency condition based on completion of run_mpas_init.csh
	# set pp_fcst = `qsub -N "run_mpas_${DATE}" -A "$mpas_account" -J ${ie}-${ENS_SIZE} -W depend=afterok:${pp_init} \
	set pp_fcst = `qsub -N "run_mpas_${DATE}" -A "$mpas_account" -J ${ie}-${ENS_SIZE} \
		     -q "$mpas_queue" -V \
		     -l "select=${num_needed_nodes}:ncpus=${num_mpas_procs_per_node}:mpiprocs=${num_mpas_procs_per_node}" \
		     -l walltime=${mpas_fcst_walltime}:00 ${SCRIPT_DIR}/run_mpas.csh`
        else if ( $batch_system == SBATCH ) then
          foreach mem ( `seq $ie 1 $ENS_SIZE` )
            sbatch ${SCRIPT_DIR}/run_mpas.csh $mem 
	  end
      else if ( $batch_system == none ) then
         ${SCRIPT_DIR}/run_mpas.csh 1 # really just for testing
      endif
   endif

   if ( $RUN_MPASSIT == true ) then
      if ( $batch_system == LSF ) then
	 echo "MPASSIT functionality has not been implemented on LSF batch systems."
      else if ( $batch_system == PBS ) then
	 echo "MPASSIT functionality has not been implemented on PBS batch systems."
      else if ( $batch_system == SBATCH ) then
        foreach mem ( `seq $ie 1 $ENS_SIZE` )
          sbatch ${SCRIPT_DIR}/run_mpassit.csh $mem
        end  
      endif
   endif

   if ( $RUN_UPP == true ) then
      if ( $batch_system == LSF ) then
         echo "MPASSIT functionality has not been implemented on LSF batch systems."
      else if ( $batch_system == PBS ) then
         echo "MPASSIT functionality has not been implemented on PBS batch systems."
      else if ( $batch_system == SBATCH ) then
        foreach mem ( `seq $ie 1 $ENS_SIZE` )
          sbatch ${SCRIPT_DIR}/run_upp.csh $mem
	end
      endif
   endif

   # Done with this initialization; go to next one
   setenv DATE `$TOOL_DIR/da_advance_time.exe $DATE $inc_init`

end # loop over time/initializations

exit 0
