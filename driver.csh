#!/bin/csh
#########################################################################
# Script: driver.csh
#
# Purpose: Top-level driver script to run MPAS-A forecasts
# DTC modified scripts that developed by Craigs Schwartz (NCAR/MMM)
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

# Decide which stages to run (run if true, otherwise, false; lowercase to be used):
setenv RUN_UNGRIB              false  # (true, false )
setenv RUN_MPAS_INITIALIZE     false
setenv RUN_MPAS_FORECAST       true
setenv RUN_MPASSIT             false
setenv RUN_UPP                 false

#######################################
# Directories pointing to source code and required datasets #
#######################################
# Path and directory containing all code and scripts
setenv   HOMEDIR              /glade/campaign/ral/jntp/weiweili/MMM_model/workflow/dev_mwf

# Path to all scripts (e.g., *.csh) under dtc_mpas_workflow
setenv   SCRIPT_DIR           ${HOMEDIR}/dtc_mpas_workflow

# Path to MPAS model code - could be different from MPAS_INIT_CODE_DIR (MPAS initialization code; fetch the code: git clone https://github.com/MPAS-Dev/MPAS-Model.git)
#setenv   MPAS_CODE_DIR        ${HOMEDIR}/MPAS-Model
setenv    MPAS_CODE_DIR       /glade/campaign/ral/jntp/mayfield/mpas_stoch/merge/MPAS-Model_spptint

# Path to MPAS-A initialization source code 
setenv   MPAS_INIT_CODE_DIR   $MPAS_CODE_DIR

# Path to MPASSIT code
# Optional (one need to have it built on Derecho)
setenv   MPASSIT_CODE_DIR     /scratch2/BMC/fv3lam/HWT/code/MPASSIT

# Path to UPP code
# Optional (one needs to have it built on Derecho)
setenv   UPP_CODE_DIR        /scratch2/BMC/fv3lam/HWT/code/UPP_NSSL
setenv   UPP_CODE_DIR        /scratch2/BMC/fv3lam/ajohns/mpas_pp

# Path to WPS, where ungrib.exe is located for MPAS initialization
# Not platform agnostic. Hera has it built somewhere, which is used in mpas_app.
setenv   WPS_DIR              /glade/work/wrfhelp/derecho_pre_compiled_code/wpsv4.6.0

# Path to da_advance_time.exe built from WRFDA
# So far, only point to Craig's personal dir. 
setenv   TOOL_DIR             /glade/u/home/schwartz/utils/derecho

# Path to VTables, required by ungrib.exe ($WPS_DIR); available in WPS package
setenv   VTABLE_DIR           $SCRIPT_DIR

# Path to MPAS static geographic datasets, needed for MPAS initialization
# TODO: Not platform agnostic. Can pull from https://www2.mmm.ucar.edu/projects/mpas/site/access_code/static.html and created a directory to contain these files
# TODO: not all files in WPS_GEOG are used for MPAS. Also consider renaming WPS to something more MPAS?
setenv   WPS_GEOG_DIR         /glade/u/home/wrfhelp/WPS_GEOG

############################################################
# Naming top-level directory of the experiment #
############################################################

# Self-defined experiment name
setenv EXPT      test_wl      # The experiment name that you are running

# Self-defined experiment name with mesh information (a subdirectory of $EXPT)
setenv MESH      conus_15km_ensemble

###########################
# Experiment configurations #
###########################

# Starting and Ending forecast initialization date/time (forecast cycle related)
setenv START_INIT    2022050100  # starting and ending forecast initialization times
setenv END_INIT      2022050100

# Interval (in hours) between forecast cycles/initializations, needed for running WRF-DA (da_advance_time.exe)
setenv INC_INIT      24

# Length of forecasts/simulations (in hours)
setenv FCST_RANGE              12    #Length of MPAS forecasts (hours)

# Diagnostic file output frequency (in hours)
setenv diag_output_interval    1

# Restart file output frequency (in hours)
setenv restart_output_interval 24000 # Restart file output frequency (hours)

# If SST and sea-ice fields to be periodically updated (from an external source) as the model runs (matters for longer simulations)
setenv update_sst  .false.    

# If run simulations over regional domains (aka limited-area simulation) - how mpas_init and ungrib work would be different
setenv MPAS_REGIONAL .true.

# Frequency of applying lateral boundary condition (LBC; in hours)
setenv LBC_FREQ 6          

#This is the model providing initial conditions for MPAS. Mostly needed to tell the MPAS initialization
#  how many vertical levels to expect in the GRIB files.  See run_mpas_init.csh
setenv  COLD_START_INITIAL_CONDITIONS_MODEL GEFS # (GFS, GEFS, RRFS, HRRR.pressure)
setenv  COLD_START_BOUNDARY_CONDITIONS_MODEL_CTL GFS #atj: new variable                                                                                  
setenv  COLD_START_BOUNDARY_CONDITIONS_MODEL_PERT GEFS #atj: new variable 
#setenv  COLD_START_BOUNDARY_CONDITIONS_MODEL GFS

# Ensemble size for the forecasts (note: all ensemble members are run all at once)
# ENS_SIZE = 1 to run deterministic forecast
setenv ENS_SIZE             1 # Ensemble size for the forecasts. Ensemble forecasts for all members are run all at once.

# Which index the ensemble member should start with (must >= 1 so that ${IENS}-${ENS_SIZE}>0) 
setenv IENS                 1  

########################################################################
# Parameters for preparing IC (and LBCs if running LAM) using ungrib.exe
#######################################################################

# MPAS namelist and streams tempalte (HARD-WIRED; FILLED ON-THE-FLY; DO NOT CHANGE)
setenv      NAMELIST_TEMPLATE        ${SCRIPT_DIR}/namelist_template.csh
setenv      STREAMS_TEMPLATE         ${SCRIPT_DIR}/streams_template.csh 

# Directory to IC and LBC data (must be in GRIB format), which will be ingested into ungrib.exe to generate intermediate IC and LBC files (sub-directories by ensemble member and initialization time)
# Need to manually set up ens_* and lbc_*0
setenv      GRIB_INPUT_DIR_MODEL      /glade/campaign/ral/jntp/mayfield/dtc_ncar_mpas/ic_bc_data

# As GRIB_INPUT_DIR_MODEL but for SST data
setenv      GRIB_INPUT_DIR_SST       $GRIB_INPUT_DIR_MODEL

# Prefix of intermediate IC and LBC file for running ungrib.exe (DO NOT CHANGE))
setenv      ungrib_prefx_model   "FILE"

# Prefix of intermediate SST file for running ungrib.exe (Usually "SST" but can set to "FILE" (or $ungrib_prefx_model) if using GFS SST) 
setenv      ungrib_prefx_sst     "FILE"

# Directory to contain intermediate files "FILE*" generated by WPS ungrib.exe for IC and LBC data (sub-dirs arranged by date/ensemble_member)
setenv   UNGRIB_OUTPUT_DIR_MODEL      ${HOMEDIR}/${EXPT}/ungrib_met
# Directory to contain intermediate files "FILE*" generated by WPS ungrib.exe for SST data (if SST is periodically updated) (sub-dirs arranged by date/ensemble_member)
setenv   UNGRIB_OUTPUT_DIR_SST        ${HOMEDIR}/${EXPT}/ungrib_sst


########################################################################
# Directories containing files/data ingested to and produced by MPAS-A
########################################################################
# Directory to hold MPAS initialization files (sub-dirs arranged by date/ensemble_member)
setenv   MPAS_INIT_DIR                ${HOMEDIR}/${EXPT}/${MESH}/mpas_init

# Directory to save MPAS-A simulation output (sub-dirs arranged by date/ensemble_member)
setenv   EXP_DIR                      ${HOMEDIR}/${EXPT}/${MESH}/mpas_atm

#########################################################
# Directories containing postprocessed data if appliable
#########################################################
# Directory to save all post-processed data produced by MPASSIT
setenv   MPASSIT_OUTPUT_DIR           ${HOMEDIR}/${EXPT}/${MESH}/mpassit
# Directory to save all post-processed data produced by UPP
setenv   UPP_OUTPUT_DIR_TOP           ${HOMEDIR}/${EXPT}/${MESH}/upp


############################
# MPAS mesh and static files
############################

# -------------------------------------------------------
# Settings for regional or global forecasts
# -------------------------------------------------------

# Directory containing MPAS mesh, grid and static files
# TODO: Archive all the DTC generated mesh to a generic place under JNT?
setenv    MPAS_GRID_INFO_DIR      /glade/campaign/ral/jntp/mayfield/dtc_ncar_mpas/meshes/ # Directory containing MPAS grid files, must be there

# Specify which "mesh decomposition file" (under $MPAS_GRID_INFO_DIR) to be used, usually named as *graph.info.part.* (in ASCII format), the number after this prefixdenotes an appropriate number of partitions that are equal to the number of MPI tasks that will be used
setenv    graph_info_prefx        conus_15km.graph.info.part.

# Path to SCVT mesh (under $MPAS_GRID_INFO_DIR) - usually named as *.grid.nc. (in netCDF format)
# !!CRITICAL: For running global simulation, this must exist or be pre-generated when using this workflow
# For global run, available meshes can be downloaded at https://mpas-dev.github.io/atmosphere/atmosphere_meshes.html (for NCAR HPC users, some can be found at /glade/campaign/mmm/wmr/mpas_tutorial/meshes/). 
# For running CONUS simulations, one can find a few existing meshes under /glade/campaign/ral/jntp/mayfield/dtc_ncar_mpas/meshes. These meshes can be created using tool "create_region" (https://github.com/MPAS-Dev/MPAS-Limited-Area) that uses MPAS global grid to produce a regional area grid given a region specifications
# For running other limited-area simulation, one can generate regional mesh by using "create_region" tool.
# TODO: need to have a generic place to store these pre-created meshes. Also WL don't understand "#not for regional though b/c you'll start with static.nc" 
setenv    grid_file_netcdf        ${MPAS_GRID_INFO_DIR}/conus_15km.grid.nc 

# Path to MPAS static file to be created if not exist (NOTE: it is required to initialize MPAS and is generated based on $WPS_GEOG_DIR) (under $MPAS_GRID_INFO_DIR) - usually named as *.static.nc. (in netCDF format)
# For global run, available static files can be downloaded at https://mpas-dev.github.io/atmosphere/atmosphere_meshes.html (for NCAR HPC users, some can be found at /glade/campaign/mmm/wmr/mpas_tutorial/meshes/). 
# TODO: Not clear how regional static file is genereate, which approach? 1) use the create_region tool to create a subset of an existing global "static" file for specifield region, as in the tutorial, or 2) apply init_atmosphere to regional grid $grid_file_netcdf  
setenv    mpas_static_data_file   ${MPAS_GRID_INFO_DIR}/conus_15km.static.nc


############################################################################
# Other MPAS model configurations (dimensions, physics namelist options)
# Note: These can also be hardcoded in $NAMELIST_TEMPLATE
############################################################################
# --------------------------------------------------------------------------------------------
# Vertical grid dimensions, same for both the ensemble and high-res determinsitic forecasts
# --------------------------------------------------------------------------------------------
setenv    num_mpas_vert_levels   55      # Number of vertical levels ( mass levels )
setenv    num_mpas_soil_levels   4       # Number of soil levels
#setenv    num_soilcat            16     #atj: added to be consistent with NSSL
#setenv    z_top_meters           25878.712      # MH commenting out since not an integer
#setenv    z_top_km               20      # MPAS model top (km)

setenv    time_step              60.0   # Seconds. Typically should be 4-6*dx; use closer to 4 for cycling DA
setenv    radiation_frequency    30     # Minutes. Typically the same as dx (for dx = 15 km, 15 minutes)
setenv    config_len_disp        3000. # Meters, diffusion length scale, which should be finest resolution in mesh (not needed in MPASv8.0+)
setenv    soundings_file         dum #${SCRIPT_DIR}/sounding_locations.txt   # set to a dummy to disable soundings
setenv    physics_suite          "convection_permitting" #"mesoscale_reference" #"convection_permitting_wrf390" 
#setenv    deep_conv_param       "cu_ntiedtke" #"cu_grell_freitas" # "tiedtke"

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

setenv DATE $START_INIT
while ( $DATE <= $END_INIT )

   echo "Processing $DATE"

   setenv FCST_RUN_DIR   ${EXP_DIR}/${DATE} # Where MPAS forecasts are run

   if ( $RUN_UNGRIB == true ) then
      foreach mem ( `seq $IENS 1 $ENS_SIZE` )
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
         set ii = `expr $IENS + 1`
         set last_member = "${ii}:2"
	 set pp_init = `qsub -N "mpas_init_${DATE}" -A "$mpas_account" -q "$mpas_queue" -V \
	                     -l walltime=${mpas_init_walltime}:00 -J ${IENS}-${last_member} -S "/bin/csh" \
	                     -l "select=${this_num_needed_nodes}:ncpus=${this_num_procs_per_node}:mpiprocs=${this_num_procs_per_node}" \
			      ${SCRIPT_DIR}/run_mpas_init.csh`
      else if ( $batch_system == SBATCH ) then
        foreach mem ( `seq $IENS 1 $ENS_SIZE` )
          sbatch ${SCRIPT_DIR}/run_mpas_init.csh $mem
        end
      else if ( $batch_system == none ) then
         ${SCRIPT_DIR}/run_mpas_init.csh 1 # really just for testing
      endif
   endif

   if ( $RUN_MPAS_FORECAST == true ) then
      if ( $batch_system == LSF ) then
	 set pp_fcst = `bsub -J "run_mpas_${DATE}[${IENS}-${ENS_SIZE}]" -n $NUM_PROCS_MPAS_ATM -W $mpas_fcst_walltime -q regular < ${SCRIPT_DIR}/run_mpas.csh`
      else if ( $batch_system == PBS ) then
	 set num_needed_nodes = `echo "$NUM_PROCS_MPAS_ATM / $num_mpas_procs_per_node" | bc`
	     # Use the next line if you DON'T want a dependency condition based on completion of run_mpas_init.csh
	# set pp_fcst = `qsub -N "run_mpas_${DATE}" -A "$mpas_account" -J ${IENS}-${ENS_SIZE} -W depend=afterok:${pp_init} \
#   set pp_fcst = `qsub -N "run_mpas_${DATE}" -A "$mpas_account" -J ${IENS}-${ENS_SIZE} \
         set ii = `expr $IENS + 1`
         set last_member = "${ii}:2"
         set pp_fcst = `qsub -N "run_mpas_${DATE}" -A "$mpas_account" -J ${IENS}-${last_member} \
		     -q "$mpas_queue" -V \
		     -l "select=${num_needed_nodes}:ncpus=${num_mpas_procs_per_node}:mpiprocs=${num_mpas_procs_per_node}" \
		     -l walltime=${mpas_fcst_walltime}:00 ${SCRIPT_DIR}/run_mpas.csh`
        else if ( $batch_system == SBATCH ) then
          foreach mem ( `seq $IENS 1 $ENS_SIZE` )
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
        foreach mem ( `seq $IENS 1 $ENS_SIZE` )
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
        foreach mem ( `seq $IENS 1 $ENS_SIZE` )
          sbatch ${SCRIPT_DIR}/run_upp.csh $mem
	end
      endif
   endif

   # Done with this initialization; go to next one
   setenv DATE `$TOOL_DIR/da_advance_time.exe $DATE $INC_INIT`

end # loop over time/initializations

exit 0
