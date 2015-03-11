#!/usr/bin/env python

# IMPORTANT: THIS SCRIPT HAS TO BE RUN FROM THE DIRECTORY WHERE IT RESIDES.
# The reason is that we use relative paths, because otherwise the filenames will be too long to
# be processed!!!  It may have something to do with the variable NAME_MAX from <limits.h>
# that we use to determine the maximum length of a filename inside 'common.h".
# If you use absolute paths, you end up with Zelda thinking that it was fed two params.ini files!


from zelda import *

# =======================
# = General information =
# =======================
# Location of Zelda executable
zelda_original_bin = "/Users/coccoinomane/Dropbox/my_projects/alcock_paczynski_working_copy/Code/zelda"
if (not(os.path.isfile(zelda_original_bin))):
    print "Zelda executable not found.  File %s does not exist." % zelda_original_bin
    sys.exit(1)

# We shall run the script using a copy of the original executable, since we may want to
# modify the original in the mean while.  (This script may require many hours to execute!)
zelda_bin = os.path.join(os.path.dirname(zelda_original_bin), "zelda_temp")
shutil.copy2(zelda_original_bin, zelda_bin) 

# Directory where result files will be stored.  
results_dir = "../../Results/results_from_script"

# Directory with galaxy catalogues
data_dir = "/Users/coccoinomane/data/ap/galaxies/"

# Set 'n_gals_per_block' to something very low (e.g. 10) to make the script very quick;
# set it to zero for a full run
n_gals_per_block = 0

# Choice of binning
binning_mode = "log"

# Return True if the number of galaxies that will survive the cuts are
# is going to very small. This is not an exact procedure, but we know
# it works for the Guo catalog.  This function is useful to determine
# whether to use the split or full catalogue when performing the
# neighbours cut.
def drastic_cut(catalog_part, cut_label, cut, n_neighbours_max=np.inf): 
    (n_neighbours_max<=5) or (cut_label=="np" and cut[0]>=50) or (cut_label=="rMag" and cut[1]<=-16) or (catalog_part=="guo2010a"  and cut_label=="stellarMass" and cut[0]>=0.028) or (catalog_part=="font2008a" and cut_label=="stellarMass" and cut[0]>=0.028*10**10)


### Choose catalogue to use.  Choose between 'guo2010a' and 'font2008a'
catalog_part = "guo2010a"  
# catalog_part = "font2008a"  


# ==============================================
# = Set variables accordingly to the catalogue =
# ==============================================

### Names of the file/s containing the catalogue.  If the 'n_grid' parameters 
### is greater than 1, then use a catalogue split into blocks.
def catalog_input_format(redshift_part, n_grid):

    if   catalog_part=="guo2010a":
        if n_grid==1:
            return os.path.join(data_dir, "galaxies_" + catalog_part + "_" + z_part + "_x_y_z_velX_velY_velZ_type_np_mvir_stellarMass_rmag_rDust.dat")
        elif n_grid>1:
            return os.path.join(data_dir, "blocks/galaxies_" + catalog_part + "_" + z_part + "_x_y_z_velX_velY_velZ_type_np_mvir_stellarMass_rmag_rDust_ngrid" + str(n_grid) + "/gal_blocks_%d_%d_%d.dat")
        else:
            print func_name() + ": 'n_grid' must be positive, not %d" & n_grid
            return "FAILURE"            
            
    elif catalog_part=="font2008a":
        if n_grid==1:
            return os.path.join(data_dir, "galaxies_" + catalog_part + "_" + z_part + "_x_y_z_vx_vy_vz_type_mhalo_stellarMass_r_SDSS.dat")
        elif n_grid>1:
            return os.path.join(data_dir, "blocks/galaxies_" + catalog_part + "_" + z_part + "_x_y_z_vx_vy_vz_type_mhalo_stellarMass_r_SDSS_ngrid" + str(n_grid) + "/gal_blocks_%d_%d_%d.dat")
        else:
            print func_name() + ": 'n_grid' must be positive, not %d" & n_grid
            return "FAILURE"




### Names & numbers of the columns (starting from 1) in the files that will be read
if   catalog_part=="guo2010a":

    # The Guo catalogue has the very useful 'np' column that stands for the number of
    # dark matter particles in the subHalo that hosts the galaxy.
    # Reference: http://gavo.mpa-garching.mpg.de/MyMillennium/Help?page=databases/guo2010a/mr

    cols = {
        "np"            :   8,
        "stellarMass"   :   10,             # Units of 10^10 Msun/h
        "rMag"          :   11              # Rest frame total absolute magnitude SDSS r-band
    }

elif catalog_part=="font2008a":

    # Reference: http://galaxy-catalogue.dur.ac.uk:8080/MyMillennium/Help?page=databases/dgalaxies/font2008a

    cols = {
        "stellarMass"   :   9,              # Units of Msun/h
        "rMag"          :   10              # Absolute rest frame SDSS r-band magnitude (Vega) of the galaxy
    }



### Redshifts to be analyzed, and relative cuts.
if catalog_part=="guo2010a": 
    
    # The cuts in stellarMass and rMag for the Guo catalogue were computed using the script 'get_equiv_cuts.py'.    
    
    runs_dict = {
        # "z0.5085"    :
        #     {
        #         "np"            :   [[0,np.inf], [50,np.inf], [100,np.inf], [500,np.inf], [1000,np.inf], [5000,np.inf]],
        #         "stellarMass"   :   [[0.02981907,np.inf], [0.14472459,np.inf], [1.75603616,np.inf], [2.9314303,np.inf], [6.1917910,np.inf]],
        #         "rMag"          :   [[-np.inf,-17.075248], [-np.inf,-18.7884140], [-np.inf,-21.283510], [-np.inf,-21.773817], [-np.inf,-22.479660]]
        #     } ,
        "z0"    :
            {
                "np"            :   [[0,np.inf], [50,np.inf], [100,np.inf], [500,np.inf], [1000,np.inf], [5000,np.inf], [0,25], [0,30], [0,35], [0,40], [0,45], [0,50], [0,100], [0,500]],
                # "stellarMass"   :   [[0.029173592,np.inf], [0.14342079,np.inf], [1.84043060,np.inf], [3.2296590,np.inf], [7.3064480,np.inf]],
                # "rMag"          :   [[-np.inf,-16.482755], [-np.inf,-18.1953770], [-np.inf,-20.702795], [-np.inf,-21.244935], [-np.inf,-22.059221]]
            } ,
        "z0.989"  :
            {
                "np"            :   [[0,np.inf], [50,np.inf], [100,np.inf], [500,np.inf], [1000,np.inf], [5000,np.inf], [0,25], [0,30], [0,35], [0,40], [0,45], [0,50], [0,100], [0,500]],
                # "stellarMass"   :   [[0.029780893,np.inf], [0.14095889,np.inf], [1.5919234,np.inf], [2.6027846,np.inf], [5.3558598,np.inf]],
                # "rMag"          :   [[-np.inf,-17.569386], [-np.inf,-19.270105], [-np.inf,-21.719086], [-np.inf,-22.17075], [-np.inf,-22.85519]]
            }
        # "z1.504"    :
        #     {
        #         "np"            :   [[0,np.inf], [50,np.inf], [100,np.inf], [500,np.inf], [1000,np.inf], [5000,np.inf]],
        #         "stellarMass"   :   [[0.02913675,np.inf], [0.13258092,np.inf], [1.38748586,np.inf], [2.2828626,np.inf], [4.7192716,np.inf]],
        #         "rMag"          :   [[-np.inf,-18.0100917], [-np.inf,-19.6758385], [-np.inf,-22.08109474], [-np.inf,-22.5187435], [-np.inf,-23.2183723]]
        #     }
    }

elif catalog_part=="font2008a":

    # The cuts for the Font catalog are just the same as the Guo ones, but without the 'np' cut and with the correct
    # units of measure for 'stellarMass'.

    runs_dict = {
        # "z0.5085"    :
        #     {
        #         "stellarMass"   :   10**10*np.array([[0.02981907,np.inf], [0.14472459,np.inf], [1.75603616,np.inf], [2.9314303,np.inf], [6.1917910,np.inf]]),
        #         "rMag"          :   [[-np.inf,-17.075248], [-np.inf,-18.7884140], [-np.inf,-21.283510], [-np.inf,-21.773817], [-np.inf,-22.479660]]
        #     } ,
        # "z0"    :
        #     {
        #         "stellarMass"   :   10**10*np.array([[0.029173592,np.inf], [0.14342079,np.inf], [1.84043060,np.inf], [3.2296590,np.inf], [7.3064480,np.inf]]),
        #         "rMag"          :   [[-np.inf,-16.482755], [-np.inf,-18.1953770], [-np.inf,-20.702795], [-np.inf,-21.244935], [-np.inf,-22.059221]]
        #     } ,
        "z0.989"  :
            {
                "stellarMass"   :   10**10*np.array([[0,np.inf], [0.029780893,np.inf], [0.14095889,np.inf], [1.5919234,np.inf], [2.6027846,np.inf], [5.3558598,np.inf]]),
                "rMag"          :   [[-np.inf,-17.569386], [-np.inf,-19.270105], [-np.inf,-21.719086], [-np.inf,-22.17075], [-np.inf,-22.85519]]
            },
        # "z1.504"    :
        #   {
        #       "stellarMass"   :   10**10*np.array([[0.02913675,np.inf], [0.13258092,np.inf], [1.38748586,np.inf], [2.2828626,np.inf], [4.7192716,np.inf]]),
        #       "rMag"          :   [[-np.inf,-18.0100917], [-np.inf,-19.6758385], [-np.inf,-22.08109474], [-np.inf,-22.5187435], [-np.inf,-23.2183723]]
        #   }
    }








# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  START OF SCRIPT  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





# ************************  TOTAL ISOLATION CRITERION, ALL CUTS, ALL REDSHIFTS  *******************************
r_iso = 4
n_neighbours_max = 1
n_separation_bins = 20
n_grid = 1

print
print "*******************************************************"
print "************** WITH ISOLATION CRITERION ***************"
print "*******************************************************"
print


for z_part,cuts_dict in runs_dict.items():

    print
    print "     ************** Now analysing redshift %s **************" % z_part
    print
    
    # Build the input filename/s
    input_format = catalog_input_format(z_part, n_grid)

    # Generate the parameters for the Zelda run.  This won't be the actual version that will be executed,
    # since the cuts are missing.  It is just a version compatible with the isolation-criterion modality.
    params = generate_params(input_format, isolation=True, n_separation_bins=n_separation_bins, r_iso=r_iso, n_neighbours_max=n_neighbours_max,
                           binning_mode=binning_mode, n_grid=n_grid, n_gals_per_block=n_gals_per_block)

    
    for cut_label,cuts in cuts_dict.items():

        print
        print "         ************** Now performing cuts on %s **************" % cut_label
        print

        # Build the output directory for the catalog & redshift
        out_dir = os.path.join(results_dir, catalog_part, z_part, cut_label + "_cut", binning_mode)

        for cut in cuts[::-1]:

            print
            print "                  > %s in (%s,%s) **************" % (cut_label, cut[0], cut[1])
            print

            # Cut parameters to be fed to Zelda
            params.update({
                "columns_to_cut"    :   cols[cut_label],
                "cuts"              :   str(cut[0]) + "," + str(cut[1])
                })
        
            # Build the part of the filename referring to the cuts & catalog
            cut_part = cut_label + str(cut[0]) + "," + str(cut[1])
            extra_naming = "_".join([catalog_part, z_part, cut_part])
    
            # Run Zelda with such cuts
            run = zelda_run(zelda_bin,
                out_dir=out_dir,
                custom_params=params,
                extra_naming=extra_naming)
            run.run()




# # ************************  ALL PAIRS, ALL CUTS, ALL REDSHIFTS  *******************************
r_max = 70
n_separation_bins = 100
n_grid = 8  # Set to one to debug the script, it will be much quicker

print
print "*******************************************************"
print "*************** CONSIDERING ALL PAIRS *****************"
print "*******************************************************"
print

for z_part,cuts_dict in runs_dict.items():

    print
    print "     ************** Now analysing redshift %s **************" % z_part
    print

    # Build the input filename/s
    input_format = catalog_input_format(z_part, n_grid)

    # Generate the parameters for the Zelda run.  This won't be the actual version that will be executed,
    # since the cuts are missing.  It is just a version compatible with the all-pairs modality.
    params = generate_params(input_format, isolation=False, n_separation_bins=n_separation_bins, r_max=r_max, n_grid=n_grid,
                           binning_mode=binning_mode, n_gals_per_block=n_gals_per_block)

    for cut_label,cuts in cuts_dict.items():

        print
        print "         ************** Now performing cuts on %s **************" % cut_label
        print

        # Build the output directory for the catalog & redshift
        out_dir = os.path.join(results_dir, catalog_part, z_part, cut_label + "_cut", binning_mode)

        for cut in cuts[::-1]:
            
            print
            print "                  > %s in (%s,%s) **************" % (cut_label, cut[0], cut[1])
            print

            # Cut parameters to be fed to Zelda
            params.update({
                "columns_to_cut"    :   cols[cut_label],
                "cuts"              :   str(cut[0]) + "," + str(cut[1])
                })

            # Build the part of the filename referring to the cuts & catalog
            cut_part = cut_label + str(cut[0]) + "," + str(cut[1])
            extra_naming = "_".join([catalog_part, z_part, cut_part])

            # Run Zelda with such cuts
            run = zelda_run(zelda_bin,
                out_dir=out_dir,
                custom_params=params,
                extra_naming=extra_naming)
            run.run()



# # ************************  TOTAL ISOLATION CRITERION, VARYING R_ISO, ALL REDSHIFTS  *******************************
# r_iso_run = [0.1, 0.3, 0.6, 1.0, 2.0, 4.0]
# z_part_run = [
#     # "z0",
#     #     "z0.5085",
#     "z0.989"# ,
#     #     "z1.504"
#     ]
# 
# sub_dir = "r_iso_cut"
# n_neighbours_max = 1
# n_separation_bins = 20
# n_grid = 1
#
# # We further subdivide each run into different values of "np"
# cut_label = "np"
# cuts = [[0,np.inf], [500,np.inf], [1000,np.inf]]
# 
# 
# print
# print "******************************************************************"
# print "************** VARYING ISOLATION CRITERION (R_ISO) ***************"
# print "******************************************************************"
# print
# 
# 
# for z_part in z_part_run:
# 
#     print
#     print "     ************** Now analysing redshift %s **************" % z_part
#     print
#     
#     # Build the input filename/s
#     input_format = catalog_input_format(z_part, n_grid)    
#     # input_format = os.path.join(data_dir, "galaxies_" + catalog_part + "_" + z_part + "_x_y_z_velX_velY_velZ_type_np_mvir_stellarMass_rmag_rDust.dat")
# 
#     for cut in cuts[::-1]:
# 
#         print
#         print "             > %s in (%s,%s) **************" % (cut_label, cut[0], cut[1])
#         print
# 
#         # Build the part of the filename referring to the cuts & catalog
#         cut_part = cut_label + str(cut[0]) + "," + str(cut[1])
#         extra_naming = "_".join([catalog_part, z_part, cut_part])
# 
#         # Build the output directory for the catalog & redshift
#         out_dir = os.path.join(results_dir, catalog_part, z_part, sub_dir, cut_part, binning_mode)
# 
#         for r_iso in r_iso_run:
# 
#             print
#             print "                     > r_isolation = %g **************" % r_iso
#             print
# 
#             # Generate the parameters for the Zelda run
#             params = generate_params(input_format, isolation=True, n_separation_bins=n_separation_bins, r_iso=r_iso, n_neighbours_max=n_neighbours_max,
#                                    binning_mode=binning_mode, n_grid=1, n_gals_per_block=n_gals_per_block)
# 
#             # Cut parameters to be fed to Zelda
#             params.update({
#                 "columns_to_cut"    :   cols[cut_label],
#                 "cuts"              :   str(cut[0]) + "," + str(cut[1])
#                 })
#     
#             # Run Zelda with such cuts
#             run = zelda_run(zelda_bin,
#                 out_dir=out_dir,
#                 custom_params=params,
#                 extra_naming=extra_naming)
#             run.run()




# ************************  ISOLATION CRITERION, VARYING N_NEIGHBOURS_MAX, ALL REDSHIFTS  *******************************
# print
# print "****************************************************************"
# print "*************** CONSIDERING VARYING NEIGHBOURS *****************"
# print "****************************************************************"
# print
# 
# 
# r_iso = 4
# n_neighbours_max_run = [1,2,5,10,20,50,100,500,1000,5000,10000]
# n_separation_bins = 20
# z_part_run = [
#     "z0"
#     # "z0.5085",
#     # "z0.989",
#     # "z1.504"
#     ]
# sub_dir = "neigh_cut"
# 
# # We further subdivide each run into different values of "np"
# cut_label="np"
# cuts = [[0,np.inf],[500,np.inf], [1000,np.inf]]
# 
# 
# for z_part in z_part_run:
# 
#     print
#     print "     ************** Now analysing redshift %s **************" % z_part
#     print
# 
#     # Build the bits that will go in the name of the output file
#     extra_naming = "_".join([catalog_part, z_part])
# 
#     for n_neighbours_max in n_neighbours_max_run:
# 
#         print
#         print "               > n_neighbours_max = %d **************" % n_neighbours_max
#         print
# 
#         for cut in cuts[::-1]:
# 
#             print
#             print "                  > %s in (%s,%s) **************" % (cut_label, cut[0], cut[1])
#             print
#                                                     
#             if (drastic_cut(catalog_part,cut_label,cut,n_neighbours_max)):
#                 # If we are cutting the number of galaxies by a great share, then we can get away
#                 # without using the split catalogue because the execution will be fast enough
#                 n_grid = 1
# 
#             else:
#                 # When you allow too many neighbours for each galaxy, or you don't cut the number of
#                 # galaxies enough, the number of pairs becomes huge and the algorithm becomes ~ N^2, so
#                 # we need to use the split catalogue
#                 n_grid = 8
# 
#             input_format = catalog_input_format(z_part, n_grid)
# 
#             # Build the part of the output filenames referring to the cuts & catalog
#             cut_part = cut_label + str(cut[0]) + "," + str(cut[1])
#             extra_naming = "_".join([catalog_part, z_part, cut_part])
# 
#             # Build the output directory for the catalog & redshift
#             out_dir = os.path.join(results_dir, catalog_part, z_part, sub_dir, cut_part, binning_mode)
#     
#             params = generate_params(input_format, isolation=True, n_separation_bins=n_separation_bins, r_iso=r_iso, n_neighbours_max=n_neighbours_max,
#                                     binning_mode=binning_mode, n_grid=n_grid, n_gals_per_block=n_gals_per_block)    
# 
#             # Cut-parameters to be fed to Zelda
#             params.update({
#                 "columns_to_cut"    :   cols[cut_label],
#                 "cuts"              :   str(cut[0]) + "," + str(cut[1])
#                 })
#     
#             run = zelda_run(zelda_bin,
#                 out_dir=out_dir,
#                 custom_params=params,
#                 extra_naming=extra_naming)
#             run.run()








# ============
# = Clean up =
# ============
os.remove(zelda_bin)



