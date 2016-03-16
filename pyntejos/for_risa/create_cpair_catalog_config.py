#CONFIGURATION FILE FOR RUNNING create_cpair_catalog.py, FOR ANALYSIS
#OF DIFFUSE IGM IN HIGH DENSITY ENVIRONMENTS TRACED BY GALAXY
#CLUSTERS-PAIRS. (Nicolas Tejos, 2016)

param = dict(
# Which catalog
dataset = 'des_yr1', # nick name of dataset
filename = './data/redmapper_des_yr1.fits',  # location of the redmapper catalog

#Parameters for clusters
zmin    = 0.1,  # minimum redshift of the experiment
zmax    = 0.5,  # maximum redshift of the experiment
richmin = 20,   # minimum richness for clusters

#Parameters for cluster-pairs
max_tsep   = 20,     # (in Mpc) maximum transverse separation between the clusters.
max_dv     = 2000,   # (in km/s) maximum velocity difference allowed for cluster members.
z_err_max1 = 0.0002, # maximum redshift uncertainty for at least one of the clusters.
z_err_max2 = 0.02,   # maximum redshift uncertainty for both members of the cluster-pair.

#Parameters for grouping cluster-pairs if they are too close in redshift
group_dv = 1000     # (in km/s) velocity for grouping cluster-pairs as a single structure
)