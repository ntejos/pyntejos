#CONFIGURATION FILE FOR RUNNING cross_match_cpairs_qsos.py, FOR
#ANALYSIS OF DIFFUSE IGM IN HIGH DENSITY ENVIRONMENTS TRACED BY GALAXY
#CLUSTERS-PAIRS. (Nicolas Tejos, 2016)

#Dataset used for QSOs
params = dict(
dataset = 'unobserved', # choose from: [unobserved, gabor, newz1, xmq, lehner, legacy, uvqs]

#Parameters for QSOs
mag_fuv_min = 19.0, # Minimum magnitude in the FUV

#Parameters for cluster-pairs
richmin  = 20, # minimum richness
max_tsep = 20, # (in Mpc) maximum transverse separation between the clusters.

#Parameters for the cross-match
max_dsep = 3,     # (in Mpc) maximum impact parameter between the cluster pair axis and the QSO sightline.
min_xsep = 3,     # (in Mpc) minimum distance between the closest cluster and the QSO sightline.
dv_group = 1000,  # (in km/s) velocity for grouping cluster-pairs as an independent system
vel_prox = 3000,  # (in km/s) velocity window for excluding zones close to the QSOs systemic redshift.

#parameters for observations
sn_min = 10     # minimum sn needed per resolution element (either for HI or OVI)
)
