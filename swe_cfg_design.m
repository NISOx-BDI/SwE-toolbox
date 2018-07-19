function design = swe_cfg_design
% Data & design configuration file
% This builds the SwE.mat data and design structure.

% Written by Bryan Guillaume

% ---------------------------------------------------------------------
% dir Directory
% ---------------------------------------------------------------------
dir         = cfg_files;
dir.tag     = 'dir';
dir.name    = 'Directory';
dir.help    = {' '
    'Select a directory where the SwE.mat file containing the specified design matrix will be written.'};
dir.filter = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];

% ---------------------------------------------------------------------
% scans Scans
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans';
scans.help    = {' '
    'Select the images (".img" or ".nii") or a single ".mat" object containing the matrix of data with the scans in rows.'
    'The images must all have the same image dimensions, orientation, voxel size etc.'
    'The order of images or of the rows in the matrix in the ".mat" file does not matter, but the indicator vectors "Subjects", "Visits" (only for the modified SwE) and "Groups" (only for the modified SwE) must reflect this order.'};
scans.filter = {'image', 'mat'};
scans.ufilter = '.*';
scans.num     = [1 Inf];

% ---------------------------------------------------------------------
% groups Groups
% ---------------------------------------------------------------------
groups         = cfg_entry;
groups.tag     = 'groups';
groups.name    = 'Groups';
groups.help    = {' '
             'Vector of groups of length equal to the number of scans.'
             'Enter the groups in the ordering of the scans.'
             'The groups correspond to subjects sharing a common covariance matrix.'
             'The common covariance matrices are allowed to be different across groups.'
            
             ' '}';
groups.strtype = 'e';
groups.num     = [Inf 1];


% ---------------------------------------------------------------------
% subjects Subjects
% ---------------------------------------------------------------------
subjects         = cfg_entry;
subjects.tag     = 'subjects';
subjects.name    = 'Subjects';
subjects.help    = {' '
             'Vector of subjects of length equal to the number of scans.'
             'Enter the subjects in the ordering of the scans'
             ' '}';
subjects.strtype = 'e';
subjects.num     = [Inf 1];

% ---------------------------------------------------------------------
% visits Visits
% ---------------------------------------------------------------------
visits         = cfg_entry;
visits.tag     = 'visits';
visits.name    = 'Visits';
visits.help    = {' '
             'Vector of visit categories of length equal to the number of scans.'
             'Enter the visit categories in the ordering of the scans'
             'The visit categories must be consistent across subjects belonging to the same group.'
             ' '}';
visits.strtype = 'e';
visits.num     = [Inf 1];


% ---------------------------------------------------------------------
% c Vector
% ---------------------------------------------------------------------
c         = cfg_entry;
c.tag     = 'c';
c.name    = 'Vector';
c.help    = {' '
             'Vector of covariate values of length equal to the number of scans.'
             'Enter the covariate values in the ordering of the scans'
             ' ' }';
c.strtype = 'e';
c.num     = [Inf 1];
% ---------------------------------------------------------------------
% cname Name
% ---------------------------------------------------------------------
cname         = cfg_entry;
cname.tag     = 'cname';
cname.name    = 'Name';
cname.help    = {'Name of the covariate'};
cname.strtype = 's';
cname.num     = [1 Inf];
% ---------------------------------------------------------------------
% iCFI Interactions
% ---------------------------------------------------------------------
iCFI         = cfg_menu;
iCFI.tag     = 'iCFI';
iCFI.name    = 'Interactions';
iCFI.help    = {
                'For each covariate you have defined, there is an opportunity to create an additional regressor that is the interaction between the covariate and a chosen experimental factor. '
                ''
}';
iCFI.labels = {
               'None'
               'With Factor 1'
               'With Factor 2'
               'With Factor 3'
}';
iCFI.values = {1 2 3 4};
iCFI.val    = {1};
% ---------------------------------------------------------------------
% iCC Centering
% ---------------------------------------------------------------------
iCC         = cfg_menu;
iCC.tag     = 'iCC';
iCC.name    = 'Centering';
iCC.help    = {
               'The appropriate centering option is usually the one that corresponds to the interaction chosen, and ensures that main effects of the interacting factor aren''t affected by the covariate. You are advised to choose this option, unless you have other modelling considerations. '
               ''
}';
iCC.labels = {
              'Overall mean'
              'Factor 1 mean'
              'Factor 2 mean'
              'Factor 3 mean'
              'No centering'
              'User specified value'
              'As implied by ANCOVA'
              'GM'
}';
iCC.values = {1 2 3 4 5 6 7 8};
iCC.val    = {5};
% ---------------------------------------------------------------------
% cov Covariate
% ---------------------------------------------------------------------
cov         = cfg_branch;
cov.tag     = 'cov';
cov.name    = 'Covariate';
cov.val     = {c cname };
%cov.val     = {c cname iCFI iCC };
cov.help    = {'Add a new covariate to your design.'
               'Please note that no covariates is added per default. Thus, all the model covariates must be added by the user.'
               'Note also that a function called swe_splitCovariate can be used to split time-varying covariates into a cross-sectional component and a longitudinal component. This may be useful to take apart these two different mode of variation'};

% ---------------------------------------------------------------------
% generic Covariates
% ---------------------------------------------------------------------
generic         = cfg_repeat;
generic.tag     = 'generic';
generic.name    = 'Covariates';
generic.help    = {
                   'This option allows for the specification of covariates variables.'
                   ''
}';
generic.values  = {cov };
generic.num     = [1 Inf];

% ---------------------------------------------------------------------
% tm_none None
% ---------------------------------------------------------------------
tm_none         = cfg_const;
tm_none.tag     = 'tm_none';
tm_none.name    = 'None';
tm_none.val     = {1};
tm_none.help    = {'No threshold masking'};

% ---------------------------------------------------------------------
% athresh Threshold
% ---------------------------------------------------------------------
athresh         = cfg_entry;
athresh.tag     = 'athresh';
athresh.name    = 'Threshold';
athresh.help    = {
                   'Enter the absolute value of the threshold.'
                   ''
}';
athresh.strtype = 'e';
athresh.num     = [1 1];
athresh.val     = {100};
% ---------------------------------------------------------------------
% tma Absolute
% ---------------------------------------------------------------------
tma         = cfg_branch;
tma.tag     = 'tma';
tma.name    = 'Absolute';
tma.val     = {athresh };
tma.help    = {
               'Images are thresholded at a given value and only voxels at which all images exceed the threshold are included. '
               ''
               'This option allows you to specify the absolute value of the threshold.'
               ''
}';
% ---------------------------------------------------------------------
% rthresh Threshold
% ---------------------------------------------------------------------
rthresh         = cfg_entry;
rthresh.tag     = 'rthresh';
rthresh.name    = 'Threshold';
rthresh.help    = {
                   'Enter the threshold as a proportion of the global value'
                   ''
}';
rthresh.strtype = 'e';
rthresh.num     = [1 1];
rthresh.val     = {.8};
% ---------------------------------------------------------------------
% tmr Relative
% ---------------------------------------------------------------------
tmr         = cfg_branch;
tmr.tag     = 'tmr';
tmr.name    = 'Relative';
tmr.val     = {rthresh };
tmr.help    = {
               'Images are thresholded at a given value and only voxels at which all images exceed the threshold are included. '
               ''
               'This option allows you to specify the value of the threshold as a proportion of the global value. '
               ''
}';
% ---------------------------------------------------------------------
% tm Threshold masking
% ---------------------------------------------------------------------
tm         = cfg_choice;
tm.tag     = 'tm';
tm.name    = 'Threshold masking';
tm.val     = {tm_none };
tm.help    = {
              'Images are thresholded at a given value and only voxels at which all images exceed the threshold are included. '
              ''
}';
tm.values  = {tm_none tma tmr };
% ---------------------------------------------------------------------
% im Implicit Mask
% ---------------------------------------------------------------------
im         = cfg_menu;
im.tag     = 'im';
im.name    = 'Implicit Mask';
im.help    = {
              'An "implicit mask" is a mask implied by a particular voxel value. Voxels with this mask value are excluded from the analysis. '
              ''
              'For image data-types with a representation of NaN (see spm_type.m), NaN''s is the implicit mask value, (and NaN''s are always masked out). '
              ''
              'For image data-types without a representation of NaN, zero is the mask value, and the user can choose whether zero voxels should be masked out or not.'
              ''
              'By default, an implicit mask is used. '
}';
im.labels = {
             'Yes'
             'No'
}';
im.values = {1 0};
im.val    = {1};
% ---------------------------------------------------------------------
% em Explicit Mask
% ---------------------------------------------------------------------
em         = cfg_files;
em.tag     = 'em';
em.name    = 'Explicit Mask';
em.val     = {{''}};
em.help    = {
              'Explicit masks are other images containing (implicit) masks that are to be applied to the current analysis.'
              ''
              'All voxels with value NaN (for image data-types with a representation of NaN), or zero (for other data types) are excluded from the analysis. '
              ''
              'Explicit mask images can have any orientation and voxel/image size. Nearest neighbour interpolation of a mask image is used if the voxel centers of the input images do not coincide with that of the mask image.'
              ''
              'Note that, when the data is specified as a matrix save as a ".mat" file, the mask is expected to be specified as a mask vector saved as a ".mat" file.'

}';
em.filter = {'image', 'mat'};
em.ufilter = '.*';
em.num     = [0 1];

% ---------------------------------------------------------------------
% masking Masking
% ---------------------------------------------------------------------
masking         = cfg_branch;
masking.tag     = 'masking';
masking.name    = 'Masking';
masking.val     = {tm im em };
masking.help    = {
                   'The mask specifies the voxels within the image volume which are to be assessed. SPM supports three methods of masking (1) Threshold, (2) Implicit and (3) Explicit. The volume analysed is the intersection of all masks.'
                   ''
}';

% ---------------------------------------------------------------------
% ss Small sample adjustments
% ---------------------------------------------------------------------
ss         = cfg_menu;
ss.tag     = 'ss';
ss.name    = 'Small sample adjustments';
ss.labels  = { 'type 0' 'type 1' 'type 2' 'type 3' 'type C2' 'type C3' };

ss.values  = {0 1 2 3 4 5};
ss.val     = { 4 };
ss.help    = {  ' '
                'type 0: no small sample adjustment is used.'
                '             It tends to be biased and generally leads to overconfident inference in small samples.'
                'type 1: the residuals used in the SwE estimation are multiplied by sqrt(n/(n-p)).'
                '             It tends to correct for the small sample bias, but simulations seem to show that it still may lead to liberal inference in small samples.'
                'type 2: the residuals used in the SwE estimation are multiplied by 1/sqrt(1-h_ik).'
                '             It tends to correct for the small sample bias, but, even if it generally performs better than the "type 1" adjustment, simulations seems to show that it still may lead to liberal inferences in small samples.'
                'type 3: the residuals used in the SwE estimation are multiplied by 1/(1-h_ik).'
                '             It tends to correct fo the small sample bias, but simulations seem to show that it may lead to conservative inferences in small samples.'
                'type C2: the subject residuals used in the SwE estimation are multiplied by (I-H_ii)^-0.5.'
                '             Simulations seem to show that it is the best correction and removes correctly the bias in many scenarios'
                'type C3: the subject residuals used in the SwE estimation are multiplied by (I-H_ii)^-1.'
                '             Simulations seem to show that it over-correct and yield conservative inferences'
                'h_ik is the diagonal element of the hat matrix H=X''(X''X)^(-1)X corresponding to subject i and visit k.'
                'H_ii is the sub-matrix of the hat matrix H=X''(X''X)^(-1)X corresponding to subject i.'
                ' '
                }';
% ---------------------------------------------------------------------
% dof_cl Degrees of freedom type
% ---------------------------------------------------------------------
dof_cl         = cfg_menu;
dof_cl.tag     = 'dof_cl';
dof_cl.name    = 'Degrees of freedom type';
dof_cl.labels  = { 'naive' 'approx I'};
%dof_cl.labels  = { 'naive' 'approx I' 'approx II'};
dof_cl.values  = { 0 1  }; % dof_cl.values  = { 0 1 2 };
dof_cl.val     = { 0 };
dof_cl.help    = {  ' '
                'naive: naive estimation of the degrees of freedom by the total number of subjects belonging to the unseparable sub-design matrices involved in the contrast tested minus the number of non-zero pure between covariates present in these unseparable sub-design matrices.'
                '             This choice tends to overestimate the degrees of freedom, but reduce the quantity of images saved and the computation time.'
                'approx I: degrees of freedom estimation with the estimate proposed in Guillaume et al. (2014).'
                '             This choice is not recommended for the classic SwE as, with this SwE version, it generally underestimate  the degrees of freedom in small samples and a large amount of variances/covariances images (sum_i n_i*(n_i+1)/2 images) need to be saved.'
                'approx II: degrees of freedom estimation with the estimate proposed in Guillaume (2015).'
                '             This choice is not recommended for the classic SwE as, with this SwE version, it generally overestimate the degrees of freedom in small samples and a large amount of variances/covariances images (sum_i n_i*(n_i+1)/2 images) need to be saved.'
                'Note that "approx II" is not yet implemented for the classic SwE. It seems that to get accurate inferences with the classic SwE in small samples, it is preferable to consider a non-parametric analysis with the Wild Bootstrap.'
                ' '
                }';
% ---------------------------------------------------------------------
% dof_mo Small Degrees of freedom type
% ---------------------------------------------------------------------
dof_mo         = cfg_menu;
dof_mo.tag     = 'dof_mo';
dof_mo.name    = 'Degrees of freedom type';
dof_mo.labels  = { 'naive' 'approx I' 'approx II' 'approx III'};
dof_mo.values  = { 0 1 2 3};
dof_mo.val     = { 3 };
dof_mo.help    = {  ' '
                'naive: naive estimation of the degrees of freedom by the total number of subjects belonging to the unseparable sub-design matrices involved in the contrast tested minus the number of non-zero pure between covariates present in these unseparable sub-design matrices.'
                '             This choice tends to overestimate the degrees of freedom in some designs, but reduce the quantity of images saved and the computation time.'
                'approx I: degrees of freedom estimation with the estimate proposed in Guillaume et al. (2014).'
                '             This estimate assumes no missing data and does not correct for the presence of a small sample bias and a missing data bias.'
                '             Simulations seems to show that the estimate approx II and approx III are better choices (see below).' 
                'approx II: degrees of freedom estimation with an alternative estimate proposed in Guillaume (2015).'
                '             The estimate accounts partially for the presence of missing data and for a small-sample bias, but does not account for a missing data bias.'
                '             Simulations seems to indicate that it performs better than approx I, but should be used only under no missing data'
                'approx III: degrees of freedom estimation with an alternative estimate proposed in Guillaume (2015).'
                '             The estimate accounts for the presence of missing data (the missing data bias included), but not for the small-sample bias.'
                '             Simulations seems to indicate that it systematically performs better than approx II under missing data, but seems slightly less well (slightly conservative) than approx II under no missing data.'
                '             That is the recommended choice by default. Nevertheless, if there is no missing data, approx II could be selected instead.'
                ' '
               }';        
% ---------------------------------------------------------------------
% modified Modified
% ---------------------------------------------------------------------
modified         = cfg_branch;
modified.tag     = 'modified';
modified.name    = 'Modified';
modified.val     = {groups visits ss dof_mo};
modified.help    = {''};


% ---------------------------------------------------------------------
% classic Classic
% ---------------------------------------------------------------------
classic         = cfg_branch;
classic.tag     = 'classic';
classic.name    = 'Classic';
classic.val     = {ss dof_cl};
classic.help    = {''};

% ---------------------------------------------------------------------
% type SwE type
% ---------------------------------------------------------------------
type         = cfg_choice;
type.tag     = 'type';
type.name    = 'SwE type';
type.val     = {modified };
type.help    = {' '
                'Classic: classic SwE assuming that the subjects have different covariance matrices'
                'Modified: modified SwE assuming that the subjects can be divided into groups in which the subjects share a common covariance matrix'
                ' '};
type.values  = {classic modified};
% ---------------------------------------------------------------------
% g_omit Omit
% ---------------------------------------------------------------------
g_omit         = cfg_const;
g_omit.tag     = 'g_omit';
g_omit.name    = 'Omit';
g_omit.val     = {1};
g_omit.help    = {'Omit'};
% ---------------------------------------------------------------------
% global_uval Global values
% ---------------------------------------------------------------------
global_uval         = cfg_entry;
global_uval.tag     = 'global_uval';
global_uval.name    = 'Global values';
global_uval.help    = {
                       'Enter the vector of global values'
                       ''
}';
global_uval.strtype = 'e';
global_uval.num     = [Inf 1];
% ---------------------------------------------------------------------
% g_user User
% ---------------------------------------------------------------------
g_user         = cfg_branch;
g_user.tag     = 'g_user';
g_user.name    = 'User';
g_user.val     = {global_uval };
g_user.help    = {
                  'User defined  global effects (enter your own '
                  'vector of global values)'
}';
% ---------------------------------------------------------------------
% g_mean Mean
% ---------------------------------------------------------------------
g_mean         = cfg_const;
g_mean.tag     = 'g_mean';
g_mean.name    = 'Mean';
g_mean.val     = {1};
g_mean.help    = {
                  'SPM standard mean voxel value'
                  ''
                  'This defines the global mean via a two-step process. Firstly, the overall mean is computed. Voxels with values less than 1/8 of this value are then deemed extra-cranial and get masked out. The mean is then recomputed on the remaining voxels.'
                  ''
}';
% ---------------------------------------------------------------------
% globalc Global calculation
% ---------------------------------------------------------------------
globalc         = cfg_choice;
globalc.tag     = 'globalc';
globalc.name    = 'Global calculation';
globalc.val     = {g_omit };
globalc.help    = {
                   'This option is only used for PET data.'
                   ''
                   'There are three methods for estimating global effects (1) Omit (assuming no other options requiring the global value chosen) (2) User defined (enter your own vector of global values) (3) Mean: SPM standard mean voxel value (within per image full-mean/8 mask) '
                   ''
}';
globalc.values  = {g_omit g_user g_mean };
% ---------------------------------------------------------------------
% gmsca_no No
% ---------------------------------------------------------------------
gmsca_no         = cfg_const;
gmsca_no.tag     = 'gmsca_no';
gmsca_no.name    = 'No';
gmsca_no.val     = {1};
gmsca_no.help    = {'No overall grand mean scaling'};
% ---------------------------------------------------------------------
% gmscv Grand mean scaled value
% ---------------------------------------------------------------------
gmscv         = cfg_entry;
gmscv.tag     = 'gmscv';
gmscv.name    = 'Grand mean scaled value';
gmscv.help    = {
                 'The default value of 50, scales the global flow to a physiologically realistic value of 50ml/dl/min.'
                 ''
}';
gmscv.strtype = 'e';
gmscv.num     = [Inf 1];
gmscv.val     = {50};
% ---------------------------------------------------------------------
% gmsca_yes Yes
% ---------------------------------------------------------------------
gmsca_yes         = cfg_branch;
gmsca_yes.tag     = 'gmsca_yes';
gmsca_yes.name    = 'Yes';
gmsca_yes.val     = {gmscv };
gmsca_yes.help    = {
                     'Scaling of the overall grand mean simply scales all the data by a common factor such that the mean of all the global values is the value specified. For qualitative data, this puts the data into an intuitively accessible scale without altering the statistics. '
                     ''
}';
% ---------------------------------------------------------------------
% gmsca Overall grand mean scaling
% ---------------------------------------------------------------------
gmsca         = cfg_choice;
gmsca.tag     = 'gmsca';
gmsca.name    = 'Overall grand mean scaling';
gmsca.val     = {gmsca_no };
gmsca.help    = {
                 'Scaling of the overall grand mean simply scales all the data by a common factor such that the mean of all the global values is the value specified. For qualitative data, this puts the data into an intuitively accessible scale without altering the statistics. '
                 ''
                 'When proportional scaling global normalisation is used each image is separately scaled such that it''s global value is that specified (in which case the grand mean is also implicitly scaled to that value). So, to proportionally scale each image so that its global value is e.g. 20, select <Yes> then type in 20 for the grand mean scaled value.'
                 ''
                 'When using AnCova or no global normalisation, with data from different subjects or sessions, an intermediate situation may be appropriate, and you may be given the option to scale group, session or subject grand means separately. '
                 ''
}';
gmsca.values  = {gmsca_no gmsca_yes };
% ---------------------------------------------------------------------
% glonorm Normalisation
% ---------------------------------------------------------------------
glonorm         = cfg_menu;
glonorm.tag     = 'glonorm';
glonorm.name    = 'Normalisation';
glonorm.help    = {
                   'Global nuisance effects are usually accounted for either by scaling the images so that they all have the same global value (proportional scaling), or by including the global covariate as a nuisance effect in the general linear model (AnCova). Much has been written on which to use, and when. Basically, since proportional scaling also scales the variance term, it is appropriate for situations where the global measurement predominantly reflects gain or sensitivity. Where variance is constant across the range of global values, linear modelling in an AnCova approach has more flexibility, since the model is not restricted to a simple proportional regression. '
                   ''
                   '''Ancova by subject'' or ''Ancova by effect'' options are implemented using the ANCOVA options provided where each experimental factor (e.g. subject or effect), is defined. These allow e.g. different subjects to have different relationships between local and global measurements. '
                   ''
                   'Since differences between subjects may be due to gain and sensitivity effects, AnCova by subject could be combined with "grand mean scaling by subject" (an option also provided where each experimental factor is originally defined) to obtain a combination of between subject proportional scaling and within subject AnCova. '
                   ''
}';
glonorm.labels = {
                  'None'
                  'Proportional'
                  'ANCOVA'
}';
glonorm.values = {1 2 3};
glonorm.val    = {1};
% ---------------------------------------------------------------------
% globalm Global normalisation
% ---------------------------------------------------------------------
globalm         = cfg_branch;
globalm.tag     = 'globalm';
globalm.name    = 'Global normalisation';
globalm.val     = {gmsca glonorm };
globalm.help    = {
                   'This option is only used for PET data.'
                   ''
                   'Global nuisance effects are usually accounted for either by scaling the images so that they all have the same global value (proportional scaling), or by including the global covariate as a nuisance effect in the general linear model (AnCova). Much has been written on which to use, and when. Basically, since proportional scaling also scales the variance term, it is appropriate for situations where the global measurement predominantly reflects gain or sensitivity. Where variance is constant across the range of global values, linear modelling in an AnCova approach has more flexibility, since the model is not restricted to a simple proportional regression. '
                   ''
                   '''Ancova by subject'' or ''Ancova by effect'' options are implemented using the ANCOVA options provided where each experimental factor (e.g. subject or effect), is defined. These allow e.g. different subjects to have different relationships between local and global measurements. '
                   ''
                   'Since differences between subjects may be due to gain and sensitivity effects, AnCova by subject could be combined with "grand mean scaling by subject" (an option also provided where each experimental factor is originally defined) to obtain a combination of between subject proportional scaling and within subject AnCova. '
                   ''
}';
% ---------------------------------------------------------------------
% WB_no No
% ---------------------------------------------------------------------
WB_no         = cfg_const;
WB_no.tag     = 'WB_no';
WB_no.name    = 'No';
WB_no.val     = {0};
WB_no.help    = {''
  'Only a "standard" parametric SwE analysis is considered'};

% ---------------------------------------------------------------------
% WB_SwE type of SwE
% ---------------------------------------------------------------------
WB_SwE         = cfg_menu;
WB_SwE.tag     = 'WB_SwE';
WB_SwE.name    = 'Type of SwE';
WB_SwE.labels     = {'U-SwE' 'R-SwE'};
WB_SwE.values     = {0 1};
WB_SwE.val     = {0};
WB_SwE.help    = {''
  'U-SwE: unrestricted SwE which is obtained using the residuals of the unrestricted model (not imposing the null hypothesis)'
  ''
  'R-SwE: restricted SwE which is obtained using the residuals of the restricted model (imposing the null hypothesis)'
  ''
  'The R-SwE is sometimes considered in the Wild Bootstrap literature. However, in our Monte Carlo simulations (see Guillaume, 2015), no appreciable differences have been observed between this two versions when they are used with the R-WB, indicating that they could be both considered in practice. Nevertheless, it is clear that the R-SwE is generally a biased estimator of the true covariance matrix of the parameters, making it a "not-so-good" candidate for a "standard" parametric SwE analysis. The latter observation may be an argument in favour of the U-SwE, particularly for cluster analysis where a primary cluster threshold need to be defined.'
  ''
  };
% ---------------------------------------------------------------------
% WB_T_con WB contrast for T-scores
% ---------------------------------------------------------------------
WB_T_con         = cfg_entry;
WB_T_con.tag     = 'WB_T_con';
WB_T_con.name    = 'Contrast';
WB_T_con.help    = {' '
             'Contrast vector of size 1 x p to be tested, where p is the number of covariates in the unrestricted model.'
             'Enter the values in the ordering of the covariates'
             ' ' }';
WB_T_con.strtype = 'e';
WB_T_con.num     = [1 Inf];

% ---------------------------------------------------------------------
% WB_F_con WB contrast for F-scores
% ---------------------------------------------------------------------
WB_F_con         = cfg_entry;
WB_F_con.tag     = 'WB_F_con';
WB_F_con.name    = 'Contrast';
WB_F_con.help    = {' '
             'Contrast vector or matrix of size q x p to be tested, where p is the number of covariates in the unrestricted model.'
             'Enter the values in the ordering of the covariates'
             ' ' }';
WB_F_con.strtype = 'e';
WB_F_con.num     = [Inf Inf];

% ---------------------------------------------------------------------
% WB_T T-scores
% ---------------------------------------------------------------------
WB_T         = cfg_branch;
WB_T.tag     = 'WB_T';
WB_T.name    = 'T';
WB_T.val     = {WB_T_con};
WB_T.help    = {''
                     'T-scores are considered.'
                     ''
}';

% ---------------------------------------------------------------------
% WB_T T-scores
% ---------------------------------------------------------------------
WB_F         = cfg_branch;
WB_F.tag     = 'WB_F';
WB_F.name    = 'F';
WB_F.val     = {WB_F_con};
WB_F.help    = {''
                     'F-scores are considered.'
                     ''
}';
% ---------------------------------------------------------------------
% WB_Stat  Statistic type
% ---------------------------------------------------------------------
WB_stat         = cfg_choice;
WB_stat.tag     = 'WB_stat';
WB_stat.name    = 'Statistic type';
WB_stat.values     = {WB_T WB_F};
WB_stat.val     = {WB_T};
WB_stat.help    = {''
  'T: T-scores are considered. Both positive (+) and negative (-) effects p-values will be computed'
  ''
  'F: F-scores are considered.'
  ''};

% ---------------------------------------------------------------------
% WB_ss Small sample adjustments WB
% ---------------------------------------------------------------------
WB_ss         = cfg_menu;
WB_ss.tag     = 'WB_ss';
WB_ss.name    = 'Small sample adjustments for the resampling';
WB_ss.labels  = { 'type 0' 'type 1' 'type 2' 'type 3' 'type C2' 'type C3' };

WB_ss.values  = {0 1 2 3 4 5};
WB_ss.val     = { 4 };
WB_ss.help    = {  ' '
                'Small sample adjustment of the residuals to be resampled'
                ''
                'It may be different from the one used for the SwE, but it seems better to consider by default the "type C2" version for both the SwE and the WB resampling as this is the correction which seems to adjust the best the residuals to make their variances/covariances close to the one the errors'
                ''
                'Note also that, if the R-WB is considered, the toolbox accounts for the use of a restricted model to adjust the residuals, and the description below needs to be adapted accordingly (e.g., replace the unrestricted design matrix X, by the restricted design matrix X_R).'
                ''
                'Below, the description of the adjustments for an U-WB:'
                'type 0: no small sample adjustment is used.'
                '             It tends to be biased and generally leads to overconfident inference in small samples.'
                'type 1: the residuals used in the SwE estimation are multiplied by sqrt(n/(n-p)).'
                '             It tends to correct for the small sample bias, but simulations seem to show that it still may lead to liberal inference in small samples.'
                'type 2: the residuals used in the SwE estimation are multiplied by 1/sqrt(1-h_ik).'
                '             It tends to correct for the small sample bias, but, even if it generally performs better than the "type 1" adjustment, simulations seems to show that it still may lead to liberal inferences in small samples.'
                'type 3: the residuals used in the SwE estimation are multiplied by 1/(1-h_ik).'
                '             It tends to correct fo the small sample bias, but simulations seem to show that it may lead to conservative inferences in small samples.'
                'type C2: the subject residuals used in the SwE estimation are multiplied by (I-H_ii)^-0.5.'
                '             Simulations seem to show that it is the best correction and removes correctly the bias in many scenarios'
                'type C3: the subject residuals used in the SwE estimation are multiplied by (I-H_ii)^-1.'
                '             Simulations seem to show that it over-correct and yield conservative inferences'
                'h_ik is the diagonal element of the hat matrix H=X''(X''X)^(-1)X corresponding to subject i and visit k.'
                'H_ii is the sub-matrix of the hat matrix H=X''(X''X)^(-1)X corresponding to subject i.'
                ' '}';

% ---------------------------------------------------------------------
% WB_cluster_no No
% ---------------------------------------------------------------------
WB_cluster_no         = cfg_const;
WB_cluster_no.tag     = 'WB_cluster_no';
WB_cluster_no.name    = 'No';
WB_cluster_no.val     = {0};
WB_cluster_no.help    = {''
  'No cluster-wise inference will be performed'
  ''};
% ---------------------------------------------------------------------
% WB_cluster_yes Yes
% ---------------------------------------------------------------------
WB_cluster_yes         = cfg_entry;
WB_cluster_yes.tag     = 'WB_cluster_yes';
WB_cluster_yes.name    = 'Yes for image input, set the cluster-forming threshold now';
WB_cluster_yes.val     = {0.001};
WB_cluster_yes.help    = {''
                     'A cluster-wise inference will be performed alongside the voxel-wise inference. The cluster-forming threshold needs to be set now (p=0.001 per default)'
''}';
WB_cluster_yes.strtype = 'e';
WB_cluster_yes.num     = [1 1];

% ---------------------------------------------------------------------
% WB_cluster_yes_mat_clusP cluP
% ---------------------------------------------------------------------
WB_cluster_yes_mat_clusP         = cfg_entry;
WB_cluster_yes_mat_clusP.tag     = 'WB_cluster_yes_mat_clusP';
WB_cluster_yes_mat_clusP.name    = 'Set the cluster-forming threshold now for ".mat" input';
WB_cluster_yes_mat_clusP.val     = {0.001};
WB_cluster_yes_mat_clusP.help    = {''
                     'A cluster-wise inference will be performed alongside the voxel-wise inference. The cluster-forming threshold needs to be set now (p=0.001 per default)'
''}';
WB_cluster_yes_mat_clusP.strtype = 'e';
WB_cluster_yes_mat_clusP.num     = [1 1];

% ---------------------------------------------------------------------
% WB_cluster_yes_mat_type cluP
% ---------------------------------------------------------------------
WB_cluster_yes_mat_type         = cfg_menu;
WB_cluster_yes_mat_type.tag     = 'WB_cluster_yes_mat_type';
WB_cluster_yes_mat_type.name    = 'Select the type of ".mat" inputs';
WB_cluster_yes_mat_type.val     = {};
WB_cluster_yes_mat_type.labels  = { 'Volumetric (voxels)' 'Surface (vertices)'};
WB_cluster_yes_mat_type.help    = {''
                     'Select the type of ".mat" inputs. Either volumetric or surface data'
''}';
WB_cluster_yes_mat_type.values  = {0 1};

% ---------------------------------------------------------------------
% WB_cluster_yes_mat_loc loc
% ---------------------------------------------------------------------
WB_cluster_yes_mat_loc         = cfg_files;
WB_cluster_yes_mat_loc.tag     = 'WB_cluster_yes_mat_loc';
WB_cluster_yes_mat_loc.name    = 'Spatial information';
WB_cluster_yes_mat_loc.help    = {''
                     'For volumetric data, the 3D location of voxels is expected in voxel coordinates (XYZ_vox).'
                     ''
                     'For surface data, the faces (or triangles) information is expected in vertex coordinates.'
                     ''
                     'Note that this spatial information needs to be saved in a ".mat" file as a matrix of size 3 x nVoxels (or its transposed) or nFaces x 3 (or its transposed)'
}';
WB_cluster_yes_mat_loc.filter = {'mat'};
WB_cluster_yes_mat_loc.ufilter = '.*';
WB_cluster_yes_mat_loc.num     = [1 1];

% ---------------------------------------------------------------------
% WB_cluster_yesMat Yes
% ---------------------------------------------------------------------
WB_cluster_yes_mat         = cfg_branch;
WB_cluster_yes_mat.tag     = 'WB_cluster_yes_mat';
WB_cluster_yes_mat.name    = 'Yes for ".mat" input';
WB_cluster_yes_mat.val     = {WB_cluster_yes_mat_clusP WB_cluster_yes_mat_type WB_cluster_yes_mat_loc};
WB_cluster_yes_mat.help    = {''
                     'A cluster-wise inference will be performed alongside the voxel-wise inference. The cluster-forming threshold needs to be set now (p=0.001 per default) and some spatial information needs to be specified in order to form clusters.'
}';

% ---------------------------------------------------------------------
% WB_nB nB
% ---------------------------------------------------------------------
WB_nB         = cfg_entry;
WB_nB.tag     = 'WB_nB';
WB_nB.name    = 'Number of bootstraps';
WB_nB.val     = {999};
WB_nB.help    = {''
                 'This sets the number of bootstraps (nB). This will notably set the number of possible FWER-corrected p-values to nB+1 and set the minimum possible FWER-corrected p-value to 1/(nB +1).'
''}';
WB_nB.strtype = 'i';
WB_nB.num     = [1 1];


% ---------------------------------------------------------------------
% WB_cluster WB cluster-wise inference 
% ---------------------------------------------------------------------
WB_cluster         = cfg_choice;
WB_cluster.tag     = 'WB_cluster';
WB_cluster.name    = 'Cluster-wise inference';
WB_cluster.values  = {WB_cluster_no WB_cluster_yes WB_cluster_yes_mat};
WB_cluster.val     = {WB_cluster_no};
WB_cluster.help    = {''
  'No: no cluster-wise inference will be performed'
  ''
  'Yes: a cluster-wise inference will be performed alongside the voxel-wise inference. If this option is selected, the U-SwE will be automatically used in order to produce "meaningful" parametric p-value bootstrap images to be thresholded by the specified cluster-forming threshold.'
  ''};

% ---------------------------------------------------------------------
% WB_yes Yes
% ---------------------------------------------------------------------
WB_yes         = cfg_branch;
WB_yes.tag     = 'WB_yes';
WB_yes.name    = 'Yes';
WB_yes.val     = {WB_ss WB_nB WB_SwE WB_stat WB_cluster};
WB_yes.help    = {''
                     'A non-parametric Wild Bootstrap procedure is considered to analyse the data (see Guillaume, 2015)'
}';

% ---------------------------------------------------------------------
% WB Wild Bootstrap
% ---------------------------------------------------------------------
WB         = cfg_choice;
WB.tag     = 'WB';
WB.name    = 'Non-parametric Wild Bootstrap';
WB.help    = {''
              'Yes: a non-parametric Wild Bootstrap procedure is considered to analyse the data (see Guillaume, 2015)'
              ''
              'No: only a "standard" parametric SwE analysis is considered (default)'
}';
WB.values = {WB_no WB_yes};
WB.val    = {WB_no};


% ---------------------------------------------------------------------
% data Data & Design
% ---------------------------------------------------------------------
design        = cfg_exbranch;
design.tag    = 'design';
design.name   = 'Data & Design';
design.val    = {dir scans type subjects generic masking WB globalc globalm};
design.help   = {' '
                 'Module of the SwE toolbox allowing the specification of the data and design.'};
design.prog   = @swe_run_design;
design.vout   = @vout_data;

%------------------------------------------------------------------------
% Output function
%------------------------------------------------------------------------
function cdep = vout_data(job)
% Specifies the output from this modules, i.e. the filename of the mat file

cdep(1)            = cfg_dep;
cdep(1).sname      = 'SwE.mat file';
cdep(1).src_output = substruct('.','files');
cdep(1).tgt_spec   = cfg_findspec({{'filter','mat','strtype','e'}});
%------------------------------------------------------------------------
