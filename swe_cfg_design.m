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
dir.help    = {'Select a directory where the SwE.mat file containing the specified design matrix will be written.'};
dir.filter = 'dir';
dir.ufilter = '.*';
dir.num     = [1 1];

% ---------------------------------------------------------------------
% scans Scans
% ---------------------------------------------------------------------
scans         = cfg_files;
scans.tag     = 'scans';
scans.name    = 'Scans';
scans.help    = {'Select the images.  They must all have the same image dimensions, orientation, voxel size etc.'};
scans.filter = 'image';
scans.ufilter = '.*';
scans.num     = [1 Inf];

% ---------------------------------------------------------------------
% groups Groups
% ---------------------------------------------------------------------
groups         = cfg_entry;
groups.tag     = 'groups';
groups.name    = 'Groups';
groups.help    = {
             'Vector of groups.'
             'Enter the groups vector in the ordering of the scans'
             'The groups are used to define subjects sharing a common covariance matrix'
}';
groups.strtype = 'e';
groups.num     = [Inf 1];


% ---------------------------------------------------------------------
% subjects Subjects
% ---------------------------------------------------------------------
subjects         = cfg_entry;
subjects.tag     = 'subjects';
subjects.name    = 'Subjects';
subjects.help    = {
             'Vector of subjects.'
             'Enter the subjects vector in the ordering of the scans'
}';
subjects.strtype = 'e';
subjects.num     = [Inf 1];

% ---------------------------------------------------------------------
% visits Visits
% ---------------------------------------------------------------------
visits         = cfg_entry;
visits.tag     = 'visits';
visits.name    = 'Visits';
visits.help    = {
             'Vector of visits.'
             'Enter the visits vector in the ordering of the scans'
             'The visits has to be consistent accross subjects'
}';
visits.strtype = 'e';
visits.num     = [Inf 1];


% ---------------------------------------------------------------------
% c Vector
% ---------------------------------------------------------------------
c         = cfg_entry;
c.tag     = 'c';
c.name    = 'Vector';
c.help    = {
             'Vector of covariate values.'
             'Enter the covariate values in the ordering of the scans'
}';
c.strtype = 'e';
c.num     = [Inf 1];
% ---------------------------------------------------------------------
% cname Name
% ---------------------------------------------------------------------
cname         = cfg_entry;
cname.tag     = 'cname';
cname.name    = 'Name';
cname.help    = {'Name of covariate'};
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
cov.help    = {'Add a new covariate to your design'};

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
generic.num     = [0 Inf];

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
              ''
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
}';
em.filter = 'image';
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
% ss Small samples adjustments
% ---------------------------------------------------------------------
ss         = cfg_menu;
ss.tag     = 'ss';
ss.name    = 'Small samples adjustments';
ss.help    = {''};
ss.labels  = { 'type 0' 'type 1' 'type 2' 'type 3' };

ss.values  = {0 1 2 3};
ss.val     = { 3 };
ss.help    = {  ' '
                'type 0: no small sample adjustment is used.'
                '             It tends to be biased and generally leads to overconfident inference in small sample.'
                'type 1: the errors used in the SwE estimation are multiplied by sqrt(n/(n-p)).'
                '             It tends to correct for the bias in small samples, but simulations seems to show that it still may leads to overconfident inference.'
                'type 2: the errors used in the SwE estimation are multiplied by 1/sqrt(1-h_ik).'
                '             It tends to correct for the bias in small samples, but, even if it generally performs better than the "type 1" adjustment, simulations seems to show that it still may leads to overconfident inference.'
                'type 3: the errors used in the SwE estimation are multiplied by 1/(1-h_ik).'
                '             It tends to correct for the bias in small samples, but simulations seems to show that it may leads to overconsevative inference.'
                'h_ik is the diagonal element of the hat matrix H=X''(X''X)^(-1)X.'
                ' '
                }';
% ---------------------------------------------------------------------
% dof_cl Degrees of freedom type
% ---------------------------------------------------------------------
dof_cl         = cfg_menu;
dof_cl.tag     = 'dof_cl';
dof_cl.name    = 'Degrees of freedom type';
dof_cl.labels  = { 'naive' 'approx' };
dof_cl.values  = { 0 1 };
dof_cl.val     = { 0 };
dof_cl.help    = {  ' '
                'naive: naive estimation of the degrees of freedom by m-p_B (the number of subject minus the number of pure between covariates).'
                '             this choice tends to overestimate the degrees of freedom, but reduce the quantity of images saved and the computation time.'
                'approx: degrees of freedom estimation with the estimate proposed in Guillaume et al. (in preparation).'
                '             This choice is not recommended for the classic SwE as, with this SwE version, it is generally overconservative and a large amount of variances/covariances images (sum_i n_i*(n_i+1)/2 images) need to be saved.'
                ' '
                }';
% ---------------------------------------------------------------------
% dof_mo Small Degrees of freedom type
% ---------------------------------------------------------------------
dof_mo         = cfg_menu;
dof_mo.tag     = 'dof_mo';
dof_mo.name    = 'Degrees of freedom type';
dof_mo.labels  = { 'naive' 'approx' };
dof_mo.values  = { 0 1 };
dof_mo.val     = { 1 };
dof_mo.help    = {  ' '
                'naive: naive estimation of the degrees of freedom by m-p_B (the number of subject minus the number of pure between covariates).'
                '             this choice tends to overestimate the degrees of freedom, but reduce the quantity of images saved and the computation time.'
                'approx: degrees of freedom estimation with the estimate proposed in Guillaume et al. (in preparation).'
                '             This choice is not recommended for the classic SwE as, with this SwE version, it is generally overconservative and a large amount of variances/covariances images (sum_i n_i*(n_i+1)/2 images) need to be saved.'
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
type.help    = {''};
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
                   'There are three methods for estimating global effects (1) Omit (assumming no other options requiring the global value chosen) (2) User defined (enter your own vector of global values) (3) Mean: SPM standard mean voxel value (within per image fullmean/8 mask) '
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
                 'When proportional scaling global normalisation is used each image is separately scaled such that it''s global value is that specified (in which case the grand mean is also implicitly scaled to that value). So, to proportionally scale each image so that its global value is eg. 20, select <Yes> then type in 20 for the grand mean scaled value.'
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
                   '''Ancova by subject'' or ''Ancova by effect'' options are implemented using the ANCOVA options provided where each experimental factor (eg. subject or effect), is defined. These allow eg. different subjects to have different relationships between local and global measurements. '
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
                   '''Ancova by subject'' or ''Ancova by effect'' options are implemented using the ANCOVA options provided where each experimental factor (eg. subject or effect), is defined. These allow eg. different subjects to have different relationships between local and global measurements. '
                   ''
                   'Since differences between subjects may be due to gain and sensitivity effects, AnCova by subject could be combined with "grand mean scaling by subject" (an option also provided where each experimental factor is originally defined) to obtain a combination of between subject proportional scaling and within subject AnCova. '
                   ''
}';
% ---------------------------------------------------------------------
% data Data & Design
% ---------------------------------------------------------------------
design        = cfg_exbranch;
design.tag    = 'design';
design.name   = 'Data & Design';
design.val    = {dir scans type subjects generic masking globalc globalm};
design.help   = {'Specify the data and design.'};
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
