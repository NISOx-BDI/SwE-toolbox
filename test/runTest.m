function result = runTest(pOrWb, inferenceType, tOrF, matNiiGiiOrCii)
  % This function runs an octave test case for the SwE toolbox.
  % ====================================================================
  % FORMAT result=runTest(pOrWb, inferenceType, tOrF, matNiiGiiOrCii)
  % --------------------------------------------------------------------
  % Inputs:
  %
  %   - pOrWb:    			'p' for parametric or 'wb' for wild bootstrap.
  %		- inferenceType		'vox' for voxel-wise, 'dat' for .mat element-wise
  %											'dpx' for CIfTI/GIfTI element-wise, 'clus' for clusterwise or
  %											'tfce' for TFCE
  %   - tOrF:     			't' for a t test or 'f' for an f test.
  %   - matNiiGiiOrCii: 'mat' for '.mat' inputs, 'nii' for NIfTI inputs,
  %                				'gii' for GIfTI inputs or 'cii' for CIfTI inputs
  % ====================================================================
  % Tom Maullin and Bryan Guillaume
  % Version Info:  $Format:%ci$ $Format:%h$

  % Turn off warnings.
  %
  % (Footnote: These warnings occur as we cannot open the displays and octave reads
  % throws warnings about how it reads in '.nii' files.)
  warning('off','SPM:noDisplay');
  warning('off','Octave:abbreviated-property-match');
  warning('off','Octave:num-to-str');

  % Disable random number seeding
  global SwEdefs
  SwEdefs.shuffle_seed = false;

  % Work out which test we are running.
  testname = [pOrWb '_' inferenceType '_' tOrF '_' matNiiGiiOrCii];

  % Tell the user which test case we are running.
  disp(sprintf('\n=============================================================='))
  disp(['Test case running: ' testname])
  disp(sprintf('==============================================================\n'))

  % Test setup
  testSetup(pOrWb, inferenceType, tOrF, matNiiGiiOrCii)

  % Generate the results for the test.
  generateData(pOrWb, inferenceType, tOrF, matNiiGiiOrCii);


  % Tell the user we have run the test.
  disp(sprintf('\n=============================================================='))
  disp(['Test case ' testname ' has been run.'])
  disp('--------------------------------------------------------------')
  disp('Verifying test results')
  disp(sprintf('==============================================================\n'))

  % Compare test results to ground truth.
  result = verifyMapsUnchanged(pOrWb, inferenceType, tOrF, matNiiGiiOrCii);

  % Tell the user whether the tests passed.
  disp(sprintf('\n=============================================================='))
  if ~result
    disp('A test has failed!')
    error(['Test ' testname ' has failed.'])
  else
    disp('All tests pass!!')
  end 
  disp(sprintf('==============================================================\n'))

  % Teardown method
  testTearDown(pOrWb, inferenceType, tOrF, matNiiGiiOrCii);

end

function testSetup(pOrWb, inferenceType, tOrF, matNiiGiiOrCii)

  % Move into the test folder and add the path to tests.
  cd(['/swe/test/data/test_' pOrWb '_' inferenceType '_' tOrF '_' matNiiGiiOrCii]);
  ls;
  addpath('/swe');
  addpath('/swe/test');
    
  % Run teardown method just in case some of the files from the previous 
  % run managed to get cached.
  testTearDown(pOrWb, inferenceType, tOrF, matNiiGiiOrCii);

  % Set RNG seed to fixed value
  load('/swe/test/data/seed.mat');
  swe_seed(seed)

  % Make a copy of the original xSwE object for future runs.
  if exist('xSwE.mat')~=0
    copyfile('xSwE.mat', 'xSwE_orig.mat')
  end

end

function generateData(pOrWb, inferenceType, tOrF, matNiiGiiOrCii)

  design = createDesignFromTemplate(pOrWb, inferenceType, tOrF, matNiiGiiOrCii);
  swe_run_smodel(design);

  if strcmp(pOrWb, 'wb')

    % Load the generated SwE file and run it.
    load('SwE.mat');
    swe_cp_WB(SwE);

  else

    % Load the generated SwE file and run it.
    load('SwE.mat');
    swe_cp(SwE);

    % Define a contrast.
    load('SwE.mat');
    
    load('/swe/test/data/xCon.mat');
    
    if strcmp(tOrF, 't')
      xCon.name = 'T contrast 1';
      xCon.STAT = 'T';
      xCon.c = [1 0 0]'; 
    else
      xCon.name = 'F contrast 1';
      xCon.STAT = 'F';
      xCon.c = eye(3); 
    end

    SwE.xCon = xCon;
    
    save('SwE.mat', 'SwE');
    
    % Run swe_getSPM(). For the nii
    % case we use the xSwE object to
    % avoid user input.
    if strcmp(matNiiGiiOrCii, 'nii') || strcmp(matNiiGiiOrCii, 'gii') || strcmp(matNiiGiiOrCii, 'cii')
      
      load('/swe/test/data/xSwE.mat')
      
      xSwE.swd = SwE.swd;
      xSwE.STAT = upper(tOrF);
      xSwE.S = SwE.xVol.S;
      xSwE.M = SwE.xVol.M;
      xSwE.iM = SwE.xVol.iM;
      xSwE.DIM = SwE.xVol.DIM;
      M    = SwE.xVol.M(1:3,1:3);
      VOX  = sqrt(diag(M'*M))';
      xSwE.VOX = VOX;
      
      swe_getSPM(xSwE);
    
    else
    
      swe_getSPM();
    
    end

  end

end

function design = createDesignFromTemplate(pOrWb, inferenceType, tOrF, matNiiGiiOrCii)
  
  % Load template design
  load('/swe/test/data/design_wb_clus_f_mat.mat');
  
  design.dir = {['/swe/test/data/test_' pOrWb '_' inferenceType '_' tOrF '_' matNiiGiiOrCii]};
  
  if strcmp(matNiiGiiOrCii, 'nii')
    fileExtension = 'img';
  elseif strcmp(matNiiGiiOrCii, 'cii')
    fileExtension = 'dscalar.nii';
  else
    fileExtension = matNiiGiiOrCii;
  end

  if ~strcmp(matNiiGiiOrCii, 'mat')
    for i = 1:36
      design.scans{i} = sprintf('/swe/test/data/%s_input/con_%04i.%s', matNiiGiiOrCii, i+2, fileExtension);
    end
  end

  design.giftiAdditionalInfo = struct;
  design.giftiAdditionalInfo.areaFileForGiftiInputs = '';

  if strcmp(matNiiGiiOrCii, 'cii')
    design.ciftiAdditionalInfo = struct;
    
    design.ciftiAdditionalInfo.ciftiGeomFile(1) = struct;
    design.ciftiAdditionalInfo.ciftiGeomFile(1).brainStructureLabel = 'CIFTI_STRUCTURE_CORTEX_LEFT';
    design.ciftiAdditionalInfo.ciftiGeomFile(1).geomFile = '/swe/test/data/gii_input/L.sphere.4k_fs_LR.surf.gii';
    
    design.ciftiAdditionalInfo.volRoiConstraint = 1;
  end

  if strcmp(pOrWb, 'p')
    design.WB = rmfield(design.WB, 'WB_yes');
    design.WB.WB_no = 0;
  else
    
    if strcmp(inferenceType, 'clus')
      
      if ~strcmp(matNiiGiiOrCii, 'mat')
        design.WB.WB_yes.WB_infType.WB_clusterwise.WB_inputType = ...
          rmfield(design.WB.WB_yes.WB_infType.WB_clusterwise.WB_inputType, 'WB_mat');
        design.WB.WB_yes.WB_infType.WB_clusterwise.WB_inputType.WB_img = 0;
      end

    elseif strcmp(inferenceType, 'tfce')
      
      design.WB.WB_yes.WB_infType = rmfield(design.WB.WB_yes.WB_infType, 'WB_clusterwise');
      design.WB.WB_yes.WB_infType.WB_TFCE = struct;
      design.WB.WB_yes.WB_infType.WB_TFCE.WB_TFCE_E = 0.5;
      design.WB.WB_yes.WB_infType.WB_TFCE.WB_TFCE_H = 2;

    else

      design.WB.WB_yes.WB_infType = rmfield(design.WB.WB_yes.WB_infType, 'WB_clusterwise');
      design.WB.WB_yes.WB_infType.WB_voxelwise = 0;

    end

    if strcmp(tOrF, 't')
      design.WB.WB_yes.WB_stat = rmfield(design.WB.WB_yes.WB_stat, 'WB_F');
      design.WB.WB_yes.WB_stat.WB_T = struct;
      design.WB.WB_yes.WB_stat.WB_T.WB_T_con = [1, 0, 0];
    end

  end

end

function mapsEqual = verifyMapsUnchanged(pOrWb, inferenceType, tOrF, matNiiGiiOrCii)
  
  % List all files for testing
  if strcmp(matNiiGiiOrCii, 'nii')
    files = ls("*.nii");
    filetype = 'nii';
  elseif strcmp(matNiiGiiOrCii, 'gii')
    files = ls("*.gii");
    filetype = 'gii';
  elseif strcmp(matNiiGiiOrCii, 'cii')
    files = ls("*.nii");
    filetype = 'cii';  
  else
    files = ls("swe_*.mat");
    filetype = 'mat';
  end

  % store whether maps are equal.
  equalMaps = [];

  % Compare each file to ground truth
  for i = 1:size(files, 1)

    % Get the filenames
    file = files(i, :);
    gt_file = ['ground_truth' filesep file];
    disp(['Testing file: ' file])

    if ~strcmp(filetype, 'mat')

      % Read in the volumes
      file = swe_data_hdr_read(file);
      file = swe_data_read(file);
      gt_file = swe_data_hdr_read(gt_file);
      gt_file = swe_data_read(gt_file);

    else

      % Read in the surface data.
      file = load(strrep(file, " ", ""));
      gt_file = load(strrep(gt_file, " ", ""));

      % Retrieve field name.
      fieldname = fieldnames(file){1};

      % Retrieve data
      file = getfield(file, fieldname);
      gt_file = getfield(gt_file, fieldname);

    end

    % Remove NaN and zero values
    file = file(~isnan(file) & file ~= 0);
    gt_file = gt_file(~isnan(gt_file) & gt_file ~= 0);
    
    % Check whether the remaining values are equal.
    %
    % Footnote: There is some form of machine tolerance
    % error with the mat cases. Even with all seeds
    % reset, differences in voxel values of around ~e16
    % ~e16 can occur occasionally on runs. These errors 
    % can propogate and grow as large as  ~e8. To 
    % counter this, in these cases we just look to see
    % if we are within e-10 of the ground truth.
    % Suspected cause: The `beta` and `betainc` 
    % functions which are only run in these cases and 
    % are built on old fortran numeric approximations 
    % in octave. (Tom Maullin 07/06/2018)
    if strcmp(matNiiGiiOrCii, 'mat')
      result = ~any(abs(file-gt_file) > 10^(-7));
      indexWrong = abs(file-gt_file) > 10^(-7);
    else
      result = ~any(file~=gt_file)
      indexWrong = file~=gt_file;
    end
    
    % Useful for debugging.
    if ~result
      disp('Length file: ')
      disp(size(file))
      
      disp('Length gt_file: ')
      disp(size(gt_file))
      
      disp('size disagreement values: ')
      disp(sum(indexWrong))
      
      disp('sum of disagreement (file): ')
      disp(sum(file(indexWrong)))
      
      disp('sum of disagreement (gt_file): ')
      disp(sum(gt_file(indexWrong)))
      
      disp('disagreement values (file)')
      d1 = file(indexWrong);
      disp(sprintf('%.9f', d1(1)))
      
      disp('disagreement values (gt_file)')
      d2 = gt_file(indexWrong);
      disp(sprintf('%.9f', d2(1)))
      
      disp('diff')
      disp(d2(1)-d1(1))
      
    end

    if result
      disp('PASS');
    else
      disp('FAIL');
    end

    % Record result
    equalMaps = [equalMaps result];

    end

    % return true if all maps were equal false otherwise.

    mapsEqual = ~any(~equalMaps);

end

function testTearDown(pOrWb, inferenceType, tOrF, matNiiGiiOrCii)
  
  % Delete all files from this run.
  delete('swe_*');
  delete('SwE*');

  if exist('xSwE_orig.mat')~=0
    delete('xSwE.mat');
    rename('xSwE_orig.mat', 'xSwE.mat');
  end

end