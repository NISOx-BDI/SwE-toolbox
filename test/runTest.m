function result=runTest(porwb, inftype, torf, matorimg)
% This function runs an octave test case for the SwE toolbox.
% ====================================================================
% FORMAT result=runTest(porwb, inftype, torf, matorimg)
% --------------------------------------------------------------------
% Inputs:
%
%   - porwb:    'p' for parametric or 'wb' for wild bootstrap.
%   - torf:     't' for a t test or 'f' for an f test.
%   - matorimg: 'mat' for '.mat' file surface input pr 'img' for 
%                '.img' and '.hdr' volumetric input.
% ====================================================================
% Author: Tom Maullin (08/06/20178) 

	% Turn off warnings.
	%
	% (Footnote: These warnings occur as we cannot open the displays and octave reads
	% throws warnings about how it reads in '.img' files.)
	warning('off','SPM:noDisplay');
	warning('off','Octave:abbreviated-property-match');
	warning('off','Octave:num-to-str');

	% Disable random number seeding
	global SwEdefs
	SwEdefs.shuffle_seed = false;

	% Work out which test we are running.
	testname = [porwb '_' inftype '_' torf '_' matorimg];

	% Tell the user which test case we are running.
	disp(sprintf('\n=============================================================='))
	disp(['Test case running: ' testname])
	disp(sprintf('==============================================================\n'))

	% Test setup
	testSetup(porwb, inftype, torf, matorimg)

	% Generate the results for the test.
	generateData(porwb, inftype, torf, matorimg);


	% Tell the user we have run the test.
	disp(sprintf('\n=============================================================='))
	disp(['Test case ' testname ' has been run.'])
	disp('--------------------------------------------------------------')
	disp('Verifying test results')
	disp(sprintf('==============================================================\n'))

	% Compare test results to ground truth.
	result = verifyMapsUnchanged(porwb, inftype, torf, matorimg);

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
	testTearDown(porwb, inftype, torf, matorimg)

end

function testSetup(porwb, inftype, torf, matorimg)

	% Move into the test folder and add the path to tests.
	cd(['/swe/test/data/test_' porwb '_' inftype '_' torf '_' matorimg]);
	ls;
	addpath('/swe');
	addpath('/swe/test');
    
    % Run teardown method just in case some of the files from the previous 
    % run managed to get cached.
    testTearDown(porwb, inftype, torf, matorimg);

	% Set RNG seed to fixed value
	load('/swe/test/data/seed.mat');
	swe_seed(seed)

	% Make a copy of the original xSwE object for future runs.
	if exist('xSwE.mat')~=0
		copyfile('xSwE.mat', 'xSwE_orig.mat')
	end

end

function generateData(porwb, inftype, torf, matorimg)


	% Load the test design and run it.
	load('design.mat');
	swe_run_smodel(design);

	if strcmp(porwb, 'wb')

		% Load the generated SwE file and run it.
		load('SwE.mat');
		swe_cp_WB(SwE);

	else

		% Load the generated SwE file and run it.
		load('SwE.mat');
		swe_cp(SwE);

		% Define a contrast.
		load('SwE.mat');
		load('xCon.mat');
		SwE.xCon = xCon;
	    save('SwE.mat', 'SwE');
	    
	    % Run swe_getSPM(). For the img
	    % case we use the xSwE object to
	    % avoid user input.
	    if strcmp(matorimg, 'img')
	    	load('xSwE.mat')
	    	swe_getSPM(xSwE);
	    else
	    	swe_getSPM();
	    end

	end

end

function mapsEqual = verifyMapsUnchanged(porwb, inftype, torf, matorimg)
	
	% List all files for testing
	if strcmp(matorimg, 'img')
		files = ls("*.nii");
		filetype = 'nii';
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
		gt_file = strrep(gt_file, '_dat_', '_vox_');
		disp(['Testing file: ' file])

		if strcmp(filetype, 'nii')

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
		if strcmp(matorimg, 'mat')
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

function testTearDown(porwb, inftype, torf, matorimg)
	
	% Delete all files from this run.
	delete('swe_*')
	delete('SwE*')

	if exist('xSwE_orig.mat')~=0
		delete('xSwE.mat')
		rename('xSwE_orig.mat', 'xSwE.mat')
	end

end