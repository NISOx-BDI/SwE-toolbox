function test_suite=test_swe_run2
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;
end

function mapsEqual = verifyMapsUnchanged()
	
	% List all niftis
	files = ls("*.nii");

	% store whether maps are equal.
	equalMaps = [];

	% Compare each file to ground truth
	for i = 1:size(files, 1)

		% Get the filenames
		file = files(1, :);
		gt_file = ['ground_truth' filesep file];
		disp(['Testing file: ' file])

		% Read in the volumes
		file = spm_vol(file);
		file = spm_read_vols(file);
		gt_file = spm_vol(gt_file);
		gt_file = spm_read_vols(gt_file);

		% Remove NaN values
		file = file(~isnan(file));
		gt_file = gt_file(~isnan(gt_file));

		% Check whether the remaining values are equal.
		result = ~any(file~=gt_file);

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

function test_wb_t_img()
	
	disp('Test case running: wb_t_img')

	% Move into the test folder and add the path to it.
	cd('/swe/test/data/test_wb_t_img');
	addpath(genpath('/swe/'));

	% Reset the seed
	load('/swe/test/data/seed.mat');
	rand('state',seed);

	% Load the test design and run it.
	load('design.mat');
	swe_run_design(design);

	% Load the generated SwE file and run it.
	load('SwE.mat');
	swe_cp_WB(SwE);

	% Check against ground truth.
	mapsUnchanged = verifyMapsUnchanged('/swe/test/data/test_wb_t_img');
	assertEqual(mapsUnchanged, true);
	
end