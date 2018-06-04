function test_suite=test_swe_run2
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;
end

function mapsEqual = verifyMapsUnchanged()
	
	% List all files for testing
	if ~isempty(strfind(pwd, 'img'))
		files = ls("*.nii");
		filetype = 'nii';
	else
		files = ls("swe_*.mat");
		filetype = 'mat';
	end
	
	disp('filetype: ')
	disp(filetype)
	disp('files: ')
	disp(files)

	% store whether maps are equal.
	equalMaps = [];

	% Compare each file to ground truth
	for i = 1:size(files, 1)

		% Get the filenames
		file = files(i, :);
		gt_file = ['ground_truth' filesep file];
		disp(['Testing file: ' file])

		if strcmp(filetype, 'nii')

			% Read in the volumes
			file = spm_vol(file);
			file = spm_read_vols(file);
			gt_file = spm_vol(gt_file);
			gt_file = spm_read_vols(gt_file);

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

		% Remove NaN values
		file = file(~isnan(file));
		gt_file = gt_file(~isnan(gt_file));

		% Check whether the remaining values are within 
		% machine tolerance.
		result = ~any(abs(file-gt_file) > 5*eps);
		
		% Useful for debugging.
		if ~result
			disp('Length file: ')
			disp(size(file))
			
			disp('Length gt_file: ')
			disp(size(gt_file))
			
			disp('size disagreement values: ')
			disp(sum(file~=gt_file))
			
			disp('sum of disagreement (file): ')
			disp(sum(file(file~=gt_file)))
			
			disp('sum of disagreement (gt_file): ')
			disp(sum(gt_file(file~=gt_file)))
			
			disp('disagreement values (file)')
			d1 = file(file~=gt_file);
			disp(sprintf('%.9f', d1(1)))
			
			disp('disagreement values (gt_file)')
			d2 = gt_file(file~=gt_file);
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

function test_wb_t_img()
	
	disp('Test case running: wb_t_img')

	% Move into the test folder and add the path to tests.
	cd('/swe/test/data/test_wb_t_img');

	% Check against ground truth.
	mapsUnchanged = verifyMapsUnchanged();
	assertEqual(mapsUnchanged, true);
	
end

function test_wb_f_img()
	
	disp('Test case running: wb_f_img')

	% Move into the test folder and add the path to tests.
	cd('/swe/test/data/test_wb_f_img');

	% Check against ground truth.
	mapsUnchanged = verifyMapsUnchanged();
	assertEqual(mapsUnchanged, true);
	
end

function test_wb_t_mat()

	disp('Test case running: wb_t_mat')

	% Move into the test folder and add the path to tests.
	cd('/swe/test/data/test_wb_t_mat');

	% Check against ground truth.
	mapsUnchanged = verifyMapsUnchanged();
	assertEqual(mapsUnchanged, true);
	
end

function test_wb_f_mat()
	
	disp('Test case running: wb_f_mat')

	% Move into the test folder and add the path to tests.
	cd('/swe/test/data/test_wb_f_mat');

	% Check against ground truth.
	mapsUnchanged = verifyMapsUnchanged();
	assertEqual(mapsUnchanged, true);
	
end
