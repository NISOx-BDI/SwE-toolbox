function result=runTest(porwb, torf, matorimg)
	
	% These warnings occur as we cannot open the displays.
	warning('off', 'SPM:noDisplay');
	warning('off','Octave:abbreviated-property-match');

	% Work out which test we are running.
	testname = [porwb '_' torf '_' matorimg];
	disp('==============================================================')
	disp(['Test case running: ' testname])
	disp('==============================================================')

	% Generate the results for the test.
	generateData(porwb, torf, matorimg);
	disp('==============================================================')
	disp(['Test case ' testname ' has been run.'])
	disp('--------------------------------------------------------------')
	disp('Verifying test results')
	disp('==============================================================')

	% Compare test results to ground truth.
	result = verifyMapsUnchanged();


	disp('==============================================================')
	if ~result
		disp('A test has failed!')
		error(['Test ' testname ' has failed.'])
	else
		disp('All tests pass!!')
	end
	disp('==============================================================')

end


function generateData(porwb, torf, matorimg)

	% Move into the test folder and add the path to tests.
	cd(['/swe/test/data/test_' porwb '_' torf '_' matorimg]);
	addpath('/swe');
	addpath('/swe/test');

	% Reset all seeds (in octave these are all different!!).
	load('/swe/test/data/seed.mat');
	rand('state',seed);
	randn('state', seed);
	randp('state', seed);
	randg('state',seed);
	rande('state',seed);

	% Load the test design and run it.
	load('design.mat');
	swe_run_design(design);

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

function mapsEqual = verifyMapsUnchanged()
	
	% List all files for testing
	if ~isempty(strfind(pwd, 'img'))
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