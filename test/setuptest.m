function setuptest(porwb, torf, matorimg)

	testname = [porwb '_' torf '_' matorimg]

	disp(['Test case generating: ' testname])

	% Move into the test folder and add the path to tests.
	cd(['/swe/test/data/test_' testname]);
	addpath('/swe');
	addpath('/swe/test');

	% Reset the seed
	load('/swe/test/data/seed.mat');
	rand('state',seed);

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