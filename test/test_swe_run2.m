function test_suite=my_test_of_abs
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;
end

function verifyMapsUnchanged(testLocation)

	

end

function test_example
	
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

	
end