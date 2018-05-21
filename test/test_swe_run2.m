function test_suite=my_test_of_abs
    try % assignment of 'localfunctions' is necessary in Matlab >= 2016
        test_functions=localfunctions();
    catch % no problem; early Matlab versions can use initTestSuite fine
    end
    initTestSuite;
end

function test_example
	addpath(genpath('/swe/'));
	load('test/data/seed.mat');
	rand('state',seed);
	load('test/data/design.mat');
	swe_run_design(design);
	load('SwE.mat');
	swe_cp_WB(SwE);

	
end