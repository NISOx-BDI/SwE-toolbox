%==========================================================================
%Unit tests for testing whether all use cases run in the SwE toolbox. To 
%run the below run the runTest function. 
%
%Author: Tom Maullin (02/05/2018).
%==========================================================================
classdef test_swe_run < matlab.unittest.TestCase
    
    properties
       
        img_data
        mat_data
        
    end
    
    methods
        
        function delete_prev_run(testCase, data_path)
            
        end
        
    end
    
    methods(TestClassSetup)
        
        function getDataLocation(testCase)
            
            % Save the locations of the mat and img data.
            testDir = fileparts(which('test_swe_run.m'));
            mat_datapath = fullfile(testDir,'data',...
                                    'mat_input');
            img_datapath = fullfile(testDir,'data',...
                                    'img_input');
                                
            % List the files.
            testCase.img_data = arrayfun(@(x) fullfile(img_datapath,...
                    ['con_00' sprintf('%02d',x) '.img']),3:38,...
                'UniformOutput',false);
            testCase.mat_data ={fullfile(mat_datapath, 'subj_data.mat'),...
                                fullfile(mat_datapath, 'tris.mat')};
            
        end
        
    end
    
    methods(Test)
       
        function test_blah_runs(testCase)
            
            % The below is the job design.
            mkdir('C:\Users\TomM\Documents\Repositorys\SwE-toolbox\test\data\mat_input\temp')
            design.dir = {'C:\Users\TomM\Documents\Repositorys\SwE-toolbox\test\data\mat_input\temp'};
            design.scans = testCase.img_data;
            design.type.modified.groups = ones(36,1);
            design.type.modified.visits = kron((1:3),ones(1,12))';
            design.type.modified.ss = 3;
            design.type.modified.dof_mo = 3;
            design.subjects = kron(ones(1,3),1:12)';
            design.cov(1).c = kron([1 0 0],ones(1,12))';
            design.cov(1).cname = 'Cov1';
            design.cov(2).c = kron([0 1 0],ones(1,12))';
            design.cov(2).cname = 'Cov2';
            design.cov(3).c = kron([0 0 1],ones(1,12))';
            design.cov(3).cname = 'Cov3';
            design.masking.tm.tm_none = 1;
            design.masking.im = 1;
            design.masking.em = {''};
            design.WB.WB_yes.WB_type = 1;
            design.WB.WB_yes.WB_ss = 4;
            design.WB.WB_yes.WB_nB = 20;
            design.WB.WB_yes.WB_SwE = 0;
            design.WB.WB_yes.WB_stat.WB_T.WB_T_con = [1 0 0];
            design.WB.WB_yes.WB_cluster.WB_cluster_yes = 0.001;
            design.globalc.g_omit = 1;
            design.globalm.gmsca.gmsca_no = 1;
            design.globalm.glonorm = 1;
            
            % We now run the design.
            swe_run_design(design);
            load('C:\Users\TomM\Documents\Repositorys\SwE-toolbox\test\data\mat_input\temp\SwE.mat');
            swe_cp(SwE);
            
        end
        
    end
end