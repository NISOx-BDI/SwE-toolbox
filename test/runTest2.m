%==========================================================================
%This function runs the test_swe_run tests, which test whether various use
%cases of the SwE toolbox run without error.
%
%Author: Thomas Maullin.
%==========================================================================
function runTest()
    
    setup();
    
end

function setup()
    
    % Temporarily change path (in case currently in data directory).
    current = pwd;
    cd(fileparts(mfilename('fullpath')));

    % Download the '.img' test data and unzip it.
    zipped_data= urlwrite(['http://www.fil.ion.ucl.ac.uk/spm/download/data'...
                        '/face_rfx/face_rfx.zip'],...
                        fullfile(fileparts(mfilename('fullpath')),...
                            'data', 'temp.zip'));
    unzip(zipped_data, fileparts(mfilename('fullpath')));
    delete(zipped_data);
    
    % Remove unnecessary test data.
    unzipped_data = fullfile(fileparts(mfilename('fullpath')),'face_rfx');
    rmdir(fullfile(unzipped_data, 'cons_can'), 's');
    rmdir(fullfile(unzipped_data,'cons_fir'), 's');
    
    % We won't use all maps from con_informed.
    for i = 39:74
        delete(fullfile(unzipped_data,'cons_informed',...
                        ['con_00' num2str(i) '.hdr']));
        delete(fullfile(unzipped_data,'cons_informed',...
                        ['con_00' num2str(i) '.img']));
    end
    
    % We won't use any of the '.mat' files from here.
    for i = 3:74
        delete(fullfile(unzipped_data, 'cons_informed',...
                        ['con_00' sprintf('%02d',i) '.mat']));
    end
    
    % Move into a file named data for cleaner pathnames.
    if exist(fullfile(unzipped_data, '..', 'data'), 'dir')
        rmdir(fullfile(unzipped_data, '..', 'data'), 's');
	mkdir(fullfile(unzipped_data, '..', 'data'));
    end
    movefile(fullfile(unzipped_data, 'cons_informed'),...
             fullfile(unzipped_data, '..', 'data', 'img_input'));
    
    % Download the '.mat' test data (timeout for download might have to be
    % set a little higher than normal).
    if ~exist('OCTAVE_VERSION', 'builtin')
       weboptions('Timeout',10);
    end
    mkdir(fullfile(unzipped_data, '..', 'data', 'mat_input'));
    urlwrite(['https://drive.google.com/uc?export=download&id=1RXHFtnB1'...
             'N14-FcOuda8139zgR5dnjcly'],...
             fullfile(unzipped_data, '..', 'data', 'mat_input',...
                     'subj_data.mat'));
    urlwrite(['https://drive.google.com/uc?export=download&id=1vuHGkPvQMulp'...
         'nqnRF0L-BOQV0RyF3dlD'],...
         fullfile(unzipped_data, '..', 'data', 'mat_input',...
                 'tri_data.mat'));
     
     % Remove thd old data location and change back to working directory.
     rmdir(unzipped_data, 's');
     cd(current);
    
end
