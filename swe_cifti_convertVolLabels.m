function convertedLabels = swe_cifti_convertVolLabels(originalLabels, abbreviate)
  % Convert CIFTI struture labels from long labels to abbreviations or convert abbreviations to long labels
  % FORMAT convertedLabels = swe_cifti_convertVolLabels(originalLabels, abbreviate)
  % originalLabels            - a cell array of cifti volume labels or a single cifti volume label
  % abbreviate   							- a boolean indicating if the labels needs to be abbreviated (true) or if the full labels need to be given (false)
  % convertedLabels      			- a cell array of converted cifti volume labels or a single cifti converted volume label
	ciftiVolumeLabels = {...
				'CIFTI_STRUCTURE_ACCUMBENS_LEFT', 'Acc_L';
				'CIFTI_STRUCTURE_ACCUMBENS_RIGHT', 'Acc_R';
				'CIFTI_STRUCTURE_ALL_WHITE_MATTER', 'WM';
				'CIFTI_STRUCTURE_ALL_GREY_MATTER', 'GM';
				'CIFTI_STRUCTURE_AMYGDALA_LEFT', 'Amg_L';
				'CIFTI_STRUCTURE_AMYGDALA_RIGHT', 'Amg_R';
				'CIFTI_STRUCTURE_BRAIN_STEM', 'BrStm';
				'CIFTI_STRUCTURE_CAUDATE_LEFT', 'Cau_L';
				'CIFTI_STRUCTURE_CAUDATE_RIGHT','Cau_R';
				'CIFTI_STRUCTURE_CEREBELLAR_WHITE_MATTER_LEFT', 'Cer_WM_L';
				'CIFTI_STRUCTURE_CEREBELLAR_WHITE_MATTER_RIGHT', 'Cer_WM_R';
				'CIFTI_STRUCTURE_CEREBELLUM', 'Cer';
				'CIFTI_STRUCTURE_CEREBELLUM_LEFT', 'Cer_L';
				'CIFTI_STRUCTURE_CEREBELLUM_RIGHT', 'Cer_R';
				'CIFTI_STRUCTURE_CEREBRAL_WHITE_MATTER_LEFT', 'Cbl_WM_L';
				'CIFTI_STRUCTURE_CEREBRAL_WHITE_MATTER_RIGHT', 'Cbl_WM_R';
				'CIFTI_STRUCTURE_CORTEX', 'Ctx';
				'CIFTI_STRUCTURE_CORTEX_LEFT', 'Ctx_L';
				'CIFTI_STRUCTURE_CORTEX_RIGHT', 'Ctx_R';
				'CIFTI_STRUCTURE_DIENCEPHALON_VENTRAL_LEFT', 'Dnc_V_L';
				'CIFTI_STRUCTURE_DIENCEPHALON_VENTRAL_RIGHT', 'Dnc_V_R';
				'CIFTI_STRUCTURE_HIPPOCAMPUS_LEFT', 'Hip_L';
				'CIFTI_STRUCTURE_HIPPOCAMPUS_RIGHT', 'Hip_R';
				'CIFTI_STRUCTURE_OTHER', 'Oth';
				'CIFTI_STRUCTURE_OTHER_GREY_MATTER', 'Oth_GM';
				'CIFTI_STRUCTURE_OTHER_WHITE_MATTER', 'Oth_WM';
				'CIFTI_STRUCTURE_PALLIDUM_LEFT', 'Pal_L';
				'CIFTI_STRUCTURE_PALLIDUM_RIGHT', 'Pal_R';
				'CIFTI_STRUCTURE_PUTAMEN_LEFT', 'Put_L';
				'CIFTI_STRUCTURE_PUTAMEN_RIGHT', 'Put_R';
				'CIFTI_STRUCTURE_THALAMUS_LEFT', 'Thl_L';
				'CIFTI_STRUCTURE_THALAMUS_RIGHT', 'Thl_R'};

originalClass = class(originalLabels);
if ~iscellstr(originalLabels)
	originalLabels = {originalLabels};
end
convertedLabels = cell(size(originalLabels));
if abbreviate
	for i = 1:numel(originalLabels)
		convertedLabels{i} = char(ciftiVolumeLabels(strcmpi(ciftiVolumeLabels(:,1), originalLabels(i)), 2));
	end
else
	for i = 1:numel(originalLabels)
		convertedLabels{i} = char(ciftiVolumeLabels(strcmpi(ciftiVolumeLabels(:,2), originalLabels(i)), 1));
	end
end

if strcmp(originalClass, 'char')
	convertedLabels = char(convertedLabels(1));
end
