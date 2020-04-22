function clusterStatistics = swe_getClusterStatistics(dataType, locationActivatedElements, dataTypeSpecificInformation, giftiAreaFile)
  % Compute the cluster statitics of surviving data elements
  % FORMAT clusterStatistics = swe_cifti_clusters(dataType, locationActivatedElements, dataTypeSpecificInformation, giftiAreaFile)
  % dataType      							- one of the type of data defined in swe_DataType
  % locationActivatedElements 	- an array of locations
	% dataTypeSpecificInformation - a variable containing additional information specific to each type of data as described below:
	%		- For swe_DataType.Nifti or swe_DataType.VolumeMat, this variable is not used
	%		- For swe_DataType.Gifti or swe_DataType.SurfaceMat, a [mx3] faces array or a patch structure
	%		- For swe_DataType.Cifti, a variable containing cifti specific information
	% giftiAreaFile								- the name of a GIfTI file containing the surface area of each vertex (optional and only used for swe_DataType.Gifti or swe_DataType.SurfaceMat)
	% clusterStatistics						- a struct containing several cluster statistics
  % =========================================================================
  % Bryan Guillaume
	% Version Info:  $Format:%ci$ $Format:%h$

	clusterStatistics = struct;

	switch dataType

		case {swe_DataType.Nifti, swe_DataType.VolumeMat}

			clusterStatistics.clusterAssignment = spm_clusters(locationActivatedElements);

		case {swe_DataType.Gifti, swe_DataType.SurfaceMat}

			[clusterStatistics.clusterAssignment, ~, clusterAreas] = swe_mesh_clusters(dataTypeSpecificInformation, locationActivatedElements(1,:), giftiAreaFile);
			
			if ~isnan(clusterStatistics.clusterAreas)
				clusterStatistics.clusterAreas = clusterAreas;
				clusterStatistics.maxClusterArea = max(clusterStatistics.clusterAreas);
			end

		case swe_DataType.Cifti
			
			[clusterStatistics.clusterAssignment, ...
				clusterStatistics.clusterSizesInSurfaces, ...
				clusterStatistics.clusterSizesInVolume] = ...
				swe_cifti_clusters(dataTypeSpecificInformation, locationActivatedElements(1,:));

			if isempty(clusterStatistics.clusterSizesInSurfaces)
				clusterStatistics.maxClusterSizeInSurfaces = 0;
			else 
				clusterStatistics.maxClusterSizeInSurfaces = max(clusterStatistics.clusterSizesInSurfaces);
			end

			if isempty(clusterStatistics.clusterSizesInVolume)
				clusterStatistics.maxClusterSizeInVolume = 0;
			else 
				clusterStatistics.maxClusterSizeInVolume = max(clusterStatistics.clusterSizesInVolume);
			end

		otherwise
			error('The data type is not recognised');
	end

	if isnan(clusterStatistics.clusterAssignment)
		clusterStatistics.clusterAssignment = [];
	end

	clusterStatistics.nCluster = max(clusterStatistics.clusterAssignment);
	clusterStatistics.clusterSize = histc(clusterStatistics.clusterAssignment, 1:clusterStatistics.nCluster);
	clusterStatistics.maxClusterSize = max(clusterStatistics.clusterSize);

end