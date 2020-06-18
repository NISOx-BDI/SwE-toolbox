function [clusterAssignment, clusterNumbersOfVertices, clusterAreas] = swe_mesh_clusters(M, activatedVertices, areaFile)
	% Label connected components of surface mesh data
	% FORMAT [clusterAssignment, clusterExtentInVertices, clusterAreas] = swe_mesh_clusters(M, activatedVertices, areaFile)
	% M        						- a [m x 3] faces array or a patch structure
	% activatedVertices 	- a [nVertices x 1] data vector (using NaNs or logicals)
	% areaFile 					 	- an optional file containing the surface area of each vertex
	%
	% clusterAssignment 				- a [1 x nActivatedVertices] vector of cluster indices
	% clusterNumbersOfVertices 	- a [1 x nClusters] size of clusters {in vertices}
	% clusterAreas 							- a [1 x nClusters] vecot of cluster areas
	% =========================================================================
	% Bryan Guillaume
	% Version Info:  $Format:%ci$ $Format:%h$

	[clusterAssignment, clusterNumbersOfVertices] = spm_mesh_clusters(M, activatedVertices);
	clusterAssignment = clusterAssignment(activatedVertices)';
	clusterNumbersOfVertices = clusterNumbersOfVertices';

	if nargin > 2 && ~isempty(areaFile)
		area = swe_data_read(areaFile, find(activatedVertices));
		nClusters = numel(clusterExtents);
		clusterAreas = nan(1, nClusters)
		for iCluster = 1:nClusters
			clusterAreas(iCluster) = sum(area(clusterAssignment == iCluster));
		end
	else
		clusterAreas = NaN;
	end

end
