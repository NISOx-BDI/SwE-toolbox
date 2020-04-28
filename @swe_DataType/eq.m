function out = eq(a,b)
	% Overload the '==' operator for we_DataType objects
  % =========================================================================
  % Bryan Guillaume
	% Version Info:  $Format:%ci$ $Format:%h$
	
	out = strcmp(a.value, b.value);

end