function boxCoxTransformedData = swe_boxCoxTransform(data, boxCoxLambda)
  % Transform the data using a Box-cox transformation with parameter boxCoxLambda
  %
  % FORMAT boxCoxTransformedData = swe_boxCoxTransform(data, boxCoxLambda)
  %
  % data          - data to be Box-Cox transformed (must be positive)
  % boxCoxLambda  - lambda parameter of the Box-Cox transformation
  %
  % boxCoxTransformedData - Box-Cox transformed data
  %
  % =========================================================================
  % Bryan Guillaume
  % Version Info:  $Format:%ci$ $Format:%h$
  
  if min(data) < 0
    error('The data must be positive.')
  end

  if (boxCoxLambda == 0)
    boxCoxTransformedData = log(data);
  else
    boxCoxTransformedData = (data .^ boxCoxLambda - 1) / boxCoxLambda;
  end

  boxCoxTransformedData(data == 0) = -Inf;