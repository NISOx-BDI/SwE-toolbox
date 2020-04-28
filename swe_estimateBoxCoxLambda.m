function boxCoxLambda = swe_estimateBoxCoxLambda(data)
  % Estimate the lambda of the Box-cox transformation of the data supplied
  %
  % FORMAT boxCoxLambda = swe_estimateBoxCoxLambda(data)
  %
  % data          - data to be Box-Cox transformed (must be positive)
  %
  % boxCoxLambda  - lambda parameter of the Box-Cox transformation
  %
  % =========================================================================
  % Bryan Guillaume
  % Version Info:  $Format:%ci$ $Format:%h$
  
  if min(data) <= 0
    error('The data must be positive.')
  end

  startingLambdaValue = 0;
  boxCoxLambda = fminsearch( @(lambda) getMinusBoxCoxLogLikelihood(data, lambda), startingLambdaValue);

  function minusBoxCoxLogLikelihood = getMinusBoxCoxLogLikelihood(data, lambda)    

    nData = numel(data);
    
    transformedData = swe_boxCoxTransform(data, lambda);

    minusBoxCoxLogLikelihood = 0.5 * nData * log(var(transformedData, 0)) + (1 - lambda) * sum(log(data));
    
  end
end