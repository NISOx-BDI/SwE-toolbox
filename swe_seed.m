function swe_seed(varargin)
% Sets the random number generator seed
% =========================================================================
% Randomly or deterministically seed the random number generator
% =========================================================================
% FORMAT swe_seed
% -------------------------------------------------------------------------
% Randomly seed the random number generator. If SwE default "shuffle_seed"
% is false no action is taken.
% =========================================================================
% FORMAT swe_seed(Seed)
% -------------------------------------------------------------------------
% Inputs:
%   - Seed: Value for random number generator seed.
% Deterministically seed the random number generator with value Seed.  
% Useful for creating reproducible random numbers for testing; then make 
% sure that "shuffle_seed" default is false to prevent re-seeding.
% =========================================================================
% T. Nichols
% Version Info:  $Format:%ci$ $Format:%h$

isOctave = exist('OCTAVE_VERSION','builtin');

if nargin==0

  % Random seeding
  %----------------------------------------------------------------------
  if swe_get_defaults('shuffle_seed')
    if isOctave
      rand( 'state','reset');
      randn('state','reset');
      rande('state','reset');
      randg('state','reset');
      randp('state','reset');
    else
      rng('default');
      rng('shuffle');
    end
  end

elseif nargin==1

  % Deterministic seeding
  %----------------------------------------------------------------------
  Seed=varargin{1};
  if isOctave
    rand( 'state',Seed);
    randn('state',Seed);
    rande('state',Seed);
    randg('state',Seed);
    randp('state',Seed);
  else
    rng(Seed(1));
  end

else

  error('swe_seed: Wrong number of arguments')

end

