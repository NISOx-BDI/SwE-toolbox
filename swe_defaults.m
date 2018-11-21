%
% FORMAT swe_defaults
%_______________________________________________________________________
%
% This file is intended to be customised for the site.
% Individual users can make copies which can be stored in their own
% matlab subdirectories. If ~/matlab is ahead of the SwE directory
% in the MATLABPATH, then the users own personal defaults are used.
%
% Care must be taken when modifying this file
%_______________________________________________________________________
% Based on snpm_defaults.m
% Thomas Nichols
% Version Info:  $Format:%ci$ $Format:%h$

global SwEdefs

% If true, shuffles the seed of the random number generator to get 
% different results every time. Use false, if you want to specify your own 
% seed, for instance to insure that results can be replicated or when using 
% a high performance cluster.
%------------------------------------------------------------------------
SwEdefs.shuffle_seed = true; 

