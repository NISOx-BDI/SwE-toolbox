function varargout = swe_get_defaults(defstr, varargin)
% Get/set the defaults values associated with an identifier
% FORMAT defval = swe_get_defaults(defstr)
% Return the defaults value associated with identifier "defstr".
% Currently, this is a '.' subscript reference into the global
% "defaults" variable defined in swe_defaults.m.
%
% FORMAT swe_get_defaults(defstr, defval)
% Sets the defaults value associated with identifier "defstr". The new
% defaults value applies immediately to:
% * new modules in batch jobs
% * modules in batch jobs that have not been saved yet
% This value will not be saved for future sessions of SwE. To make
% persistent changes, edit swe_defaults.m.
%__________________________________________________________________________
% Based on swe_get_defaults.m
% Thomas Nichols
% Version Info:  $Format:%ci$ $Format:%h$

global SwEdefs;
if isempty(SwEdefs)
    swe_defaults;
end

% % construct subscript reference struct from dot delimited tag string
% tags = textscan(defstr,'%s', 'delimiter','.');
% subs = struct('type','.','subs',tags{1}');
subs = defstr;

if nargin == 1
    varargout{1} = getfield(SwEdefs, subs);
else
    SwEdefs = setfield(SwEdefs, subs, varargin{1});
end
