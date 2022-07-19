

function [T,t,Err] = chopper(rg,evtt,fs,maxN) 

% [T,t,Err] = chopper(rg,evtt,[fs],[maxN]) 
%
% Creates a matrix to segment a signal sampled at fs into windows ranging from rg(1) to
% rg(end) around event times specified in evtt.
%
% INPUT VARIABLES:
%
%  rg: 2 element array containing start and end times (e.g. -0.5 to 1.0) in seconds
%
%  evtt:  time stamps (in s)
%
%  fs:  sampling rate. Default is 1.
%
%  maxN: Maximum index. If this value is greater than 0, then indices are wrapped
%        as mod(T-1,maxN)-1. Note that non-positive indices are therefore wrapped to the
%        end of the array while indices greater than maxN are wrapped to the beginning. 
%        If maxN is negative, then segments with out-of-range indices are
%        simply discarded. Default is 0 (do nothing).
%
% OUTPUT VARIABLE:
%
%  T: matrix of indices into a signal sampled at fs. Each column contains sample indices for 
%     the window surrounding the corresponding event time in evtt, as specified by rg. 
%     To segment the signal in the vector 'x' into a time-by-event matrix, use X = x(T), 
%     after correcting for any overlap with the ends of x. One way to implement such a 
%     correction is
%		T(T<1)=1; 
%		T(T>length(x)) = length(x);
%     which clamps values in X falling outside the time range of x to the first and last 
%     values of x, respectively.
%
%  t: vector, timestamps for the samples relative to evtt for the rows of T.
%  
%  Err: Difference between the sampled times and the exact window time.
% 

% ----------- SVN REVISION INFO ------------------
% $URL$
% $Revision$
% $Date$
% $Author$
% ------------------------------------------------

if nargin < 4
    maxN = false;
end

if isstruct(rg)
    if isfield(rg,'fs')
        fs = rg.fs;
    end
    if isfield(rg,'wint')
        evtt = rg.wint;
    end
    if isfield(rg,'maxN')
        maxN = rg.maxN;
    end
    rg = rg.Trange;
end    

t = (rg(1):1/fs:rg(2))';

T = fs*( repmat(t,1,length(evtt)) + repmat(evtt(:)',length(t),1))+1;



if nargout > 2
   Err = (round(T)-T)./fs;  %%% Return the rounding error 
end

T = round(T);

if maxN > 0
    T = mod(T-1,maxN)+1;
elseif maxN < 0
    T = T(:,all(T>=1) & all(T<= -maxN));
end