% Some notes on functions:
% which would you rather read: analyze_bandit.m vs avgspecgram.m
% it's not only that the latter is shorter -- it's documented, commented,
% and better-named

% why write functions? the same reason you want to buy generic lego blocks
% the benefit is composability and black-boxing
% a function lets us replace a block of code with a meaningful name and
% *encapsulate* functionality -- mentall easier

% in matlab:
% only one (main) function per m-file (boo!)
% function setup

function [output1, output2] = myfunction(input1, input2, etc)
%code goes here

% somewhere in the function, output1 and output2 must be set

end

% think of function as a "clean room" named variables passed in through
% slots (parameters), passed out through outputs
% everything else stays inside (scope!)

% exercise -- turn our psth code into a function
% can we turn it into smaller functions? why would we want to do this?
%  - reuse those functions in other functions
%  - swap those pieces out with something better later on

% advanced exercises: work through getting_started_spikes and
% getting_started_lfp
% look at code for functions called