%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% dost_tutorial (script)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% this script shows how to use the dost and the idost and how to visualize 
% them.
% For more informations on the dost, the idost and their visualization see 
% the comments on their code.
%
% Code by: U. Battisti and L. Riba 
% July, 02 2014
%
% References:
% [1]   R.G. Stockwell, "Why use the S-Transform", Pseudo-differential 
%       operators partial differential equations and time-frequency 
%       analysis, vol. 52 Fields Inst. Commun., pages 279--309, 
%       Amer. Math. Soc., Providence, RI 2007;
% [2]   R.G. Stockwell, "A basis for efficient representation of the
%       S-transform", Digital Signal Processing, 17: 371--393, 2007;
% [3]   Y. Wang and J. Orchard, "Fast-discrete orthonormal 
%       Stockwell transform", SISC: 31:4000--4012, 2009;
% [4]   Y. Wang, "Efficient Stockwell transform with applications to 
%       image processing", PhD thesis, University of Waterloo, 
%       Ontario Canada, 2011;
% [5]   U. Battisti, L.Riba, "Window-dependent bases for efficient 
%       representation of the Stockwell transform", 2014
%       http://arxiv.org/abs/1406.0513.
%
%
% Additional details:
% Copyright (c) by U. Battisti and L. Riba
% $Revision: 1.0 $  
% $Date: 2014/07/02  $
%
%
% Summary:
% 1) create a test-signal (we use the same one used in the References); 
% 2) compute the dost of the signal;
% 3) compute the inverse dost of the dost to reconstruct the signal;
% 4) plot the signal, the dost coefficients, the reconstructed signal and
%    the difference between reconstructed and original signal.


clear all
close all

% create the test-signal, ns is the number of samples
ns = 1024;
t = linspace(0,1024,ns);
test_signal  = cos( 2.*pi .* (204 + 6 .* ( ns./(1.+t) ).* ...
       cos(2.* pi .*(3 .* t./ns) )) .* t./ns );
% we want the test_signal to be a column vector   
test_signal = test_signal';   

% compute the dost coefficients of the signal "test_signal"
dost_test_signal = dost(test_signal);
% rearrange the dost coefficients in a more readable way using rearrange_dost. 
rearranged_dost_coefficients = rearrange_dost(dost_test_signal);

% reconstruct the signal using the idost
reconstructed_test_signal = idost(dost_test_signal);

figure
subplot(4,1,1)
plot(test_signal)
title('test signal')
axis tight
subplot(4,1,2)
imagesc([1 1024],[-.5 .5], flipud(abs(rearranged_dost_coefficients))) 
set(gca,'YDir','normal') 
title(' (norm) of the rearranged dost coefficients')
axis tight
subplot(4,1,3)
plot(real(reconstructed_test_signal))
title('reconstructed test signal')
axis tight
subplot(4,1,4)
plot(abs(reconstructed_test_signal-test_signal))
title('error')
axis tight

