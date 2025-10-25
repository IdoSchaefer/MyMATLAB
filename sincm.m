function y=sincm(x)
% This is a slight modification of the sinc function of MATLAB to its
% common definition in mathematics. Only a matter of convenience.
%SINC Sin(x)/(x) function.
%   SINC(X) returns a matrix whose elements are the sinc of the elements 
%   of X, i.e.
%        y = sin(x)/(x)          if x ~= 0
%          = 1                   if x == 0
%   where x is an element of the input matrix and y is the resultant
%   output element.
%
%   % Example of a sinc function for a linearly spaced vector:
%   t = linspace(-5,5);
%   y = sinc(t);
%   plot(t,y);
%   xlabel('Time (sec)');ylabel('Amplitude'); title('Sinc Function')
%
%   See also SQUARE, SIN, COS, CHIRP, DIRIC, GAUSPULS, PULSTRAN, RECTPULS,
%   and TRIPULS.

%   Author(s): T. Krauss, 1-14-93
%   Copyright 1988-2004 The MathWorks, Inc.
%   $Revision: 1.7.4.1 $  $Date: 2004/08/10 02:11:27 $

i=find(x==0);                                                              
x(i)= 1;      % From LS: don't need this is /0 warning is off                           
y = sin(x)./(x);                                                     
y(i) = 1;   

