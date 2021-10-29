%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code has been adapted from the ALEX toolbox and incorporated into this m-file.  
% The original free source code, which is copywritten by Claudio Rivetti and Mark Young 
% for ALEX is available at www.mathtools.net.
%
%This function/script is authorized for use in government and academic
%research laboratories and non-profit institutions only. Though this
%function has been tested prior to its posting, it may contain mistakes or
%require improvements. In exchange for use of this free product, we 
%request that its use and any issues relating to it be reported to us. 
%Comments and suggestions are therefore welcome and should be sent to 
%Prof. Robert Carpick <carpick@engr.wisc.edu>, Engineering Physics 
%Department, UW-Madison.
%Date posted: 7/8/2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function val = extract_num(str)

%Ascii table of relevant numbers
%character    ascii code
%  e          101
%  E          69
%  0          48
%  1          49
%  2          50
%  3          51 
%  4          52
%  5          53
%  6          54 
%  7          55 
%  8          56
%  9          57

eos = 0;
R = str;

while(~eos)
   
   [T,R] = strtok(str);
   if( length(R) == 0) eos = 1; end
   I = find( (T>=48) & (T<=57) | 101==T | 69==T | T==173 | T== 45 | T==46 | T==40);
  
   LT = length(T);
   LI = length(I);
   
   if( LI == LT )
      J = find(T~='(');
      val = str2num(T(J));
      break
   end
      str =R;
end


