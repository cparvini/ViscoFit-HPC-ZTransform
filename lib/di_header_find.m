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

% [ position ] = DI_header_find( file_name ,  find_string);
%
% Function DI_header_find finds all occurences of "find_string" in a Digital 
% Instruments file header from "file_name"
% Example 
% [location] = DI_header_find( '033115500.001' , 'Offline Planefit')



function [position ] = di_header_find( file_name ,  find_string);

% Define End of file identifier



% Opend the file given in argument and reference as
% fid.  Also if there was an error output error
% number and error message to screen
fid = fopen(file_name,'r');
[message,errnum] = ferror(fid);
if(errnum)
   fprintf(1,'I/O Error %d \t %s',[errnum,message]);
%   break
end


header_end=0;
eof = 0;
counter = 1;


byte_location = 0;


while( and( ~eof, ~header_end ) )
   
   byte_location = ftell(fid);
   line = fgets(fid);
   extract_num(int16(line))
   
   if( (-1)==line )
      eof  = 1;
      break
   end
     
   
   if( length( findstr(line,find_string) ) )
      position(counter) = byte_location;
      counter = counter + 1;
   end
   
      
   if length( findstr( line, '\*File list end' ) )
      header_end = 1;
   end
   
      
end


fclose(fid);

