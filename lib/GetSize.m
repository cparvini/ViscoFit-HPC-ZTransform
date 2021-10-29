function totSize = GetSize(this) 
   props = whos('this'); 
   totSize = props.bytes;
%    fprintf(1, '%d bytes\n', totSize); 
end