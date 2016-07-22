%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Input Program - To Read the Image & write the R/G/B-pixel Buffer, HxW
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

A       = imread('F:\Test_edge4.bmp');
fid     = fopen('F:\Info_R.bin', 'w');

R       = A(:,:,1);
G       = A(:,:,2);
B       = A(:,:,3);

[Xwidth, Xheight]          = size(R);
for n = 1:Xheight
    for m = 1:Xwidth
     fprintf(fid, '%4d', R(n, m));    
    end    
    fprintf(fid, '\n');
end

fclose(fid);
%%
%%
