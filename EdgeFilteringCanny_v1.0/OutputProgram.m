%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  Output Program - To write & render the Output to a .bmp file
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Xwidth  = 690;%225;%480;%633;
Xheight = 690;%225;%360;%490;


fid     = fopen('F:\Info_R1.bin', 'r');
R       = uint8(zeros(Xwidth, Xheight));

     R = fscanf(fid, '%4d', [Xwidth, Xheight]);  
     R = uint8(R);
     R = R';
     fclose(fid);

imshow(R);
imwrite(R, 'F:\R1_new.bmp', 'bmp');