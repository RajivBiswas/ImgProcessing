To Run the Sources, plz follow the steps as follows:   

>1 Extract the sources.   


>2 Open MATLAB and first run "InputProgram.m", which will read an Input Image and extract the R, G, B Pixel Buffer and save it 
in "Info_R.bin", the Input Image is "Test_edge4.bmp", user needs to specify in the MATLAB Code the location from where it will 
pick the Input Images to Read and where the extracted Pixel Buffer would be stored.     



>3 Open Visual Studio and import the Visual Studio Project, "CannyEdgeFilter.vcproj". Change the values of " IWidth" and "IHeight" 
as per the Width and Height of the Input Image. 

Also, change the File Pointers to the Input Pixel Buffer extracted using the MATLAB program, "InputProgram.m", in the following line, as per user preference, "fp = fopen("F:\Info_R.bin", "rb")". Similarly, change the Output File Pointer, in the following line, as per user preference, "fp1 = fopen("F:\Info_R1.bin", "w")". This is the Output Pixel Buffer after Cannny Algorithm in Visual Studio, written in Info_R1.bin.    


>4 In MATLAB now run the "OutputProgram.m" and the Output is saved in the form of .bmp file. Nothing much happens here except rendering the Output using 'imshow' and a simple file write to write the Output in some specific Image format. Please take care in specifying the correct Width and Height of the Image, which is same as the Width and Height of the Input Image given.  

For any issues, plz mail me at: rajiv.biswas55@gmail.com

Regards, 
Rajiv.
