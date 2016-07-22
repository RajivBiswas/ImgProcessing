
/***************************************************************
 *
 * This file is part of Image Processing Projects, by Rajiv  
 * Biswas, 
 *
 * Copyright (C) 2016  Rajiv Biswas
 *
 * This program is free software: you can redistribute it and/or 
 * modify it under the terms of the GNU General Public License as 
 * published by the Free Software Foundation, either version 3 of 
 * the License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be  
 * useful,but WITHOUT ANY WARRANTY; without even the implied 
 * warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
 * PURPOSE.  See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public 
 * License along with this program.  If not, see 
 *      <http://www.gnu.org/licenses/>.
 *
 *
 **************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <malloc.h>

#define IWidth  690//480//225//480//633
#define IHeight 690//480//225//360//490
#define FFABS(a) ((a) >= 0 ? (a) : (-(a)))

enum {
    DIRECTION_45UP,
    DIRECTION_45DOWN,
    DIRECTION_HORIZONTAL,

    DIRECTION_VERTICAL,
};

static int get_rounded_direction(int gx, int gy)
{
    /* reference angles:
     *   tan( pi/8) = sqrt(2)-1
     *   tan(3pi/8) = sqrt(2)+1
     * Gy/Gx is the tangent of the angle (theta), so Gy/Gx is compared against
     * <ref-angle>, or more simply Gy against <ref-angle>*Gx
     *
     * Gx and Gy bounds = [-1020;1020], using 16-bit arithmetic:
     *   round((sqrt(2)-1) * (1<<16)) =  27146
     *   round((sqrt(2)+1) * (1<<16)) = 158218
     */

    if (gx) {
        int tanpi8gx, tan3pi8gx;

        if (gx < 0)
            gx = -gx, gy = -gy;
        gy <<= 16;
        tanpi8gx  =  27146 * gx;
        tan3pi8gx = 158218 * gx;
        if (gy > -tan3pi8gx && gy < -tanpi8gx)  return DIRECTION_45UP;
        if (gy > -tanpi8gx  && gy <  tanpi8gx)  return DIRECTION_HORIZONTAL;
        if (gy >  tanpi8gx  && gy <  tan3pi8gx) return DIRECTION_45DOWN;
    }
    return DIRECTION_VERTICAL;
}

static void non_maximum_suppression(int w, int h,
                                          unsigned int  **dst,
                                    const int  **dir,
                                    const unsigned int **src)
{
    int i, j;

#define COPY_MAXIMA(ay, ax, by, bx) do {                \
    if (src[j][i] > src[j+ay][i+ax] &&     \
        src[j][i] > src[j+by][i+bx])       \
		dst[j][i] = (src[j][i] & (~0xFF)) ? ((-src[j][i])>>31) : src[j][i];                 \
} while (0)		

    for (j = 1; j < h - 1; j++) {
        for (i = 1; i < w - 1; i++) {
            switch (dir[j][i]) {
            case DIRECTION_45UP:        COPY_MAXIMA( 1, -1, -1,  1); break;
            case DIRECTION_45DOWN:      COPY_MAXIMA(-1, -1,  1,  1); break;
            case DIRECTION_HORIZONTAL:  COPY_MAXIMA( 0, -1,  0,  1); break;
            case DIRECTION_VERTICAL:    COPY_MAXIMA(-1,  0,  1,  0); break;
            }
        }
    }
}

static void sobel(int w, int h,
                       unsigned int **dst,
                         int **dir,
                  unsigned int **src)
{
    int i, j;

    for (j = 1; j < h-1; j++)
    {
	for (i = 1; i < w-1; i++)
	{
            int gx =
		-1*src[j-1][i-1] + 1*src[j-1][i+1]
		-2*src[j][i-1] + 2*src[j][i+1]
		-1*src[j+1][i-1] + 1*src[j+1][i+1];
            int gy =
		-1*src[j-1][i-1] + 1*src[j+1][i-1]
		-2*src[j-1][i] + 2*src[j+1][i]
		-1*src[j-1][i+1] + 1*src[j+1][i+1];

            dst[j][i] = FFABS(gx) + FFABS(gy);
	    dir[j][i] = get_rounded_direction(gx, gy);
	}
    }
}

/*

Sigma and the Gaussian Filter Kernels gnerated using MATLAB

PSF = fspecial('gaussian',5,1.0)
round(386*PSF)
     1     5     8     5     1
     5    23    38    23     5
     8    38    63    38     8
     5    23    38    23     5
     1     5     8     5     1
Sum = 383

PSF = fspecial('gaussian',5,1.2)
round(256*PSF)
     2     5     8     5     2
     5    15    21    15     5
     8    21    30    21     8
     5    15    21    15     5
     2     5     8     5     2
Sum = 254

PSF = fspecial('gaussian',5,1.4)
round(159*PSF)
     2     4     5     4     2
     4     9    12     9     4
     5    12    15    12     5
     4     9    12     9     4
     2     4     5     4     2
Sum = 159

PSF = fspecial('gaussian',5,1.6)
round(256*PSF)
     4     8     9     8     4
     8    14    17    14     8
     9    17    20    17     9
     8    14    17    14     8
     4     8     9     8     4
Sum = 260

PSF = fspecial('gaussian',5,1.8)
round(256*PSF)
     5     8    10     8     5
     8    13    15    13     8
    10    15    18    15    10
     8    13    15    13     8
     5     8    10     8     5
Sum = 254

PSF = fspecial('gaussian',5,2.0)
round(256*PSF)
     6     9    10     9     6
     9    13    14    13     9
    10    14    16    14    10
     9    13    14    13     9
     6     9    10     9     6
Sum = 260

	a11   a12  a13   a14   a15
	a21   a22  a23   a24   a25
	a31   a32  a33   a34   a35
	a41   a42  a43   a44   a45
	a51   a52  a53   a54   a55
*/

static void gaussian_blur(int w, int h,
                          unsigned int **dst,
                          unsigned int **src)
{
    int i, j, l, i32opt, i32filtker[3][3], i32Sum;
    unsigned int filt[5][5];

	
	for (i = 0; i < 2; i++){
		for (j = 0; j < w; j++){
			dst[i][j] = src[i][j];
			//printf ("%d", dst[i][j]);
		}
		//printf ("\n");
	}
	printf("\t\t\tEnter the value of Sigma :\n");
	printf("\t\t\t 1> 1.0 :\n");
	printf("\t\t\t 2> 1.2 :\n");
	printf("\t\t\t 3> 1.4 :\n");
	printf("\t\t\t 4> 1.6 :\n");
	printf("\t\t\t 5> 1.8 :\n");
	printf("\t\t\t 6> 2.0 :\n");
	printf("\t\t\t Enter one of the options, [1/2/3/4/5/6]: ");
	scanf("%d", &i32opt);

	switch(i32opt)
	{
	case 1: 
		i32filtker[0][0] = 1;
		i32filtker[0][1] = 5;
		i32filtker[0][2] = 8;

		i32filtker[1][0] = 5;
		i32filtker[1][1] = 23;
		i32filtker[1][2] = 38;

		i32filtker[2][0] = 8;
		i32filtker[2][1] = 38;
		i32filtker[2][2] = 63;
		
		i32Sum = 383;
		break;
	case 2:
		i32filtker[0][0] = 2;
		i32filtker[0][1] = 5;
		i32filtker[0][2] = 8;

		i32filtker[1][0] = 5;
		i32filtker[1][1] = 15;
		i32filtker[1][2] = 21;

		i32filtker[2][0] = 8;
		i32filtker[2][1] = 21;
		i32filtker[2][2] = 30;
		
		i32Sum = 254;
		break;
	case 3:
		i32filtker[0][0] = 2;
		i32filtker[0][1] = 4;
		i32filtker[0][2] = 5;

		i32filtker[1][0] = 4;
		i32filtker[1][1] = 9;
		i32filtker[1][2] = 12;

		i32filtker[2][0] = 5;
		i32filtker[2][1] = 12;
		i32filtker[2][2] = 15;
		
		i32Sum = 159;
		break;
	case 4:
		i32filtker[0][0] = 4;
		i32filtker[0][1] = 8;
		i32filtker[0][2] = 9;

		i32filtker[1][0] = 8;
		i32filtker[1][1] = 14;
		i32filtker[1][2] = 17;

		i32filtker[2][0] = 9;
		i32filtker[2][1] = 17;
		i32filtker[2][2] = 20;	

		i32Sum = 260;
		break;
	case 5:
		i32filtker[0][0] = 5;
		i32filtker[0][1] = 8;
		i32filtker[0][2] = 10;

		i32filtker[1][0] = 8;
		i32filtker[1][1] = 13;
		i32filtker[1][2] = 15;

		i32filtker[2][0] = 10;
		i32filtker[2][1] = 15;
		i32filtker[2][2] = 18;

		i32Sum = 254;
		break;
	case 6:
		i32filtker[0][0] = 6;
		i32filtker[0][1] = 9;
		i32filtker[0][2] = 10;

		i32filtker[1][0] = 9;
		i32filtker[1][1] = 13;
		i32filtker[1][2] = 14;

		i32filtker[2][0] = 10;
		i32filtker[2][1] = 14;
		i32filtker[2][2] = 16;

		i32Sum = 260;
		break;
	default:
		i32filtker[0][0] = 2;
		i32filtker[0][1] = 4;
		i32filtker[0][2] = 5;

		i32filtker[1][0] = 4;
		i32filtker[1][1] = 9;
		i32filtker[1][2] = 12;

		i32filtker[2][0] = 5;
		i32filtker[2][1] = 12;
		i32filtker[2][2] = 15;
		
		i32Sum = 159;
		break;
	}
	for (i = 2; i < h-2; i++){
		dst[i][0] = src[i][0];
		dst[i][1] = src[i][1];
		for (j = 2; j < w-2; j++){
			for (l=0; l<5; l++)
			{
				filt[0][l] = src[i-2][j-2+l];
				filt[1][l] = src[i-1][j-2+l];
				filt[2][l] = src[i][j-2+l];
				filt[3][l] = src[i+1][j-2+l];
				filt[4][l] = src[i+2][j-2+l];
			}

			
// Sigma 1.4
/*
dst[i][j] = ((filt[0][0] + filt[4][0]) * 2
           + (filt[0][1] + filt[4][1]) * 4
           + (filt[0][2] + filt[4][2]) * 5
           + (filt[0][3] + filt[4][3]) * 4
           + (filt[0][4] + filt[4][4]) * 2

           + (filt[1][0] + filt[3][0]) *  4
           + (filt[1][1] + filt[3][1]) *  9
           + (filt[1][2] + filt[3][2]) * 12
           + (filt[1][3] + filt[3][3]) *  9
           + (filt[1][4] + filt[3][4]) *  4

                    + filt[2][0] *  5
                    + filt[2][1] * 12
                    + filt[2][2] * 15
                    + filt[2][3] * 12
                    + filt[2][4] *  5) / 159;			

*/  
			// PSF = fspecial('gaussian',5,1.0);  159 * PSF

dst[i][j] = ((filt[0][0] + filt[4][0]) * i32filtker[0][0]
           + (filt[0][1] + filt[4][1]) * i32filtker[0][1]
           + (filt[0][2] + filt[4][2]) * i32filtker[0][2]
           + (filt[0][3] + filt[4][3]) * i32filtker[0][1]
           + (filt[0][4] + filt[4][4]) * i32filtker[0][0]

           + (filt[1][0] + filt[3][0]) *  i32filtker[1][0]
           + (filt[1][1] + filt[3][1]) *  i32filtker[1][1]
           + (filt[1][2] + filt[3][2]) *  i32filtker[1][2]
           + (filt[1][3] + filt[3][3]) *  i32filtker[1][1]
           + (filt[1][4] + filt[3][4]) *  i32filtker[1][0]

                    + filt[2][0] *  i32filtker[2][0]
                    + filt[2][1] *  i32filtker[2][1]
                    + filt[2][2] *  i32filtker[2][2]
                    + filt[2][3] *  i32filtker[2][1]
                    + filt[2][4] *  i32filtker[2][0]) / i32Sum;	
			//dst[i][j] = src[i][j];
			//printf ("%d", dst[i][j]);
		}
		//printf ("\n");

		dst[i][j]   = src[i][j];
		dst[i][j+1] = src[i][j+1];
	}

	for (l = i; l < h; l++){
		for (j = 0; j < w; j++){
			dst[l][j] = src[l][j];
			//printf ("%d", dst[i][j]);
		}
		//printf ("\n");
	}
}

static void double_threshold(int low, int high, int w, int h,
                                   unsigned int **dst,
                            	unsigned int **src)
{
    int i, j;

    for (j = 0; j < h; j++) {
        for (i = 0; i < w; i++) {
            if (src[j][i] > high) {
                dst[j][i] = src[j][i];
                continue;
            }
			

            if ((!i || i == w - 1 || !j || j == h - 1) &&
                src[j][i] > low &&
                (src[j-1][i-1] > high ||
                 src[j-1][i  ] > high ||

                 src[j-1][i+1] > high ||
                 src[j  ][i-1] > high ||
                 src[j  ][i+1] > high ||

                 src[j+1][i-1] > high ||
                 src[j+1][i  ] > high ||
                 src[j+1][i+1] > high))
                dst[j][i] = src[j][i];
            else
                dst[j][i] = 0;

			
        }

    }
}

static int filter_frame(unsigned int **pu8charDir, unsigned int **pu8charGrad, unsigned int **pu8char, int i32Width, int i32Height)
{
	FILE* fp1;
	int p, i, j, u32low=0, u32high=0;

	fp1 = fopen("F:\Info_R1.bin", "w");
	if (fp1 ==NULL){
		printf ("File Open 2 failed.\n");
		return -1;
	}else
		printf ("\t\t\t File Open 2 success.\n");

		printf("=============================================================\n");
		printf("..................1st Stage --  Gaussian Filtering -----\n\n");

		/* gaussian filter to reduce noise  */
		gaussian_blur(i32Width, i32Height,
					  pu8charGrad,
					  pu8char);

		printf("=============================================================\n");
		printf("..................2nd Stage --  Gradient & Direction using Sobel Gradient Filtering -----\n\n");
		printf("\n");
		printf("\n");
        /* compute the 16-bits gradients and directions for the next step */
        sobel(i32Width, i32Height,
              pu8char,
              pu8charDir,
              pu8charGrad,    i32Width);

        /* non_maximum_suppression() will actually keep & clip what's necessary and
         * ignore the rest, so we need a clean output buffer */        
		printf("=============================================================\n");
		printf("..................3rd Stage --  Non-Maximal Suppression Stage --\n\n");

		for (j = 0; j < i32Height; j++)
		{
			for (i =0; i < i32Width; i++)
			pu8charGrad[j][i] = 0;

		}
        non_maximum_suppression(i32Width, i32Height,
                                pu8charGrad,
                                pu8charDir,
                                pu8char);
		printf("\n");
		printf("\n");
		printf("-------------------------------------------------------------\n");
        /* keep high values, or low values surrounded by high values */
		printf("=============================================================\n");
		printf("..................4th Stage --  Hystersis Thresholding Stage --\n\n");
		printf("Enter the Lower Threshold, L [default:0]: ");
		scanf("%d", &u32low);
		printf("Enter the Higher Threshold, H [default:255]: ");
		scanf("%d", &u32high);
		printf("\n");
		printf("\n");
        double_threshold(u32low, u32high,
                         i32Width, i32Height,
                         pu8char,
                         pu8charGrad);

		printf("=============================================================\n");

		/*
			for (j = 0; j < i32Height; j++){
		for (i = 0; i < i32Width; i++){
			
			printf( "%d ", pu8charDir[j][i] );				
		}
		printf ( "%c", '\n' );
	}*/

	for (j = 0; j < i32Height; j++){
		for (i = 0; i < i32Width; i++){
			
			fprintf( fp1, "%d ", pu8char[j][i] );				
		}
		fprintf ( fp1, "%c", '\n' );
	}
	fclose(fp1);

	printf("Output Image buffer successfully written to :\n");
    printf("F:\Info_R1.bin\n\n");
}

void main()
{
	FILE*            fp;
	int              i32Height;
	int              i32Width, i32val, i32valy, u8val;
	unsigned int     **pu8char=NULL;
	unsigned int     **pu8charGrad=NULL;
	unsigned int     **pu8charDir =NULL;
		
	fp  = fopen("F:\Info_R.bin", "rb");

	if (fp ==NULL){
		printf ("File Open 1 failed.\n");
		return -1;
	}else
		printf ("File Open 1 success.\n");
	
	// Set the Width & Height of the Image Buffer
	i32Width  = IWidth;
	i32Height = IHeight;

	if ((pu8char = (unsigned char* )malloc( i32Height * sizeof( unsigned char* ))) == NULL){
		printf ("Allocation Failed for Frame\n");
		return -1;
	}

	for (i32val = 0; i32val < i32Height; i32val++){	
	    pu8char[i32val] = malloc( i32Width * sizeof (*pu8char[i32val]) );
		memset( pu8char[i32val], 0x00, (i32Width * sizeof (*pu8char[i32val])));
	}
	
	if ((pu8charGrad = (unsigned char* )malloc( i32Height * sizeof( unsigned char* ))) == NULL){
		printf ("Allocation Failed for Gradient\n");
		return -1;
	}

	for (i32val = 0; i32val < i32Height; i32val++){	
	    pu8charGrad[i32val] = malloc( i32Width * sizeof (*pu8charGrad[i32val]) );
		memset( pu8charGrad[i32val], 0x00, (i32Width * sizeof (*pu8charGrad[i32val])));
	}

	if ((pu8charDir = (unsigned char* )malloc( i32Height * sizeof( unsigned char* ))) == NULL){
		printf ("Allocation Failed for Direction\n");
		return -1;
	}

	for (i32val = 0; i32val < i32Height; i32val++){	
	    pu8charDir[i32val] = malloc( i32Width * sizeof (*pu8charDir[i32val]) );
		memset( pu8charDir[i32val], 0x00, (i32Width * sizeof (*pu8charDir[i32val])));
	}

	for (i32val = 0; i32val < i32Height; i32val++){	
		for (i32valy = 0; i32valy < i32Width; i32valy++)
			fscanf( fp, "%d", &pu8char[i32val][i32valy] );				
	}
   
	fseek(fp, 0L, SEEK_SET);

	for (i32val = 0; i32val < i32Height; i32val++){	
		for (i32valy = 0; i32valy < i32Width; i32valy++)
			fscanf( fp, "%d", &pu8charGrad[i32val][i32valy] );				
	}

	fseek(fp, 0L, SEEK_SET);

	for (i32val = 0; i32val < i32Height; i32val++){	
		for (i32valy = 0; i32valy < i32Width; i32valy++)
			fscanf( fp, "%d", &pu8charDir[i32val][i32valy] );				
	}

	filter_frame(pu8charDir, pu8charGrad, pu8char, i32Width, i32Height);
	//DCTCalc(pu8char, i32Width);
}
