#include <stdio.h>
#include <stdlib.h>
#include <iostream>

int main(){
    FILE *file;                                              //Define file using the FILE data type, * as it is using a pointer for some reason?
    file = fopen("test.imd","rb");                           //open .imd  file is a streamwrited in read mode, b to say it isn't a text file
    float m_s32data;
    if (file){
        unsigned short int VER,W,H;
        // read file version
        // read image width
        // read image height
        fread(&VER,sizeof(unsigned short int),1,file);
        fread(&W,sizeof(unsigned short int),1,file);
        fread(&H,sizeof(unsigned short int),1,file);
        m_s32data = new signed int[W*H];
        // retrieve the 32-bits signed int buffer
        fread(m_s32data,sizeof(signed int),W*H,file);
        fclose(file);
        // If you want to create and use a floating buffer with the floating pixel value,
        // proceed like this
        float *fdata = new float[W*H];

        for (int k=0;k<H*W;k++);
            fdata[k] = m_s32data[k]/1000;
            delete fdata;
            
            delete m_s32data;
    }

    return 0;
}

