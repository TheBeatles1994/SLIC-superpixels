#include <opencv2/highgui/highgui.hpp>
#include <opencv/cv.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
 #include <iomanip> //用于16进制输出
#include "SLIC.h"
#include "glcm.h"

using namespace cv;
using namespace std;

void testSLIC(Mat imgMat);
void testGLCM(Mat imgMat);
void testEndian();
void testBitOper();
string doubleToString(double num);

int main(int argc, char *argv[])
{
    Mat imgMat = imread("imgs/sigmoid.jpg");
    //Mat imgMat = imread("imgs/FB043_293.jpg");
    //Mat imgMat = imread("imgs/dog.jpg");

    testSLIC(imgMat);
    //testGLCM(imgMat);
    //testBitOper();

    return 0;
}
/* ===================================================================
 * @函数功能: 测试SLIC超像素分割
 * @输入参数: 无
 * @输出参数: 无
 * @注意事项:
 *      无
   ===================================================================
 */
void testSLIC(Mat imgMat)
{
    double lambda=1;
    bool newALGO = true;
#if 0
    for(int i=0;i<10;i++)
    {
        cout<<"Processing..."<<endl;
        SLIC slic;
        slic.setLambda(lambda);
        slic.runSLIC(imgMat, 25, 10,newALGO);

        Mat slicMat = slic.getSlicImg();
        //imshow("SLIC", slicMat);
        //cvWaitKey();
        imwrite("slic_new_" + doubleToString(lambda) + ".jpg", slicMat);

        lambda -= 0.001;
    }
#else
    for(int i=0;i<2;i++)
    {
        SLIC slic;
        slic.setLambda(lambda);
        slic.runSLIC(imgMat, 27, 10,newALGO);

        Mat slicMat = slic.getSlicImg();
        //imshow("SLIC", slicMat);
        //cvWaitKey();
        if(newALGO)
            imwrite("slic_new_" + doubleToString(lambda) + ".jpg", slicMat);
        else
            imwrite("slic_old.jpg", slicMat);
        newALGO = !newALGO;
    }

#endif
    cout<<"Finished!"<<endl;
}

void testGLCM(Mat imgMat)
{
    GLCM glcm;
    glcm.run(imgMat);

    TextureEValues EValues = glcm.getEValues();
    Mat imgEnergy = glcm.getEnergy();
    Mat imgContrast = glcm.getContrast();
    Mat imgHomogenity = glcm.getHomogenity();
    Mat imgEntropy = glcm.getEntropy();

    cout<<"Image's Texture EValues:"<<endl;
    cout<<"    Contrast: "<<EValues.contrast<<endl;
    cout<<"    Energy: "<<EValues.energy<<endl;
    cout<<"    EntropyData: "<<EValues.entropy<<endl;
    cout<<"    Homogenity: "<<EValues.homogenity<<endl;
#if 1
    imshow("Energy", imgEnergy);
    imshow("Contrast", imgContrast);
    imshow("Homogenity", imgHomogenity);
    imshow("Entropy", imgEntropy);
    cvWaitKey(0);
#endif
#if 0
    imwrite("Energy1.jpg", imgEnergy);
    imwrite("Contrast1.jpg", imgContrast);
    imwrite("Homogenity1.jpg", imgHomogenity);
    imwrite("Entropy1.jpg", imgEntropy);
#endif
}

void testEndian()
{
    union Test
    {
        int a;
        char b;
    };

    //本机器测试结果（教研室电脑）是小端模式
    Test t;
    t.a = 1;
    if (t.b == 1)
        cout << "Little Endian" << endl;
    else
        cout << "Big Endian" << endl;
}

void testBitOper()
{
    int i = 0;
    unsigned char x = 3, y = 11;
    unsigned char a,b;

    i |= x;
    i |= y<<8;

    cout << hex <<i <<endl;

    a = i;
    b = i/256;
    cout <<(int)a<<" "<<(int)b<<endl;
}

string doubleToString(double num)
{
    char str[256];
    sprintf(str, "%lf", num);
    string result = str;
    return result;
}
