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

int main(int argc, char *argv[])
{

    Mat imgMat = imread("imgs/dog.jpg");

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
    int* labels = new int[imgMat.rows * imgMat.cols];
    int numlabels(0);
    SLIC slic;
    slic.DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(imgMat, labels, numlabels, 400, 40);
    slic.DrawContoursAroundSegments(imgMat, labels, imgMat.cols, imgMat.rows, 0);
    if(labels) delete [] labels;

    imshow("SLIC", imgMat);
    cvWaitKey();
//    imwrite("slic.jpg", slicMat);
}

void testGLCM(Mat imgMat)
{
    GLCM glcm;
    TextureEValues EValues;

    // 纹理特征值矩阵
    // the Matrixs of Texture Features
    Mat imgEnergy, imgContrast, imgHomogenity, imgEntropy;

    Mat dstChannel;
    glcm.getOneChannel(imgMat, dstChannel, CHANNEL_B);

    // 灰度量化
    // Magnitude Gray Image
    glcm.GrayMagnitude(dstChannel, dstChannel, GRAY_8);

    // 计算整幅图像的纹理特征值图像
    // Calculate Texture Features of the whole Image
    glcm.CalcuTextureImages(dstChannel, imgEnergy, imgContrast, imgHomogenity, imgEntropy, 5, GRAY_8, true);
    glcm.CalcuTextureEValue(dstChannel, EValues, 5, GRAY_8);

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
#else
    imwrite("Energy1.jpg", imgEnergy);
    imwrite("Contrast1.jpg", imgContrast);
    imwrite("Homogenity1.jpg", imgHomogenity);
    imwrite("Entropy1.jpg", imgEntropy);
    cvWaitKey(0);
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
