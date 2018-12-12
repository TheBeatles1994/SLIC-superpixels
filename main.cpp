#include <opencv2/highgui/highgui.hpp>
#include <opencv/cv.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "SLIC.h"

using namespace cv;
using namespace std;

int main(int argc, char *argv[])
{

    Mat imgMat = imread("imgs/dog.jpg");
    //cout << imgMat.type()<<endl;
#if 1
    {
        unsigned int* img = (unsigned int *)malloc(sizeof(unsigned int) * imgMat.cols * imgMat.rows);
        for (int row = 0; row < imgMat.rows; row++)
        {
            for (int col = 0; col < imgMat.cols; col++)
            {
                //链接：https://www.cnblogs.com/klitech/p/6020491.html
                /* 注意 Mat::at 函数是个模板函数, 需要指明参数类型, 因为这张图是具有红蓝绿三通道的图,
                   所以它的参数类型可以传递一个 Vec3b, 这是一个存放 3 个 uchar 数据的 Vec(向量). 这里
                   提供了索引重载, [2]表示的是返回第三个通道, 在这里是 Red 通道, 第一个通道(Blue)用[0]返回 */
                img[row*imgMat.cols + col] = 0;
                img[row*imgMat.cols + col] |= ( (uchar)255                                 );
                img[row*imgMat.cols + col] |= ( (uchar)imgMat.at<Vec3b>(row, col)[2] <<  8 );   //Red
                img[row*imgMat.cols + col] |= ( (uchar)imgMat.at<Vec3b>(row, col)[1] << 16 );   //Green
                img[row*imgMat.cols + col] |= ( (uchar)imgMat.at<Vec3b>(row, col)[0] << 24 );   //Blue
            }
        }

        int* labels = new int[imgMat.rows * imgMat.cols];
        int numlabels(0);
        SLIC slic;
        slic.DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(img, imgMat.cols, imgMat.rows, labels, numlabels, 400, 40);
        slic.DrawContoursAroundSegments(img, labels, imgMat.cols, imgMat.rows, 0);
        if(labels) delete [] labels;

        Mat slicMat(imgMat.rows,imgMat.cols,CV_8UC3,Scalar(0,0,0));
        for (int row = 0; row < slicMat.rows; row++)
        {
            for (int col = 0; col < slicMat.cols; col++)
            {
                unsigned int x = img[row*slicMat.cols + col];
                slicMat.at<Vec3b>(row, col)[2] = ((char *)&x)[1];   //Red
                slicMat.at<Vec3b>(row, col)[1] = ((char *)&x)[2];   //Green
                slicMat.at<Vec3b>(row, col)[0] = ((char *)&x)[3];   //Blue
            }
        }

        imshow("SLIC", slicMat);
        cvWaitKey();
    }
#endif
//    imshow("Dog Modified", imgMat);
//    cvWaitKey();
    return 0;
}
