// SLIC.cpp: implementation of the SLIC class.
//
// Copyright (C) Radhakrishna Achanta 2012
// All rights reserved
// Email: firstname.lastname@epfl.ch
//////////////////////////////////////////////////////////////////////
//#include "stdafx.h"
#include <cfloat>
#include <cmath>
#include <iostream>
#include <fstream>
#include "SLIC.h"
#include "glcm.h"

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

SLIC::SLIC()
{
    m_lvec = NULL;
    m_avec = NULL;
    m_bvec = NULL;

    m_lvecvec = NULL;
    m_avecvec = NULL;
    m_bvecvec = NULL;
}

SLIC::~SLIC()
{
    if(m_lvec) delete [] m_lvec;
    if(m_avec) delete [] m_avec;
    if(m_bvec) delete [] m_bvec;


    if(m_lvecvec)
    {
        for( int d = 0; d < m_depth; d++ ) delete [] m_lvecvec[d];
        delete [] m_lvecvec;
    }
    if(m_avecvec)
    {
        for( int d = 0; d < m_depth; d++ ) delete [] m_avecvec[d];
        delete [] m_avecvec;
    }
    if(m_bvecvec)
    {
        for( int d = 0; d < m_depth; d++ ) delete [] m_bvecvec[d];
        delete [] m_bvecvec;
    }
}

void SLIC::runSLIC(Mat imgMat, int K, double compactness, bool newALGO)
{
    int* labels = new int[imgMat.rows * imgMat.cols];
    int numlabels(0);

    glcm.run(imgMat);

    DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(imgMat, labels, numlabels, K, compactness, newALGO);
    //slicImg.create(imgMat.size(), imgMat.type());
    imgMat.copyTo(slicImg);
    DrawContoursAroundSegments(slicImg, labels, imgMat.cols, imgMat.rows, 0);
    if(labels) delete [] labels;
}

//==============================================================================
///	RGB2XYZ
///
/// sRGB (D65 illuninant assumption) to XYZ conversion
//==============================================================================
void SLIC::RGB2XYZ(
        const int&		sR,
        const int&		sG,
        const int&		sB,
        double&			X,
        double&			Y,
        double&			Z)
{
    double R = sR/255.0;
    double G = sG/255.0;
    double B = sB/255.0;

    double r, g, b;

    if(R <= 0.04045)	r = R/12.92;
    else				r = pow((R+0.055)/1.055,2.4);
    if(G <= 0.04045)	g = G/12.92;
    else				g = pow((G+0.055)/1.055,2.4);
    if(B <= 0.04045)	b = B/12.92;
    else				b = pow((B+0.055)/1.055,2.4);

    X = r*0.4124564 + g*0.3575761 + b*0.1804375;
    Y = r*0.2126729 + g*0.7151522 + b*0.0721750;
    Z = r*0.0193339 + g*0.1191920 + b*0.9503041;
}

//===========================================================================
///	RGB2LAB
//===========================================================================
void SLIC::RGB2LAB(const int& sR, const int& sG, const int& sB, double& lval, double& aval, double& bval)
{
    //------------------------
    // sRGB to XYZ conversion
    //------------------------
    double X, Y, Z;
    RGB2XYZ(sR, sG, sB, X, Y, Z);

    //------------------------
    // XYZ to LAB conversion
    //------------------------
    double epsilon = 0.008856;	//actual CIE standard
    double kappa   = 903.3;		//actual CIE standard

    double Xr = 0.950456;	//reference white
    double Yr = 1.0;		//reference white
    double Zr = 1.088754;	//reference white

    double xr = X/Xr;
    double yr = Y/Yr;
    double zr = Z/Zr;

    double fx, fy, fz;
    if(xr > epsilon)	fx = pow(xr, 1.0/3.0);
    else				fx = (kappa*xr + 16.0)/116.0;
    if(yr > epsilon)	fy = pow(yr, 1.0/3.0);
    else				fy = (kappa*yr + 16.0)/116.0;
    if(zr > epsilon)	fz = pow(zr, 1.0/3.0);
    else				fz = (kappa*zr + 16.0)/116.0;

    lval = 116.0*fy-16.0;
    aval = 500.0*(fx-fy);
    bval = 200.0*(fy-fz);
}

//===========================================================================
///	DoRGBtoLABConversion
///
///	For whole image: overlaoded floating point version
//===========================================================================
void SLIC::DoRGBtoLABConversion(
        Mat                         imgMat,
        double*&					lvec,
        double*&					avec,
        double*&					bvec)
{
    int sz = m_width*m_height;
    lvec = new double[sz];
    avec = new double[sz];
    bvec = new double[sz];


    for (int row = 0; row < imgMat.rows; row++)
    {
        for (int col = 0; col < imgMat.cols; col++)
        {
            int r = (uchar)imgMat.at<Vec3b>(row, col)[2];
            int g = (uchar)imgMat.at<Vec3b>(row, col)[1];
            int b = (uchar)imgMat.at<Vec3b>(row, col)[0];

            int j = row*m_width + col;
            RGB2LAB( r, g, b, lvec[j], avec[j], bvec[j] );
        }
    }
}
//===========================================================================
///	DoRGBtoLABConversion
///
/// For whole volume
//===========================================================================
void SLIC::DoRGBtoLABConversion(
        unsigned int**&		ubuff,
        double**&					lvec,
        double**&					avec,
        double**&					bvec)
{
    int sz = m_width*m_height;
    for( int d = 0; d < m_depth; d++ )
    {
        for( int j = 0; j < sz; j++ )
        {
            int r = (ubuff[d][j] >> 16) & 0xFF;
            int g = (ubuff[d][j] >>  8) & 0xFF;
            int b = (ubuff[d][j]      ) & 0xFF;

            RGB2LAB( r, g, b, lvec[d][j], avec[d][j], bvec[d][j] );
        }
    }
}

//=================================================================================
/// DrawContoursAroundSegments
/// 内部边界
/// Internal contour drawing option exists. One only needs to comment the if
/// statement inside the loop that looks at neighbourhood.
//=================================================================================
void SLIC::DrawContoursAroundSegments(
        Mat                     imgMat,
        int*&					labels,
        const int&				width,
        const int&				height,
        const unsigned int&				color )
{
    // 八邻域
    const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
    const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

    /*	int sz = width*height;

    vector<bool> istaken(sz, false); // 若当前像素位置已经证明是边界，则此处设置为true

    int mainindex(0); //正在处理的像素的编号

    // 循环处理每一个像素点
    for( int j = 0; j < height; j++ )
    {
        for( int k = 0; k < width; k++ )
        {
            int np(0);
            // 处理该像素点的八邻域像素
            for( int i = 0; i < 8; i++ )
            {
                int x = k + dx8[i];
                int y = j + dy8[i];
                // 判断邻域像素位置是否合法
                if( (x >= 0 && x < width) && (y >= 0 && y < height) )
                {
                    int index = y*width + x;

                    if( false == istaken[index] )//comment this to obtain internal contours
                    {
                        // 若当前处理的像素与其邻域像素不在同一个超像素块中，那么可知当前处理的像素很有可能是边界
                        if( labels[mainindex] != labels[index] ) np++;
                    }
                }
            }
            if( np > 1 )//change to 2 or 3 for thinner lines
            {
                ubuff[mainindex] = color;
                istaken[mainindex] = true;
            }
            mainindex++;
        }
    }*/


    int sz = width*height;
    vector<bool> istaken(sz, false);    //记录该点是否已标记为边缘点
    vector<int> contourx(sz);   //记录x坐标
    vector<int> contoury(sz);   //记录y坐标
    int mainindex(0);int cind(0);
    for( int j = 0; j < height; j++ )
    {
        for( int k = 0; k < width; k++ )
        {
            int np(0);
            //八邻域
            for( int i = 0; i < 8; i++ )
            {
                int x = k + dx8[i];
                int y = j + dy8[i];
                //判断位置是否有效
                if( (x >= 0 && x < width) && (y >= 0 && y < height) )
                {
                    int index = y*width + x;

                    //if( false == istaken[index] )//comment this to obtain internal contours
                    {
                        if( labels[mainindex] != labels[index] ) np++;  //若邻域像素的label与本像素的label不同，则np++
                    }
                }
            }
            if( np > 1 )    //若邻域像素的label与本像素的label不同数量超过了1，则记录进contour数组，同时对应的istaken的本像素设置为true
            {
                contourx[cind] = k;
                contoury[cind] = j;
                istaken[mainindex] = true;
                //img[mainindex] = color;
                cind++;
            }
            mainindex++;
        }
    }

    //边缘像素个数
    int numboundpix = cind;//int(contourx.size());
    for( int j = 0; j < numboundpix; j++ )
    {
        imgMat.at<Vec3b>(contoury[j], contourx[j])[2] = 255;   //Red
        imgMat.at<Vec3b>(contoury[j], contourx[j])[1] = 255;   //Green
        imgMat.at<Vec3b>(contoury[j], contourx[j])[0] = 255;   //Blue

        //在白色边缘像素的邻域加上黑色框
        for( int n = 0; n < 8; n++ )
        {
            int x = contourx[j] + dx8[n];
            int y = contoury[j] + dy8[n];
            //判断位置是否有效
            if( (x >= 0 && x < width) && (y >= 0 && y < height) )
            {
                int ind = y*width + x;
                if(!istaken[ind])
                {
                    imgMat.at<Vec3b>(y, x)[2] = 0;   //Red
                    imgMat.at<Vec3b>(y, x)[1] = 0;   //Green
                    imgMat.at<Vec3b>(y, x)[0] = 0;   //Blue
                }
            }
        }
    }
}


//==============================================================================
///	DetectLabEdges
//==============================================================================
void SLIC::DetectLabEdges(
        const double*				lvec,
        const double*				avec,
        const double*				bvec,
        const int&					width,
        const int&					height,
        vector<double>&				edges)
{
    int sz = width*height;

    edges.resize(sz,0);
    for( int j = 1; j < height-1; j++ )
    {
        for( int k = 1; k < width-1; k++ )
        {
            int i = j*width+k;

            double dx = (lvec[i-1]-lvec[i+1])*(lvec[i-1]-lvec[i+1]) +
                    (avec[i-1]-avec[i+1])*(avec[i-1]-avec[i+1]) +
                    (bvec[i-1]-bvec[i+1])*(bvec[i-1]-bvec[i+1]);

            double dy = (lvec[i-width]-lvec[i+width])*(lvec[i-width]-lvec[i+width]) +
                    (avec[i-width]-avec[i+width])*(avec[i-width]-avec[i+width]) +
                    (bvec[i-width]-bvec[i+width])*(bvec[i-width]-bvec[i+width]);

            //edges[i] = fabs(dx) + fabs(dy);
            edges[i] = dx*dx + dy*dy;
        }
    }
}

//===========================================================================
///	PerturbSeeds
//===========================================================================
void SLIC::PerturbSeeds(
        vector<double>&				kseedsl,
        vector<double>&				kseedsa,
        vector<double>&				kseedsb,
        vector<double>&				kseedsx,
        vector<double>&				kseedsy,
        const vector<double>&                   edges)
{
    const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
    const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

    int numseeds = kseedsl.size();

    for( int n = 0; n < numseeds; n++ )
    {
        int ox = kseedsx[n];//original x
        int oy = kseedsy[n];//original y
        int oind = oy*m_width + ox;

        int storeind = oind;
        for( int i = 0; i < 8; i++ )
        {
            int nx = ox+dx8[i];//new x
            int ny = oy+dy8[i];//new y

            if( nx >= 0 && nx < m_width && ny >= 0 && ny < m_height)
            {
                int nind = ny*m_width + nx;
                if( edges[nind] < edges[storeind])
                {
                    storeind = nind;
                }
            }
        }
        if(storeind != oind)
        {
            kseedsx[n] = storeind%m_width;
            kseedsy[n] = storeind/m_width;
            kseedsl[n] = m_lvec[storeind];
            kseedsa[n] = m_avec[storeind];
            kseedsb[n] = m_bvec[storeind];
        }
    }
}


//===========================================================================
///	GetLABXYSeeds_ForGivenStepSize
///
/// The k seed values are taken as uniform spatial pixel samples.
//===========================================================================
void SLIC::GetLABXYSeeds_ForGivenStepSize(
        vector<double>&				kseedsl,
        vector<double>&				kseedsa,
        vector<double>&				kseedsb,
        vector<double>&				kseedsx,
        vector<double>&				kseedsy,
        const int&					STEP,
        const bool&					perturbseeds,
        const vector<double>&       edgemag)
{
    const bool hexgrid = false;
    int numseeds(0);
    int n(0);

    //int xstrips = m_width/STEP;
    //int ystrips = m_height/STEP;
    int xstrips = (0.5+double(m_width)/double(STEP));
    int ystrips = (0.5+double(m_height)/double(STEP));

    int xerr = m_width  - STEP*xstrips;if(xerr < 0){xstrips--;xerr = m_width - STEP*xstrips;}
    int yerr = m_height - STEP*ystrips;if(yerr < 0){ystrips--;yerr = m_height- STEP*ystrips;}

    double xerrperstrip = double(xerr)/double(xstrips);
    double yerrperstrip = double(yerr)/double(ystrips);

    int xoff = STEP/2;
    int yoff = STEP/2;
    //-------------------------
    numseeds = xstrips*ystrips;
    //-------------------------
    kseedsl.resize(numseeds);
    kseedsa.resize(numseeds);
    kseedsb.resize(numseeds);
    kseedsx.resize(numseeds);
    kseedsy.resize(numseeds);

    for( int y = 0; y < ystrips; y++ )
    {
        int ye = y*yerrperstrip;
        for( int x = 0; x < xstrips; x++ )
        {
            int xe = x*xerrperstrip;
            int seedx = (x*STEP+xoff+xe);
            if(hexgrid){ seedx = x*STEP+(xoff<<(y&0x1))+xe; seedx = min(m_width-1,seedx); }//for hex grid sampling
            int seedy = (y*STEP+yoff+ye);
            int i = seedy*m_width + seedx;

            kseedsl[n] = m_lvec[i];
            kseedsa[n] = m_avec[i];
            kseedsb[n] = m_bvec[i];
            kseedsx[n] = seedx;
            kseedsy[n] = seedy;
            n++;
        }
    }


    if(perturbseeds)
    {
        PerturbSeeds(kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, edgemag);
    }
}

//===========================================================================
///	GetKValues_LABXYZ
///
/// The k seed values are taken as uniform spatial pixel samples.
//===========================================================================
void SLIC::GetKValues_LABXYZ(
        vector<double>&				kseedsl,
        vector<double>&				kseedsa,
        vector<double>&				kseedsb,
        vector<double>&				kseedsx,
        vector<double>&				kseedsy,
        vector<double>&				kseedsz,
        const int&				STEP)
{
    const bool hexgrid = false;
    int numseeds(0);
    int n(0);

    int xstrips = (0.5+double(m_width)/double(STEP));
    int ystrips = (0.5+double(m_height)/double(STEP));
    int zstrips = (0.5+double(m_depth)/double(STEP));

    int xerr = m_width  - STEP*xstrips;if(xerr < 0){xstrips--;xerr = m_width - STEP*xstrips;}
    int yerr = m_height - STEP*ystrips;if(yerr < 0){ystrips--;yerr = m_height- STEP*ystrips;}
    int zerr = m_depth  - STEP*zstrips;if(zerr < 0){zstrips--;zerr = m_depth - STEP*zstrips;}

    double xerrperstrip = double(xerr)/double(xstrips);
    double yerrperstrip = double(yerr)/double(ystrips);
    double zerrperstrip = double(zerr)/double(zstrips);

    int xoff = STEP/2;
    int yoff = STEP/2;
    int zoff = STEP/2;
    //-------------------------
    numseeds = xstrips*ystrips*zstrips;
    //-------------------------
    kseedsl.resize(numseeds);
    kseedsa.resize(numseeds);
    kseedsb.resize(numseeds);
    kseedsx.resize(numseeds);
    kseedsy.resize(numseeds);
    kseedsz.resize(numseeds);

    for( int z = 0; z < zstrips; z++ )
    {
        int ze = z*zerrperstrip;
        int d = (z*STEP+zoff+ze);
        for( int y = 0; y < ystrips; y++ )
        {
            int ye = y*yerrperstrip;
            for( int x = 0; x < xstrips; x++ )
            {
                int xe = x*xerrperstrip;
                int i = (y*STEP+yoff+ye)*m_width + (x*STEP+xoff+xe);

                kseedsl[n] = m_lvecvec[d][i];
                kseedsa[n] = m_avecvec[d][i];
                kseedsb[n] = m_bvecvec[d][i];
                kseedsx[n] = (x*STEP+xoff+xe);
                kseedsy[n] = (y*STEP+yoff+ye);
                kseedsz[n] = d;
                n++;
            }
        }
    }
}

//===========================================================================
///	PerformSuperpixelSLIC
/// 注意SLIC的思路与K-mean很像，只是SLIC算法因为专注于局部，而K-mean是全部，故SLIC要更快
///	Performs k mean segmentation. It is fast because it looks locally, not
/// over the entire image.
//===========================================================================
void SLIC::PerformSuperpixelSLIC(
        vector<double>&				kseedsl,
        vector<double>&				kseedsa,
        vector<double>&				kseedsb,
        vector<double>&				kseedsx,
        vector<double>&				kseedsy,
        int*&					klabels,
        const int&				STEP,
        const vector<double>&                   edgemag,
        const double&				M)
{
    int sz = m_width*m_height;
    const int numk = kseedsl.size();	// 即种子点的个数
    //----------------
    int offset = STEP;
    //if(STEP < 8) offset = STEP*1.5;//to prevent a crash due to a very small step size
    //----------------

    vector<double> clustersize(numk, 0);
    vector<double> inv(numk, 0);//to store 1/clustersize[k] values

    vector<double> sigmal(numk, 0);
    vector<double> sigmaa(numk, 0);
    vector<double> sigmab(numk, 0);
    vector<double> sigmax(numk, 0);
    vector<double> sigmay(numk, 0);
    vector<double> distvec(sz, DBL_MAX);

    double invwt = 1.0/((STEP/M)*(STEP/M));

    int x1, y1, x2, y2;
    double l, a, b;
    double dist;
    double distxy;
    for( int itr = 0; itr < 10; itr++ )
    {
        // 每个像素到其种子点的最短距离
        distvec.assign(sz, DBL_MAX);
        for( int n = 0; n < numk; n++ )	//[0, numk)是种子点的label
        {
            y1 = max(0.0,			kseedsy[n]-offset);
            y2 = min((double)m_height,	kseedsy[n]+offset);
            x1 = max(0.0,			kseedsx[n]-offset);
            x2 = min((double)m_width,	kseedsx[n]+offset);


            // 在2S*2S的范围内的每个像素进行计算，即计算其到种子点的距离
            for( int y = y1; y < y2; y++ )
            {
                for( int x = x1; x < x2; x++ )
                {
                    // i是图像中第几个像素点
                    int i = y*m_width + x;

                    l = m_lvec[i];
                    a = m_avec[i];
                    b = m_bvec[i];

                    dist =			(l - kseedsl[n])*(l - kseedsl[n]) +
                            (a - kseedsa[n])*(a - kseedsa[n]) +
                            (b - kseedsb[n])*(b - kseedsb[n]);

                    distxy =		(x - kseedsx[n])*(x - kseedsx[n]) +
                            (y - kseedsy[n])*(y - kseedsy[n]);

                    //------------------------------------------------------------------------
                    dist += distxy*invwt;//dist = sqrt(dist) + sqrt(distxy*invwt);//this is more exact
                    //------------------------------------------------------------------------
                    if( dist < distvec[i] )
                    {
                        distvec[i] = dist;//表示该像素到种子点的最短距离
                        klabels[i]  = n;//表示该像素所属的种子点label
                    }
                }
            }
        }
        //-----------------------------------------------------------------
        // Recalculate the centroid and store in the seed values
        // 重新计算centroid图心
        //-----------------------------------------------------------------
        //instead of reassigning memory on each iteration, just reset.

        sigmal.assign(numk, 0);
        sigmaa.assign(numk, 0);
        sigmab.assign(numk, 0);
        sigmax.assign(numk, 0);
        sigmay.assign(numk, 0);
        clustersize.assign(numk, 0);
        //------------------------------------
        //edgesum.assign(numk, 0);
        //------------------------------------

        {int ind(0);
            for( int r = 0; r < m_height; r++ )
            {
                for( int c = 0; c < m_width; c++ )
                {
                    // 对所有的l a b x y进行求和
                    sigmal[klabels[ind]] += m_lvec[ind];
                    sigmaa[klabels[ind]] += m_avec[ind];
                    sigmab[klabels[ind]] += m_bvec[ind];
                    sigmax[klabels[ind]] += c;
                    sigmay[klabels[ind]] += r;
                    //------------------------------------
                    //edgesum[klabels[ind]] += edgemag[ind];
                    //------------------------------------
                    clustersize[klabels[ind]] += 1.0;
                    ind++;
                }
            }}

        {for( int k = 0; k < numk; k++ )
            {
                if( clustersize[k] <= 0 ) clustersize[k] = 1;
                inv[k] = 1.0/clustersize[k];//computing inverse now to multiply, than divide later
            }}
        //求平均，即新的l a b x y
        {for( int k = 0; k < numk; k++ )
            {
                kseedsl[k] = sigmal[k]*inv[k];
                kseedsa[k] = sigmaa[k]*inv[k];
                kseedsb[k] = sigmab[k]*inv[k];
                kseedsx[k] = sigmax[k]*inv[k];
                kseedsy[k] = sigmay[k]*inv[k];
                //------------------------------------
                //edgesum[k] *= inv[k];
                //------------------------------------
            }}
        // 重新进行新的迭代过程，一共十次
    }
}
//===========================================================================
void SLIC::newPerformSuperpixelSLIC(
        vector<double>&				kseedsl,
        vector<double>&				kseedsa,
        vector<double>&				kseedsb,
        vector<double>&				kseedsx,
        vector<double>&				kseedsy,
        int*&					klabels,
        const int&				STEP,
        const vector<double>&                   edgemag,
        const double&				M)
{
    TextureEValues EValues = glcm.getEValues();
    Mat imgEnergy = glcm.getEnergy();
    Mat imgContrast = glcm.getContrast();
    Mat imgHomogenity = glcm.getHomogenity();
    Mat imgEntropy = glcm.getEntropy();

    int sz = m_width*m_height;
    const int numk = kseedsl.size();	// 即种子点的个数
    //----------------
    int offset = STEP;
    //if(STEP < 8) offset = STEP*1.5;//to prevent a crash due to a very small step size
    //----------------

    vector<double> clustersize(numk, 0);
    vector<double> inv(numk, 0);//to store 1/clustersize[k] values

    vector<double> sigmal(numk, 0);
    vector<double> sigmaa(numk, 0);
    vector<double> sigmab(numk, 0);
    vector<double> sigmax(numk, 0);
    vector<double> sigmay(numk, 0);
    vector<double> distvec(sz, DBL_MAX);

    double invwt = 1.0/((STEP/M)*(STEP/M));

    int x1, y1, x2, y2;
    double l, a, b;
    double dist;
    double distxy;
    for( int itr = 0; itr < 10; itr++ )
    {
        // 每个像素到其种子点的最短距离
        distvec.assign(sz, DBL_MAX);
        for( int n = 0; n < numk; n++ )	//[0, numk)是种子点的label
        {
            y1 = max(0.0,			kseedsy[n]-offset);
            y2 = min((double)m_height,	kseedsy[n]+offset);
            x1 = max(0.0,			kseedsx[n]-offset);
            x2 = min((double)m_width,	kseedsx[n]+offset);

            // 在2S*2S的范围内的每个像素进行计算，即计算其到种子点的距离
            for( int y = y1; y < y2; y++ )
            {
                for( int x = x1; x < x2; x++ )
                {
                    // i是图像中第几个像素点
                    int i = y*m_width + x;

                    l = m_lvec[i];
                    a = m_avec[i];
                    b = m_bvec[i];

                    dist =			(l - kseedsl[n])*(l - kseedsl[n]) +
                            (a - kseedsa[n])*(a - kseedsa[n]) +
                            (b - kseedsb[n])*(b - kseedsb[n]);

                    distxy =		(x - kseedsx[n])*(x - kseedsx[n]) +
                            (y - kseedsy[n])*(y - kseedsy[n]);

                    //------------------------------------------------------------------------
                    dist += distxy*invwt;//dist = sqrt(dist) + sqrt(distxy*invwt);//this is more exact

                    double distEnergy = (imgEnergy.at<uchar>(kseedsy[n], kseedsx[n]) - imgEnergy.at<uchar>(y, x))*(imgEnergy.at<uchar>(kseedsy[n], kseedsx[n]) - imgEnergy.at<uchar>(y, x));
                    double distContrast = (imgContrast.at<uchar>(kseedsy[n], kseedsx[n]) - imgContrast.at<uchar>(y, x))*(imgContrast.at<uchar>(kseedsy[n], kseedsx[n]) - imgContrast.at<uchar>(y, x));
                    double distHomogenity = (imgHomogenity.at<uchar>(kseedsy[n], kseedsx[n]) - imgHomogenity.at<uchar>(y, x))*(imgHomogenity.at<uchar>(kseedsy[n], kseedsx[n]) - imgHomogenity.at<uchar>(y, x));
                    double distEntropy = (imgEntropy.at<uchar>(kseedsy[n], kseedsx[n]) - imgEntropy.at<uchar>(y, x))*(imgEntropy.at<uchar>(kseedsy[n], kseedsx[n]) - imgEntropy.at<uchar>(y, x));
                    double TextureSum = distEnergy + distContrast + distHomogenity + distEntropy;
                    dist = lambda*dist + (1-lambda)*TextureSum;
                    //------------------------------------------------------------------------
                    if( dist < distvec[i] )
                    {
                        distvec[i] = dist;//表示该像素到种子点的最短距离
                        klabels[i]  = n;//表示该像素所属的种子点label
                    }
                }
            }
        }
        //-----------------------------------------------------------------
        // Recalculate the centroid and store in the seed values
        // 重新计算centroid图心
        //-----------------------------------------------------------------
        //instead of reassigning memory on each iteration, just reset.

        sigmal.assign(numk, 0);
        sigmaa.assign(numk, 0);
        sigmab.assign(numk, 0);
        sigmax.assign(numk, 0);
        sigmay.assign(numk, 0);
        clustersize.assign(numk, 0);
        //------------------------------------
        //edgesum.assign(numk, 0);
        //------------------------------------

        {int ind(0);
            for( int r = 0; r < m_height; r++ )
            {
                for( int c = 0; c < m_width; c++ )
                {
                    // 对所有的l a b x y进行求和
                    sigmal[klabels[ind]] += m_lvec[ind];
                    sigmaa[klabels[ind]] += m_avec[ind];
                    sigmab[klabels[ind]] += m_bvec[ind];
                    sigmax[klabels[ind]] += c;
                    sigmay[klabels[ind]] += r;
                    //------------------------------------
                    //edgesum[klabels[ind]] += edgemag[ind];
                    //------------------------------------
                    clustersize[klabels[ind]] += 1.0;
                    ind++;
                }
            }}

        {for( int k = 0; k < numk; k++ )
            {
                if( clustersize[k] <= 0 ) clustersize[k] = 1;
                inv[k] = 1.0/clustersize[k];//computing inverse now to multiply, than divide later
            }}
        //求平均，即新的l a b x y
        {for( int k = 0; k < numk; k++ )
            {
                kseedsl[k] = sigmal[k]*inv[k];
                kseedsa[k] = sigmaa[k]*inv[k];
                kseedsb[k] = sigmab[k]*inv[k];
                kseedsx[k] = sigmax[k]*inv[k];
                kseedsy[k] = sigmay[k]*inv[k];
                //------------------------------------
                //edgesum[k] *= inv[k];
                //------------------------------------
            }}
        // 重新进行新的迭代过程，一共十次
    }
}

//===========================================================================
///	SaveSuperpixelLabels
///
///	Save labels in raster scan order.
//===========================================================================
void SLIC::SaveSuperpixelLabels(
        const int*&					labels,
        const int&					width,
        const int&					height,
        const string&				filename,
        const string&				path)
{
#ifdef WINDOWS
    char fname[256];
    char extn[256];
    _splitpath(filename.c_str(), NULL, NULL, fname, extn);
    string temp = fname;
    string finalpath = path + temp + string(".dat");
#else
    string nameandextension = filename;
    size_t pos = filename.find_last_of("/");
    if(pos != string::npos)//if a slash is found, then take the filename with extension
    {
        nameandextension = filename.substr(pos+1);
    }
    string newname = nameandextension.replace(nameandextension.rfind(".")+1, 3, "dat");//find the position of the dot and replace the 3 characters following it.
    string finalpath = path+newname;
#endif

    int sz = width*height;
    ofstream outfile;
    outfile.open(finalpath.c_str(), ios::binary);
    for( int i = 0; i < sz; i++ )
    {
        outfile.write((const char*)&labels[i], sizeof(int));
    }
    outfile.close();
}


//===========================================================================
///	SaveSupervoxelLabels
///
///	Save labels in raster scan order.
//===========================================================================
void SLIC::SaveSupervoxelLabels(
        const int**&				labels,
        const int&					width,
        const int&					height,
        const int&					depth,
        const string&				filename,
        const string&				path)
{
#ifdef WINDOWS
    char fname[256];
    char extn[256];
    _splitpath(filename.c_str(), NULL, NULL, fname, extn);
    string temp = fname;
    string finalpath = path + temp + string(".dat");
#else
    string nameandextension = filename;
    size_t pos = filename.find_last_of("/");
    if(pos != string::npos)//if a slash is found, then take the filename with extension
    {
        nameandextension = filename.substr(pos+1);
    }
    string newname = nameandextension.replace(nameandextension.rfind(".")+1, 3, "dat");//find the position of the dot and replace the 3 characters following it.
    string finalpath = path+newname;
#endif

    int sz = width*height;
    ofstream outfile;
    outfile.open(finalpath.c_str(), ios::binary);
    for( int d = 0; d < depth; d++ )
    {
        for( int i = 0; i < sz; i++ )
        {
            outfile.write((const char*)&labels[d][i], sizeof(int));
        }
    }
    outfile.close();
}

//===========================================================================
///	EnforceLabelConnectivity
///
///		1. finding an adjacent label for each new component at the start
///		2. if a certain component is too small, assigning the previously found
///		    adjacent label to this component, and not incrementing the label.
//===========================================================================
void SLIC::EnforceLabelConnectivity(
        const int*					labels,//input labels that need to be corrected to remove stray labels 旧标签
        const int					width,
        const int					height,
        int*&						nlabels,//new labels 新标签
        int&						numlabels,//the number of labels changes in the end if segments are removed 新标签的数量
        const int&					K) //the number of superpixels desired by the user
{
    //	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
    //	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

    const int dx4[4] = {-1,  0,  1,  0};
    const int dy4[4] = { 0, -1,  0,  1};

    const int sz = width*height;
    const int SUPSZ = sz/K;
    //nlabels.resize(sz, -1);
    for( int i = 0; i < sz; i++ ) nlabels[i] = -1;
    int label(0); // 新标签

    // 同一个标签下的像素
    int* xvec = new int[sz];
    int* yvec = new int[sz];

    int oindex(0);	//正在处理的像素的编号
    int adjlabel(0);//adjacent label
    //开始处理每一个像素 按行扫描
    for( int j = 0; j < height; j++ )
    {
        for( int k = 0; k < width; k++ )
        {
            //当该像素还没有被处理过时，开始进行处理
            if( 0 > nlabels[oindex] )
            {
                nlabels[oindex] = label;
                //--------------------
                // Start a new segment 将当前坐标加到新的超像素块坐标集合中
                //--------------------
                xvec[0] = k;
                yvec[0] = j;
                //-------------------------------------------------------
                // Quickly find an adjacent label for use later if needed
                // 先找一个临近的已经标记过的像素的标记号，当本超像素块中的像素个数偏少时会用得上
                //-------------------------------------------------------
                {for( int n = 0; n < 4; n++ )
                    {
                        int x = xvec[0] + dx4[n];
                        int y = yvec[0] + dy4[n];
                        //当临近像素的坐标合法时才进行处理
                        if( (x >= 0 && x < width) && (y >= 0 && y < height) )
                        {
                            int nindex = y*width + x;
                            //若临近像素处理过，则adjlabel值设置为临近像素的值
                            if(nlabels[nindex] >= 0) adjlabel = nlabels[nindex];
                        }
                    }}

                int count(1);//本超像素块中的像素个数
                for( int c = 0; c < count; c++ )
                {
                    //取本像素的四邻域像素
                    for( int n = 0; n < 4; n++ )
                    {
                        int x = xvec[c] + dx4[n];
                        int y = yvec[c] + dy4[n];

                        // 判断邻域像素位置是否合法
                        if( (x >= 0 && x < width) && (y >= 0 && y < height) )
                        {
                            int nindex = y*width + x;

                            // 若邻域像素未被标记过，同时邻域像素与当前处理的像素在旧标记上是同一个分类，那么添加此坐标到超像素坐标集合中
                            if( 0 > nlabels[nindex] && labels[oindex] == labels[nindex] )
                            {
                                xvec[count] = x;
                                yvec[count] = y;
                                nlabels[nindex] = label;
                                count++;
                            }
                        }
                    }
                }
                //-------------------------------------------------------
                // If segment size is less then a limit, assign an
                // adjacent label found before, and decrement label count.
                // 当本超像素块中的像素个数不足规定的一半时，将该超像素块归并到临近的像素块中，同时此次循环label标记不发生变化
                //-------------------------------------------------------
                if(count <= SUPSZ >> 2)
                {
                    for( int c = 0; c < count; c++ )
                    {
                        int ind = yvec[c]*width+xvec[c];
                        nlabels[ind] = adjlabel;
                    }
                    label--;
                }
                label++;
            }
            oindex++;
        }
    }
    numlabels = label;

    if(xvec) delete [] xvec;
    if(yvec) delete [] yvec;
}
/* ===================================================================
 * @函数功能: 新算法 - 连接超像素块
 * @输入参数: 无
 * @输出参数: 无
 * @注意事项:
 *      无
   ===================================================================
 */
void SLIC::newEnforceLabelConnectivity(
        Mat                         imgMat,
        const int*					labels,//input labels that need to be corrected to remove stray labels 旧标签
        const int					width,
        const int					height,
        int*&						nlabels,//new labels 新标签
        int&						numlabels,//the number of labels changes in the end if segments are removed 新标签的数量
        const int&					K) //the number of superpixels desired by the user
{
    //	const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
    //	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};

    const int dx4[4] = {-1,  0,  1,  0};
    const int dy4[4] = { 0, -1,  0,  1};

    const int sz = width*height;
    const int SUPSZ = sz/K;
    //nlabels.resize(sz, -1);
    for( int i = 0; i < sz; i++ ) nlabels[i] = -1;  //新的labels，全部标记为-1
    int label(0); // 新标签，从0开始计数
    //vector<vector<pair<int, int> > > labelCoordinateVec(sz, vector<pair<int, int> >());    //保存每一个标签下面的坐标值pair(x, y)
    vector<vector<pair<int, int> > > labelCoordinateVec;

    // 同一个标签下的像素
    int* xvec = new int[sz];
    int* yvec = new int[sz];

    int oindex(0);	//正在处理的像素的编号
    int adjlabel(0);//adjacent label
    vector<int> adjlabelvec;    //adjacent label vector
    Mat imgEnergy = glcm.getEnergy();
    Mat imgContrast = glcm.getContrast();
    Mat imgHomogenity = glcm.getHomogenity();
    Mat imgEntropy = glcm.getEntropy();
    //开始处理每一个像素 按行扫描
    for( int j = 0; j < height; j++ )
    {
        for( int k = 0; k < width; k++ )
        {
            //当该像素还没有被处理过时，开始进行处理
            //注意此处每次并不是只处理一个像素！而是多个！
            if( 0 > nlabels[oindex] )
            {
                vector<pair<int, int> > tempPair;
                nlabels[oindex] = label;
                //--------------------
                // Start a new segment 将当前坐标加到新的超像素块坐标集合中
                //--------------------
                xvec[0] = k;
                yvec[0] = j;
                tempPair.push_back(pair<int, int>(k, j));

                /* 注意此处将所有的同本像素一样label的像素都遍历了遍，同时可得到本超像素区域的像素总个数count */
                int count(1);//本超像素块中的像素个数
                for( int c = 0; c < count; c++ )
                {
                    //取本像素的四邻域像素
                    for( int n = 0; n < 4; n++ )
                    {
                        int x = xvec[c] + dx4[n];
                        int y = yvec[c] + dy4[n];

                        // 判断邻域像素位置是否合法
                        if( (x >= 0 && x < width) && (y >= 0 && y < height) )
                        {
                            int nindex = y*width + x;

                            // 若邻域像素未被标记过，同时邻域像素与当前处理的像素在旧标记上是同一个分类，那么添加此坐标到超像素坐标集合中，同时对其标记
                            if( 0 > nlabels[nindex] && labels[oindex] == labels[nindex] )
                            {
                                xvec[count] = x;
                                yvec[count] = y;
                                tempPair.push_back(pair<int, int>(x, y));
                                nlabels[nindex] = label;
                                count++;
                            }
                        }
                    }
                }
                //-------------------------------------------------------
                // Quickly find an adjacent label for use later if needed
                // 先对本像素的四邻域附近找一个已经标记过的与本身label不同的像素的标记号，当本超像素块中的像素个数偏少时会用得上
                //-------------------------------------------------------
                {for( int n = 0; n < 4; n++ )
                    {
                        int x = xvec[0] + dx4[n];
                        int y = yvec[0] + dy4[n];
                        //当临近像素的坐标合法时才进行处理
                        if( (x >= 0 && x < width) && (y >= 0 && y < height) )
                        {
                            int nindex = y*width + x;
                            //创新点：在此处设置一个vector，添加不同的邻域像素值，以便后面进行比较，看选择将adjlabel值设置为哪一个临近像素的label
                            if(nlabels[nindex] >= 0 && nlabels[nindex]!=label)
                            {
                                adjlabel = nlabels[nindex];
                                adjlabelvec.push_back(nlabels[nindex]); //注意此处添加的nlabels[nindex]一定是比当前label值小的
                            }
                        }
                    }}

                //-------------------------------------------------------
                // If segment size is less then a limit, assign an
                // adjacent label found before, and decrement label count.
                // 当本超像素块中的像素个数不足规定的一半时，将该超像素块归并到临近的像素块中，同时此次循环label标记不发生变化
                //-------------------------------------------------------
                if(count <= SUPSZ >> 2)
                {
                    struct TextureEValues tValues;
                    tValues.contrast = 0;
                    tValues.energy = 0;
                    tValues.entropy = 0;
                    tValues.homogenity = 0;

                    //统计本超像素块内的纹理信息
                    for(auto it=tempPair.begin();it!=tempPair.end();it++)
                    {
                        tValues.energy += imgEnergy.at<uchar>((*it).second, (*it).first);
                        tValues.contrast += imgContrast.at<uchar>((*it).second, (*it).first);
                        tValues.homogenity += imgHomogenity.at<uchar>((*it).second, (*it).first);
                        tValues.entropy += imgEntropy.at<uchar>((*it).second, (*it).first);
                    }

                    tValues.contrast /= tempPair.size();
                    tValues.energy /= tempPair.size();
                    tValues.entropy /= tempPair.size();
                    tValues.homogenity /= tempPair.size();

                    double dist = 10000000;
                    //比较本超像素块与每个临近超像素块的纹理信息，并进行最终的选择
                    for(auto ita = adjlabelvec.begin();ita!=adjlabelvec.end();ita++)
                    {
                        struct TextureEValues tempValues;
                        tempValues.contrast = 0;
                        tempValues.energy = 0;
                        tempValues.entropy = 0;
                        tempValues.homogenity = 0;
                        for(auto itb = labelCoordinateVec[*ita].begin();itb!=labelCoordinateVec[*ita].end();itb++)
                        {
                            tempValues.energy += imgEnergy.at<uchar>((*itb).second, (*itb).first);
                            tempValues.contrast += imgContrast.at<uchar>((*itb).second, (*itb).first);
                            tempValues.homogenity += imgHomogenity.at<uchar>((*itb).second, (*itb).first);
                            tempValues.entropy += imgEntropy.at<uchar>((*itb).second, (*itb).first);
                        }
                        tempValues.contrast /= labelCoordinateVec[*ita].size();
                        tempValues.energy /= labelCoordinateVec[*ita].size();
                        tempValues.entropy /= labelCoordinateVec[*ita].size();
                        tempValues.homogenity /= labelCoordinateVec[*ita].size();

                        double tempdist = (tValues.energy - tempValues.energy)*(tValues.energy - tempValues.energy)+
                                (tempValues.contrast - tempValues.contrast)*(tempValues.contrast - tempValues.contrast)+
                                (tempValues.homogenity - tempValues.homogenity)*(tempValues.homogenity - tempValues.homogenity)+
                                (tempValues.entropy - tempValues.entropy)*(tempValues.entropy - tempValues.entropy);
                        if(tempdist<dist)
                        {
                            dist = tempdist;
                            adjlabel = *ita;
                        }
                    }

                    for( int c = 0; c < count; c++ )
                    {
                        int ind = yvec[c]*width+xvec[c];
                        nlabels[ind] = adjlabel;
                    }
                    for(auto it=tempPair.begin();it!=tempPair.end();it++)
                    {
                        labelCoordinateVec[adjlabel].push_back(*it);
                    }
                    label--;
                }
                else
                    labelCoordinateVec.push_back(tempPair);
                label++;

                adjlabelvec.clear();
            }
            oindex++;
        }
    }
    numlabels = label;

    if(xvec) delete [] xvec;
    if(yvec) delete [] yvec;
}
//===========================================================================
///	DoSuperpixelSegmentation_ForGivenSuperpixelSize
///
/// The input parameter ubuff conains RGB values in a 32-bit unsigned integers
/// as follows:
///
/// [1 1 1 1 1 1 1 1]  [1 1 1 1 1 1 1 1]  [1 1 1 1 1 1 1 1]  [1 1 1 1 1 1 1 1]
///
///        Nothing              R                 G                  B
///
/// The RGB values are accessed from (and packed into) the unsigned integers
/// using bitwise operators as can be seen in the function DoRGBtoLABConversion().
///
/// compactness value depends on the input pixels values. For instance, if
/// the input is greyscale with values ranging from 0-100, then a compactness
/// value of 20.0 would give good results. A greater value will make the
/// superpixels more compact while a smaller value would make them more uneven.
///
/// The labels can be saved if needed using SaveSuperpixelLabels()
//===========================================================================
void SLIC::DoSuperpixelSegmentation_ForGivenSuperpixelSize(
        Mat                         imgMat,
        int*&						klabels,
        int&						numlabels,
        const int&					superpixelsize,
        const double&               compactness,
        bool                        newALGO)
{
    //------------------------------------------------
    const int STEP = sqrt(double(superpixelsize))+0.5;
    //------------------------------------------------
    vector<double> kseedsl(0);
    vector<double> kseedsa(0);
    vector<double> kseedsb(0);
    vector<double> kseedsx(0);
    vector<double> kseedsy(0);

    //--------------------------------------------------
    m_width  = imgMat.cols;
    m_height = imgMat.rows;
    int sz = m_width*m_height;
    //klabels.resize( sz, -1 );
    //--------------------------------------------------
    klabels = new int[sz];
    for( int s = 0; s < sz; s++ ) klabels[s] = -1;
    //--------------------------------------------------
    if(1)//LAB, the default option
    {
        DoRGBtoLABConversion(imgMat, m_lvec, m_avec, m_bvec);
    }
    else//RGB
    {
        m_lvec = new double[sz]; m_avec = new double[sz]; m_bvec = new double[sz];

        for (int row = 0; row < imgMat.rows; row++)
        {
            for (int col = 0; col < imgMat.cols; col++)
            {
                int i = row*m_width + col;
                m_lvec[i] = (uchar)imgMat.at<Vec3b>(row, col)[2];
                m_avec[i] = (uchar)imgMat.at<Vec3b>(row, col)[1];
                m_bvec[i] = (uchar)imgMat.at<Vec3b>(row, col)[0];
            }
        }
    }
    //--------------------------------------------------
    bool perturbseeds(false);//perturb seeds is not absolutely necessary, one can set this flag to false
    vector<double> edgemag(0);
    if(perturbseeds) DetectLabEdges(m_lvec, m_avec, m_bvec, m_width, m_height, edgemag);
    GetLABXYSeeds_ForGivenStepSize(kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, STEP, perturbseeds, edgemag);

    if(!newALGO)
        PerformSuperpixelSLIC(kseedsl, kseedsa, kseedsb, kseedsx, kseedsy, klabels, STEP, edgemag,compactness);
    else
        newPerformSuperpixelSLIC(kseedsl,kseedsa, kseedsb, kseedsx, kseedsy, klabels, STEP, edgemag,compactness);
    numlabels = kseedsl.size();

    int* nlabels = new int[sz];
    if(!newALGO)
        EnforceLabelConnectivity(klabels, m_width, m_height, nlabels, numlabels, double(sz)/double(STEP*STEP));
    else
        newEnforceLabelConnectivity(imgMat, klabels, m_width, m_height, nlabels, numlabels, double(sz)/double(STEP*STEP));
    {for(int i = 0; i < sz; i++ ) klabels[i] = nlabels[i];}
    if(nlabels) delete [] nlabels;
}
//===========================================================================
///	DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels
///
/// The input parameter ubuff conains RGB values in a 32-bit unsigned integers
/// as follows:
///
/// [1 1 1 1 1 1 1 1]  [1 1 1 1 1 1 1 1]  [1 1 1 1 1 1 1 1]  [1 1 1 1 1 1 1 1]
///
///        Nothing              R                 G                  B
///
/// The RGB values are accessed from (and packed into) the unsigned integers
/// using bitwise operators as can be seen in the function DoRGBtoLABConversion().
///
/// compactness value depends on the input pixels values. For instance, if
/// the input is greyscale with values ranging from 0-100, then a compactness
/// value of 20.0 would give good results. A greater value will make the
/// superpixels more compact while a smaller value would make them more uneven.
///
/// The labels can be saved if needed using SaveSuperpixelLabels()
//===========================================================================
void SLIC::DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(
        Mat                         imgMat,
        int*&						klabels,
        int&						numlabels,
        const int&					K,//required number of superpixels
        const double&               compactness,//weight given to spatial distance
        bool                        newALGO)    //new algorithm
{
    int width = imgMat.cols;
    int height = imgMat.rows;

    const int superpixelsize = 0.5+double(width*height)/double(K);
    DoSuperpixelSegmentation_ForGivenSuperpixelSize(imgMat,klabels,numlabels,superpixelsize,compactness, newALGO);
}

