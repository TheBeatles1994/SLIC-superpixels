// SLIC.h: interface for the SLIC class.
//===========================================================================
// This code implements the superpixel method described in:
//
// Radhakrishna Achanta, Appu Shaji, Kevin Smith, Aurelien Lucchi, Pascal Fua, and Sabine Susstrunk,
// "SLIC Superpixels",
// EPFL Technical Report no. 149300, June 2010.
//===========================================================================
//	Copyright (c) 2012 Radhakrishna Achanta [EPFL]. All rights reserved.
//===========================================================================
//////////////////////////////////////////////////////////////////////

#if !defined(_SLIC_H_INCLUDED_)
#define _SLIC_H_INCLUDED_


#include <vector>
#include <string>
#include <algorithm>
#include <opencv2/highgui/highgui.hpp>
#include <opencv/cv.h>

using namespace cv;
using namespace std;


class SLIC  
{
public:
    SLIC();
    virtual ~SLIC();
    //SLIC算法启动
    void runSLIC(Mat imgMat, int K, double compactness,  bool newALGO = false);
    //============================================================================
    // Superpixel segmentation for a given step size (superpixel size ~= step*step)
    //============================================================================
    void DoSuperpixelSegmentation_ForGivenSuperpixelSize(Mat                         imgMat,
            int*&						klabels,
            int&						numlabels,
            const int&					superpixelsize,
            const double&               compactness,
            bool                        newALGO);
    //============================================================================
    // Superpixel segmentation for a given number of superpixels
    //============================================================================
    void DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(
            Mat                         imgMat,
            int*&						klabels,
            int&						numlabels,
            const int&					K,//required number of superpixels
            const double&               compactness,//10-20 is a good value for CIELAB space
            bool                        newALGO);
    //============================================================================
    // Save superpixel labels in a text file in raster scan order
    // 未适应OpenCV
    //============================================================================
    void SaveSuperpixelLabels(
            const int*&					labels,
            const int&					width,
            const int&					height,
            const string&				filename,
            const string&				path);
    //============================================================================
    // Save supervoxel labels in a text file in raster scan, depth order
    // 未适应OpenCV
    //============================================================================
    void SaveSupervoxelLabels(
            const int**&				labels,
            const int&					width,
            const int&					height,
            const int&					depth,
            const string&				filename,
            const string&				path);
    //============================================================================
    // Function to draw boundaries around superpixels of a given 'color'.
    // Can also be used to draw boundaries around supervoxels, i.e layer by layer.
    //============================================================================
    void DrawContoursAroundSegments(
            Mat                     imgMat,
            int*&					labels,
            const int&				width,
            const int&				height,
            const unsigned int&				color );
    inline Mat getSlicImg(){return slicImg;}
    inline void setLambda(double lambda){this->lambda = lambda;}
private:
    //============================================================================
    // The main SLIC algorithm for generating superpixels
    //============================================================================
    void PerformSuperpixelSLIC(
            vector<double>&				kseedsl,
            vector<double>&				kseedsa,
            vector<double>&				kseedsb,
            vector<double>&				kseedsx,
            vector<double>&				kseedsy,
            int*&						klabels,
            const int&					STEP,
            const vector<double>&		edgemag,
            const double&				m = 10.0);
    void newPerformSuperpixelSLIC(
            Mat                         imgMat,
            vector<double>&				kseedsl,
            vector<double>&				kseedsa,
            vector<double>&				kseedsb,
            vector<double>&				kseedsx,
            vector<double>&				kseedsy,
            int*&						klabels,
            const int&					STEP,
            const vector<double>&		edgemag,
            const double&				m = 10.0);
    //============================================================================
    // Pick seeds for superpixels when step size of superpixels is given.
    //============================================================================
    void GetLABXYSeeds_ForGivenStepSize(
            vector<double>&				kseedsl,
            vector<double>&				kseedsa,
            vector<double>&				kseedsb,
            vector<double>&				kseedsx,
            vector<double>&				kseedsy,
            const int&					STEP,
            const bool&					perturbseeds,
            const vector<double>&		edgemag);
    //============================================================================
    // Pick seeds for supervoxels
    //============================================================================
    void GetKValues_LABXYZ(
            vector<double>&				kseedsl,
            vector<double>&				kseedsa,
            vector<double>&				kseedsb,
            vector<double>&				kseedsx,
            vector<double>&				kseedsy,
            vector<double>&				kseedsz,
            const int&					STEP);
    //============================================================================
    // Move the superpixel seeds to low gradient positions to avoid putting seeds
    // at region boundaries.
    //============================================================================
    void PerturbSeeds(
            vector<double>&				kseedsl,
            vector<double>&				kseedsa,
            vector<double>&				kseedsb,
            vector<double>&				kseedsx,
            vector<double>&				kseedsy,
            const vector<double>&		edges);
    //============================================================================
    // Detect color edges, to help PerturbSeeds()
    //============================================================================
    void DetectLabEdges(
            const double*				lvec,
            const double*				avec,
            const double*				bvec,
            const int&					width,
            const int&					height,
            vector<double>&				edges);
    //============================================================================
    // sRGB to XYZ conversion; helper for RGB2LAB()
    //============================================================================
    void RGB2XYZ(
            const int&					sR,
            const int&					sG,
            const int&					sB,
            double&						X,
            double&						Y,
            double&						Z);
    //============================================================================
    // sRGB to CIELAB conversion (uses RGB2XYZ function)
    //============================================================================
    void RGB2LAB(
            const int&					sR,
            const int&					sG,
            const int&					sB,
            double&						lval,
            double&						aval,
            double&						bval);
    //============================================================================
    // sRGB to CIELAB conversion for 2-D images
    //============================================================================
    void DoRGBtoLABConversion(
            Mat                         imgMat,
            double*&					lvec,
            double*&					avec,
            double*&					bvec);
    //============================================================================
    // sRGB to CIELAB conversion for 3-D volumes
    // 未适应OpenCV
    //============================================================================
    void DoRGBtoLABConversion(
            unsigned int**&				ubuff,
            double**&					lvec,
            double**&					avec,
            double**&					bvec);
    //============================================================================
    // Post-processing of SLIC segmentation, to avoid stray labels.
    //============================================================================
    void EnforceLabelConnectivity(
            const int*					labels,
            const int					width,
            const int					height,
            int*&						nlabels,//input labels that need to be corrected to remove stray labels
            int&						numlabels,//the number of labels changes in the end if segments are removed
            const int&					K); //the number of superpixels desired by the user
    void newEnforceLabelConnectivity(
            Mat                         imgMat,
            const int*					labels,
            const int					width,
            const int					height,
            int*&						nlabels,//input labels that need to be corrected to remove stray labels
            int&						numlabels,//the number of labels changes in the end if segments are removed
            const int&					K); //the number of superpixels desired by the user
private:
    int										m_width;
    int										m_height;
    int										m_depth;

    double*									m_lvec;
    double*									m_avec;
    double*									m_bvec;

    double**								m_lvecvec;
    double**								m_avecvec;
    double**								m_bvecvec;

    Mat slicImg;    //原始方法的结果
    double lambda = 0.8;  //lambda
};

#endif // !defined(_SLIC_H_INCLUDED_)
