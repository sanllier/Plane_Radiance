#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <process.h>
#include <Windows.h>
#include <vector>
#include <algorithm>
#include "rgbe.h"

#define PI 3.141592653589
#define DEFRESOLUTION 500
#define DEFTHREADS 2
#define LENLEAP 5

using namespace std;

typedef TripleSpace<double> Point;

/*****************************SPECTRUM*******************************/

struct Spectrum
{
    double *spec_;
    int specStart;
    int specEnd;
    int size;
    int shift;

    Spectrum():spec_(NULL), specStart(0), specEnd(0), size(0) {}
    Spectrum(string sFile, int from = 0, int to = 10000):spec_(NULL)
    { 
        readSpectrum(sFile, from, to); 
    }
    Spectrum(int start, int end):specStart(start), specEnd(end)
    {
        size = (end - start) / LENLEAP + 1;
        spec_ = (double *)calloc(size, sizeof(*spec_));
    }
    ~Spectrum() 
    { 
        free(spec_); 
    }

    void readSpectrum(string sourseFile, int from = 0, int to = 100000);
    void normilize(double total);
    inline double getValForLen(int len);
    void reverseRel(double rel);
    inline int getPosForLen(int len);
    void compShift(int len);
};

void
Spectrum::readSpectrum(string sourseFile, int from, int to)
{
    ifstream fInp(sourseFile);
    double trash = -1;
    free(spec_);
    spec_ = NULL;
    size = 0;
    specStart = -1;

    while (specStart < from) {
        fInp >> specStart >> trash;
    }   
    spec_ = (double *)(calloc(++size, sizeof(*spec_)));
    spec_[size - 1] = trash;
    while (fInp >> specEnd) {
        if (specEnd > to) break;
        spec_ = (double *)realloc(spec_, ++size * sizeof(*spec_));
        fInp >> spec_[size - 1];
    }
    specEnd = specEnd < to ? specEnd : to;

    fInp.close();
}

void
Spectrum::normilize(double total)
{
    double sum = 0.0, normConst = 0.0;

    for (int i = 0; i < size; ++i) {
        sum += spec_[i];
    }

    normConst = total / sum;

    for (int i = 0; i < size; ++i) {
        spec_[i] *= normConst;
    }
}

void
Spectrum::reverseRel(double rel)
{
    for (int i = 0; i < size; ++i) {
        spec_[i] = rel - spec_[i];
    }
}

void
Spectrum::compShift(int len) 
{
    shift = (len - specStart) / LENLEAP;
}

/****************************PARAMETERS*****************************/

struct Params
{
    string LspdFilePath;
    string AspdFilePath;
    string outFilePath;
    int N;
    int H;
    int P;
    int S;
    int cores;
    int theta;
    int nl;
    double cellSize;
    double halfCell;

    Params() {}
    Params(int argc, char **argv): N(DEFRESOLUTION), cores(DEFTHREADS), theta(0), outFilePath("img.hdr") 
    { 
        getParams(argc, argv); 
    }

    void getParams(int argc, char **argv);
};

void
Params::getParams(int argc, char **argv)
{
    if (!(argc % 2)) {
        cout << "Incorrect parameters.\n\n";
        exit(0);
    }

    bool isLspd = false, isH = false, isP = false, isAspd = false, isO = false, isAng = false;
    bool isNl = false, isRes = false, isCore = false;

    for (int pos = 1; pos < argc; pos += 2) {
        if (!strcmp(argv[pos], "-h")) {
            H = atoi(argv[pos + 1]);
            isH = true;
        }
        if (!strcmp(argv[pos], "-nl")) {
            nl = atoi(argv[pos + 1]);
            isNl = true;
        }
        if (!strcmp(argv[pos], "-res")) {
            N = atoi(argv[pos + 1]);
            isRes = true;
        }
        if (!strcmp(argv[pos], "-theta")) {
            theta = atoi(argv[pos + 1]);
            isAng = true;
        }
        if (!strcmp(argv[pos], "-core")) {
            cores = atoi(argv[pos + 1]) >= 1 ? atoi(argv[pos + 1]) : DEFTHREADS; 
            isCore = true;
        }
        if (!strcmp(argv[pos], "-S")) {
            S = atoi(argv[pos + 1]);
        }
        if (!strcmp(argv[pos], "-LSPD")) {
            LspdFilePath = argv[pos + 1];
            isLspd = true;
        }
        if (!strcmp(argv[pos], "-P")) {
            P = atoi(argv[pos + 1]);
            isP = true;
        }
        if (!strcmp(argv[pos], "-ASPD")) {
            AspdFilePath = argv[pos + 1];
            isAspd = true;
        }
        if (!strcmp(argv[pos], "-o")) {
            outFilePath = argv[pos + 1];
            isO = true;
        }
    }

    if (!isH) {
        cout << "Required parameter (-h) is not found.\n\n";
        exit(0);
    }
    if (!isP) {
        cout << "Required parameter (-P) is not found.\n\n";
        exit(0);
    }
    if (!isLspd) {
        cout << "Required parameter (-LSPD) is not found.\n\n";
        exit(0);
    }
    if (!isAspd) {
        cout << "Parameter (-ASPD) is not found.\nIdeal surface will be used.\n\n";
    }
    if (!isAng) {
        cout << "Parameter (-theta) is not found.\nUsing default value theta = 0\n\n";
    }
    if (!isNl) {
        cout << "Parameter (-nl) is not found.\nUsing 10õ10 grid of lights.\n\n";
    }
    if (!isCore) {
        cout << "Parameter (-core) is not found.\nUsing 2 threads.\n\n";
    }
    if (!isRes) {
        cout << "Parameter (-res) is not found.\nUsing 500x500 output resolution.\n\n";
    }
    if (!isO) {
        cout << "Parameter (-o) is not found.\nUsing default file name: ";
        cout << outFilePath << endl << endl;
    }

    cellSize = 2.0 * H / N; 
    halfCell = cellSize / 2.0;
}

/*******************************************************************/

struct RenderData:Params
{
    Spectrum *light;
    Spectrum *aspd;
    int from;
    int to;
    int numLightCells;
    double lightCellSize;

    RenderData():Params() {}
};

/*****************************GLOBAL********************************/

XYZ **plane;
Spectrum fx("xLambda.txt");
Spectrum fy("yLambda.txt");
Spectrum fz("zLambda.txt");

int specSpace;

/**************************************************************************/

void 
compSpecIntensity(Spectrum &light, int numCells)
{
    for (int i = 0; i < light.size; ++i) {
        light.spec_[i] /= 4 * PI * numCells * numCells;
    }
}

void
render(void *argv)
{
    RenderData *data = static_cast<RenderData *>(argv);

    double startShift = data->H - data->lightCellSize * (data->numLightCells / 2);
    startShift += data->numLightCells / 2 ? data->lightCellSize / 2.0 : 0;

    Spectrum irraddiance(data->light->specStart, data->light->specEnd);
    irraddiance.shift = data->light->shift;

    int lightShift = data->light->shift;
    int aspdShift = data->aspd != NULL ? data->aspd->shift : 0;
    int numCells = data->numLightCells;

    for (int k = 0; k < numCells; ++k) {
        for (int t = 0; t < numCells; ++t) {
            for (int i = data->from; i < data->to; ++i) {
                for (int q = 0; q < data->N; ++q) {
                    Point planePnt(i * data->cellSize + data->halfCell, 
                                   q * data->cellSize + data->halfCell);

                    double base = (startShift + t * data->lightCellSize) - data->H;
                    double tempY = base + data->H, tempZ;                    
                    if (base > 0) {
                        tempY = base * cos(data->theta * PI / 180.0);
                        tempZ = base * sin(data->theta * PI / 180.0);
                    } else {
                        tempY = base * (cos(-data->theta * PI / 180.0));
                        tempZ = base * (sin(-data->theta * PI / 180.0));
                    }

                    Point lightPnt(startShift + k * data->lightCellSize, 
                                   tempY + data->H,
                                   tempZ + data->H);

                    /**********************IRRADDIANCE************************/
                    double m1 = (lightPnt.x_ - planePnt.x_) * (lightPnt.x_ - planePnt.x_);
                    double m2 = (lightPnt.y_ - planePnt.y_) * (lightPnt.y_ - planePnt.y_);
                    double temp = sqrt(m1 + m2);

                    double pathLen =  sqrt(temp * temp + lightPnt.z_ * lightPnt.z_);
                    double cAng = (tempZ + data->H) / (pathLen * pathLen * pathLen);

                    for (int j = 0; j < specSpace ; ++j) { 
                        irraddiance.spec_[lightShift + j] = data->light->spec_[lightShift + j] * cAng;  
                    }
                    
                    /*********************************************************/

                    for (int j = 0; j < data->light->size ; ++j) {
                        double temp = irraddiance.spec_[lightShift + j] / PI;
                        if (data->aspd) temp *= data->aspd->spec_[aspdShift + j];  

                        plane[i][q].x_ += temp * fx.spec_[fx.shift + j] * LENLEAP;
                        plane[i][q].y_ += temp * fy.spec_[fy.shift + j] * LENLEAP;
                        plane[i][q].z_ += temp * fz.spec_[fz.shift + j] * LENLEAP;
                    }
                }
            }
        }
    }
}

/**************************************************************************/

int 
main(int argc, char **argv)
{
    Params args(argc, argv);
    vector<int> fMin, fMax;

    fMin.push_back(fx.specStart), fMax.push_back(fx.specEnd);
    fMin.push_back(fy.specStart), fMax.push_back(fy.specEnd);;
    fMin.push_back(fz.specStart), fMax.push_back(fz.specEnd);;

    Spectrum light(args.LspdFilePath, fx.specStart, fx.specEnd);
    light.normilize(args.P); 
    compSpecIntensity(light, args.nl);
    fMin.push_back(light.specStart), fMax.push_back(light.specEnd);

    Spectrum aspd;
    if (args.AspdFilePath.size()) {
        aspd.readSpectrum(args.AspdFilePath, fx.specStart, fx.specEnd);
        aspd.reverseRel(1.0);
        fMin.push_back(aspd.specStart), fMax.push_back(aspd.specEnd);;
    } 

    int min = *min_element(fMin.begin(), fMin.end());
    int max = *min_element(fMax.begin(), fMax.end());
    
    light.compShift(min);
    fx.compShift(min), fy.compShift(min), fz.compShift(min);
    if (args.AspdFilePath.size()) aspd.compShift(min);

    specSpace = (max - min) / LENLEAP + 1;

    plane = new XYZ *[args.N];
    for (int i = 0;  i < args.N; ++i) {
        plane[i] = new XYZ[args.N];
    } 

/***************************RENDERING******************************/

    int workStep = args.N / args.cores, start = 0;
    HANDLE *threadsHandle = new HANDLE[args.cores];
    RenderData *rData = new RenderData[args.cores];

    for (int i = 0; i < args.cores; ++i) {
        cout << "Starting thread #" << i + 1 << endl;
        *static_cast<Params *>(&rData[i]) = args;  
        rData[i].light = &light;
        rData[i].aspd = args.AspdFilePath.size() ? &aspd : NULL;
        rData[i].numLightCells = args.nl;
        rData[i].lightCellSize = sqrt(args.S) / args.nl;
        rData[i].from = start;
        rData[i].to = i == args.cores - 1 ? args.N : start += workStep;

        threadsHandle[i] = (HANDLE)_beginthread(render, 0, &rData[i]);   
    }

    WaitForMultipleObjects(args.cores, threadsHandle, true, INFINITE);
    
    delete[] threadsHandle;
    delete[] rData;

/******************************************************************/

    for (int i = 0; i < args.N; ++i) {
        for (int q = 0; q < args.N; ++q) {
            double r = (plane[i][q].x_ * 3.2404542 + 
                        plane[i][q].y_ * -1.5371385 + 
                        plane[i][q].z_ * -0.4985314); 
            double g = (plane[i][q].x_ * -0.9692660 + 
                        plane[i][q].y_ * 1.8760108 + 
                        plane[i][q].z_ * 0.0415560); 
            double b = (plane[i][q].x_ * 0.0556434 + 
                        plane[i][q].y_ * -0.2040259 + 
                        plane[i][q].z_ * 1.0572252); 
            plane[i][q].x_ = r, plane[i][q].y_ = g, plane[i][q].z_ = b;
        }
    }

    FILE *f = fopen(args.outFilePath.c_str(), "wb");
    RGBE_WriteHeader(f, args.N, args.N, NULL);
    _RGBE_WritePixels(f, plane, args.N, args.N);
    fclose(f);

    for (int i = 0;  i < args.N; ++i) {
        delete[] plane[i];
    } 
    delete[] plane;

    cout << "\nDone!\n";
    return 0;
}