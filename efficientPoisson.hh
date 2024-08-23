#include <iostream>
#include <vector>
#include <cmath>
#include <time.h>
#include <math.h>

struct PVector {
    double x, y;
    
    PVector() : x(0), y(0) {}
    PVector(double x, double y) : x(x), y(y) {}
};

double dist(double x1, double y1, double x2, double y2);

bool isValidPoint(std::vector<std::vector<PVector*>>& grid, double cellsize,
                  int gwidth, int gheight, const PVector& p, double radius, int width, int height);

void insertPoint(std::vector<std::vector<PVector*>>& grid, double cellsize, const PVector& point);

std::vector<PVector> poissonDiskSampling(double radius, int trial, int width, int height);

