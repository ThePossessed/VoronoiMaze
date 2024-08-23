#include "efficientPoisson.hh"

double dist(double x1, double y1, double x2, double y2) {
    return std::sqrt(std::pow(x2 - x1, 2) + std::pow(y2 - y1, 2));
}

bool isValidPoint(std::vector<std::vector<PVector*>>& grid, double cellsize,
                  int gwidth, int gheight, const PVector& p, double radius, int width, int height) {
    // Make sure the point is on the screen
    if (p.x < 0 || p.x >= width || p.y < 0 || p.y >= height)
        return false;

    // Check neighboring eight cells
    int xindex = std::floor(p.x / cellsize);
    int yindex = std::floor(p.y / cellsize);
    int i0 = std::max(xindex - 1, 0);
    int i1 = std::min(xindex + 1, gwidth - 1);
    int j0 = std::max(yindex - 1, 0);
    int j1 = std::min(yindex + 1, gheight - 1);

    for (int i = i0; i <= i1; i++) {
        for (int j = j0; j <= j1; j++) {
            if (grid[i][j] != nullptr) {
                if (dist(grid[i][j]->x, grid[i][j]->y, p.x, p.y) < radius)
                    return false;
            }
        }
    }

    // If we get here, return true
    return true;
}

void insertPoint(std::vector<std::vector<PVector*>>& grid, double cellsize, const PVector& point) {
    int xindex = std::floor(point.x / cellsize);
    int yindex = std::floor(point.y / cellsize);
    grid[xindex][yindex] = new PVector(point.x, point.y);
}

std::vector<PVector> poissonDiskSampling(double radius, int trial, int width, int height) {
    const int N = 2;
    std::vector<PVector> points;
    std::vector<PVector> active;

    PVector p0(rand() % width, rand() % height);
    double cellsize = std::floor(radius / std::sqrt(N));

    int ncells_width = std::ceil(width / cellsize) + 1;
    int ncells_height = std::ceil(height / cellsize) + 1;

    std::vector<std::vector<PVector*>> grid(ncells_width, std::vector<PVector*>(ncells_height, nullptr));

    insertPoint(grid, cellsize, p0);
    points.push_back(p0);
    active.push_back(p0);

    while (!active.empty()) {
        int random_index = rand() % active.size();
        PVector p = active[random_index];

        bool found = false;
        for (int tries = 0; tries < trial; tries++) {
            double theta = static_cast<double>(rand()) / RAND_MAX * 360;
            double new_radius = radius + static_cast<double>(rand()) / RAND_MAX * radius;
            double pnewx = p.x + new_radius * std::cos(theta * 3.1415 / 180.0);
            double pnewy = p.y + new_radius * std::sin(theta * 3.1415 / 180.0);
            PVector pnew(pnewx, pnewy);

            if (!isValidPoint(grid, cellsize, ncells_width, ncells_height, pnew, radius, width, height))
                continue;

            points.push_back(pnew);
            insertPoint(grid, cellsize, pnew);
            active.push_back(pnew);
            found = true;
            break;
        }

        if (!found)
            active.erase(active.begin() + random_index);
    }

    return points;
}

// int main() {
//     int width = 1000;
//     int height = 500;
//     double radius = 15;
//     int k = 3000;

//     std::vector<PVector> points = poissonDiskSampling(radius, k, width, height);

//     for (const PVector& point : points) {
//         std::cout << point.x << " " << point.y << std::endl;
//     }

//     return 0;
// }
