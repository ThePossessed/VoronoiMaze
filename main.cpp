#include <iostream>
#include <math.h>
#include <time.h>
#define WIDTH 1000
#define HEIGHT 500
#define r 20
#define n 100
#define dim 2
#define k 30


using namespace std;

bool check(int count, double x, double y, double grid[][dim]) {
    bool ans = true;
    for (int i=0; i<count; i++)
    {
        int cur_x = grid[i][0];
        int cur_y = grid[i][1];

        double distance = sqrt(pow((x - cur_x), 2) + pow((y - cur_y), 2));

        if (distance < r) {
            ans = false;
            break;
        }
    }
    return ans;
}

int main() {
    srand(time(nullptr));
    double active_grid[n][dim];

    int count = 0;
    while (count < n) {
        double x = rand() * WIDTH;
        double y = rand() * HEIGHT;

        if (check(count, x, y, active_grid)){
            active_grid[count][0] = x;
            active_grid[count][1] = y;
            count++;
        }
    }
}