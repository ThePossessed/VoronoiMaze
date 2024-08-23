#include <iostream>
#include <queue>
#include <set>
#include <math.h>
#include <time.h>
#include <vector>
#include <cmath>
#include <map>
#include <stack>
#include <float.h>
#define WIDTH 1000
#define HEIGHT 500
#define r 15
#define n 1000
#define dim 2
#define k 300

using namespace std;

// Notation for working with points
typedef std::pair<double, double> point;
#define x first
#define y second

// Arc, event, and segment datatypes
struct arc;
struct seg;

struct event {
   double x;
   point p;
   arc *a;
   bool valid;

   event(double xx, point pp, arc *aa)
      : x(xx), p(pp), a(aa), valid(true) {}
};

struct arc {
   point p;
   arc *prev, *next;
   event *e;

   seg *s0, *s1;

   arc(point pp, arc *a=0, arc *b=0)
    : p(pp), prev(a), next(b), e(0), s0(0), s1(0) {}
};

vector<seg*> output;  // Array of output segments.
map<point, std::vector<seg *>> allRooms; // Map of all points to their walls
map<seg *, std::vector<point>> allWalls; // Map of all walls to their points
map<seg *, int> includedWall;
map<point, int> visited;

struct seg {
   point start, end;
   bool done;

   seg(point p)
      : start(p), end(0,0), done(false)
   { output.push_back(this); }

   seg(point p1, point p2, bool done): start(p1), end(p2), done(done) {};
   
   // Set the end point and mark as "done."
   void finish(point p) { if (done) return; end = p; done = true; }
   bool operator==(const seg& a) {
      return (start.x == a.start.x && start.y == a.start.y)
         && (end.x == a.end.x && end.y == a.end.y);
   }
};

arc *root = 0; // First item in the parabolic front linked list

// Function declarations
bool check(int count, double x, double y, double grid[][dim]);
void process_point();
void process_event();
void front_insert(point  p);

bool circle(point a, point b, point c, double *x, point *o);
void check_circle_event(arc *i, double x0);

bool intersect(point p, arc *i, point *result = 0);
point intersection(point p0, point p1, double l);

void finish_edges();
void clean_edges();
void create_border();
void initiate_room_wall();
void maze_generation();
bool hasNeighbors(point p);
bool compare_point(point p1, point p2);
bool find_intersect(point *&p1, point *&p2, double x3, double y3, double x4, double y4);
bool intersect_segments(point *&p1, point *&p2, double x0, double y0, double x1, double y1);
void print_output();

// "Greater than" comparison, for reverse sorting in priority queue.
struct gt {
   bool operator()(point a, point b) {return a.x==b.x ? a.y>b.y : a.x>b.x;}
   bool operator()(event *a, event *b) {return a->x>b->x;}
};


// Bounding box coordinates.
double X0 = 0, X1 = 0, Y0 = 0, Y1 = 0;