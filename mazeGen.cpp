#include "mazeGen.hh"
#include "efficientPoisson.hh"

priority_queue<point,  vector<point>,  gt> points; // site events
priority_queue<event*, vector<event*>, gt> events; // circle events
double active_grid[(int) (WIDTH*HEIGHT/(r*r))][dim];
int numSite = 0;
int startCenter;
int endCenter;

int main()
{
   // Read points from input.
   point p;
   std::vector<point> sites = { };

   srand(time(nullptr));

   /* Naive Approach for Poisson Disk Distribution */
   // srand(time(nullptr));

   // int trial = 0;
   // while (numSite < n) {      
   //    double x = ((double) rand() / (RAND_MAX)) * WIDTH;
   //    double y = ((double) rand() / (RAND_MAX)) * HEIGHT;

   //    if (check(numSite, x, y, active_grid)){
   //       active_grid[numSite][0] = x;
   //       active_grid[numSite][1] = y;
   //       sites.push_back({ x, y });
   //       numSite++;
   //    } else {
   //       trial++;
   //    }

   //    if (trial > k) {
   //       break;
   //    }
   // }
   /* Naive Approach for Poisson Disk Distribution */


   std::vector<PVector> efficientPoisson = poissonDiskSampling(r, 30, WIDTH*2, HEIGHT*2);
   for (const PVector& efficientPoint : efficientPoisson) {
      if (efficientPoint.x >= WIDTH && efficientPoint.x <= WIDTH*2 && efficientPoint.y >= HEIGHT && efficientPoint.y <= HEIGHT*2){
         active_grid[numSite][0] = efficientPoint.x-WIDTH;
         active_grid[numSite][1] = efficientPoint.y-HEIGHT;
         sites.push_back({ efficientPoint.x - WIDTH, efficientPoint.y - HEIGHT });
         numSite++;
      }
   }

   for (const point& p : sites) {
      points.push(p);

      // Keep track of bounding box size.
      if (p.x < X0) X0 = p.x;
      if (p.y < Y0) Y0 = p.y;
      if (p.x > X1) X1 = p.x;
      if (p.y > Y1) Y1 = p.y;
   }
   // Add 20% margins to the bounding box.
   double dx = (X1-X0+1)/5.0; double dy = (Y1-Y0+1)/5.0;
   X0 -= dx; X1 += dx; Y0 -= dy; Y1 += dy;
   X0 = 0; Y0 = 0; X1 = WIDTH; Y1 = HEIGHT;

   // cout << "No. of Points: " << points.size() << endl;

   // Process the queues; select the top element with smaller x coordinate.
   while (!points.empty())
      if (!events.empty() && events.top()->x <= points.top().x)
         process_event();
      else
         process_point();

   // After all points are processed, do the remaining circle events.
   while (!events.empty())
      process_event();

   finish_edges(); // Clean up dangling edges.
   clean_edges();
   create_border();
   initiate_room_wall();
   maze_generation();
   // cout << "Complete Maze" << endl;
   print_output(); // Output the voronoi diagram.
}

bool check(int count, double x1, double y1, double grid[][dim]) {
    bool ans = true;
    for (int i=0; i<count; i++)
    {
        int cur_x = grid[i][0];
        int cur_y = grid[i][1];

        double distance = sqrt(pow((x1 - cur_x), 2) + pow((y1 - cur_y), 2));

        if (distance < r) {
            ans = false;
            break;
        }
    }
    return ans;
}

void process_point()
{
   // Get the next point from the queue.
   point p = points.top();
   points.pop();

   // Add a new arc to the parabolic front.
   front_insert(p);
}

void process_event()
{
   // Get the next event from the queue.
   event *e = events.top();
   events.pop();

   if (e->valid) {
      // Start a new edge.
      seg *s = new seg(e->p);

      // Remove the associated arc from the front.
      arc *a = e->a;
      if (a->prev) {
         a->prev->next = a->next;
         a->prev->s1 = s;
      }
      if (a->next) {
         a->next->prev = a->prev;
         a->next->s0 = s;
      }

      // Finish the edges before and after a.
      if (a->s0) a->s0->finish(e->p);
      if (a->s1) a->s1->finish(e->p);

      // Recheck circle events on either side of p:
      if (a->prev) check_circle_event(a->prev, e->x);
      if (a->next) check_circle_event(a->next, e->x);
   }
   delete e;
}

void front_insert(point p)
{
   if (!root) {
      root = new arc(p);
      return;
   }

   // Find the current arc(s) at height p.y (if there are any).
   for (arc *i = root; i; i = i->next) {
      point z, zz;
      if (intersect(p,i,&z)) {
         // New parabola intersects arc i.  If necessary, duplicate i.
         if (i->next && !intersect(p,i->next, &zz)) {
            i->next->prev = new arc(i->p,i,i->next);
            i->next = i->next->prev;
         }
         else i->next = new arc(i->p,i);
         i->next->s1 = i->s1;

         // Add p between i and i->next.
         i->next->prev = new arc(p,i,i->next);
         i->next = i->next->prev;

         i = i->next; // Now i points to the new arc.

         // Add new half-edges connected to i's endpoints.
         i->prev->s1 = i->s0 = new seg(z);
         i->next->s0 = i->s1 = new seg(z);

         // Check for new circle events around the new arc:
         check_circle_event(i, p.x);
         check_circle_event(i->prev, p.x);
         check_circle_event(i->next, p.x);

         return;
      }
   }

   // Special case: If p never intersects an arc, append it to the list.
   arc *i;
   for (i = root; i->next; i=i->next) ; // Find the last node.

   i->next = new arc(p,i);  
   // Insert segment between p and i
   point start;
   start.x = X0;
   start.y = (i->next->p.y + i->p.y) / 2;
   i->s1 = i->next->s0 = new seg(start);
}

// Look for a new circle event for arc i.
void check_circle_event(arc *i, double x0)
{
   // Invalidate any old event.
   if (i->e && i->e->x != x0)
      i->e->valid = false;
   i->e = NULL;

   if (!i->prev || !i->next)
      return;

   double x;
   point o;

   if (circle(i->prev->p, i->p, i->next->p, &x,&o) && x > x0) {
      // Create new event.
      i->e = new event(x, o, i);
      events.push(i->e);
   }
}

// Find the rightmost point on the circle through a,b,c.
bool circle(point a, point b, point c, double *x, point *o)
{
   // Check that bc is a "right turn" from ab.
   if ((b.x-a.x)*(c.y-a.y) - (c.x-a.x)*(b.y-a.y) > 0)
      return false;

   // Algorithm from O'Rourke 2ed p. 189.
   double A = b.x - a.x,  B = b.y - a.y,
          C = c.x - a.x,  D = c.y - a.y,
          E = A*(a.x+b.x) + B*(a.y+b.y),
          F = C*(a.x+c.x) + D*(a.y+c.y),
          G = 2*(A*(c.y-b.y) - B*(c.x-b.x));

   if (G == 0) return false;  // Points are co-linear.

   // Point o is the center of the circle.
   o->x = (D*E-B*F)/G;
   o->y = (A*F-C*E)/G;

   // o.x plus radius equals max x coordinate.
   *x = o->x + sqrt( pow(a.x - o->x, 2) + pow(a.y - o->y, 2) );
   return true;
}

// Will a new parabola at point p intersect with arc i?
bool intersect(point p, arc *i, point *result)
{
   if (i->p.x == p.x) return false;

   double a,b;
   if (i->prev) // Get the intersection of i->prev, i.
      a = intersection(i->prev->p, i->p, p.x).y;
   if (i->next) // Get the intersection of i->next, i.
      b = intersection(i->p, i->next->p, p.x).y;

   if ((!i->prev || a <= p.y) && (!i->next || p.y <= b)) {
      result->y = p.y;

      result->x = (i->p.x*i->p.x + (i->p.y-result->y)*(i->p.y-result->y) - p.x*p.x)
                / (2*i->p.x - 2*p.x);

      return true;
   }
   return false;
}

// Where do two parabolas intersect?
point intersection(point p0, point p1, double l)
{
   point res, p = p0;

   double z0 = 2*(p0.x - l);
   double z1 = 2*(p1.x - l);

   if (p0.x == p1.x)
      res.y = (p0.y + p1.y) / 2;
   else if (p1.x == l)
      res.y = p1.y;
   else if (p0.x == l) {
      res.y = p0.y;
      p = p1;
   } else {
      // Use the quadratic formula.
      double a = 1/z0 - 1/z1;
      double b = -2*(p0.y/z0 - p1.y/z1);
      double c = (p0.y*p0.y + p0.x*p0.x - l*l)/z0
               - (p1.y*p1.y + p1.x*p1.x - l*l)/z1;

      res.y = ( -b - sqrt(b*b - 4*a*c) ) / (2*a);
   }
   // Plug back into one of the parabola equations.
   res.x = (p.x*p.x + (p.y-res.y)*(p.y-res.y) - l*l)/(2*p.x-2*l);
   return res;
}

void finish_edges()
{
   // Advance the sweep line so no parabolas can cross the bounding box.
   double l = X1 + (X1-X0) + (Y1-Y0);

   // Extend each remaining segment to the new parabola intersections.
   for (arc *i = root; i->next; i = i->next)
      if (i->s1)
         i->s1->finish(intersection(i->p, i->next->p, l*2));
}

bool compare_point(point p1, point p2){
   return (std::abs(p1.x-p2.x)<1e-6 && std::abs(p1.y-p2.y)<1e-6);
}

bool find_intersect(point *&p1, point *&p2, double x3, double y3, double x4, double y4){
   double denom = (((*p1).x - (*p2).x) * (y3 - y4)) - (((*p1).y - (*p2).y) * (x3 - x4));
   // cout << "Denom: " << denom << endl;
   if (std::abs(denom) < 1e-6) {
      return false; // Parallel
   } else {
      double t = ((((*p1).x - x3) * (y3 - y4)) - (((*p1).y - y3) * (x3 - x4)))/denom;
      double u = -((((*p1).x - (*p2).x) * ((*p1).y - y3)) - (((*p1).y - (*p2).y) * ((*p1).x - x3)))/denom;
      // cout << "t: " << t << " u: " << u << endl;
      if (t >= 0 && t <= 1 && u >= 0 && u <= 1){
         // Has intersection
         double intersectX = x3 + u*(x4-x3);
         double intersectY = y3 + u*(y4-y3);
         if ((*p1).x <= X0 || (*p1).x >= X1 || (*p1).y <= Y0 || (*p1).y >= Y1){
            // p1 is outside
            point* newP = new point();
            (*newP).x = intersectX;
            (*newP).y = intersectY;
            p1 = newP;
            if ((*p2).x <= X0 || (*p2).x >= X1 || (*p2).y <= Y0 || (*p2).y >= Y1) {
               // p2 is outside as well
               return false;
            } else {
               return true;
            }
         } else {
            // p2 is outside
            point* newP = new point();
            (*newP).x = intersectX;
            (*newP).y = intersectY;
            p2 = newP;
            if ((*p1).x <= X0 || (*p1).x >= X1 || (*p1).y <= Y0 || (*p1).y >= Y1) {
               // p2 is outside as well
               return false;
            } else {
               return true;
            }
         }
      } else {
         return false;
      }
   }
}

bool intersect_segments(point *&p1, point *&p2, double x0, double y0, double x1, double y1){
   // cout << "Original point: " << (*p1).x << " " << (*p1).y << " " << (*p2).x << " " << (*p2).y << endl;
   bool r1 = find_intersect(p1, p2, x0, y0, x0, y1);
   // cout << "New point: " << (*p1).x << " " << (*p1).y << " " << (*p2).x << " " << (*p2).y << endl;
   bool r2 = find_intersect(p1, p2, x0, y1, x1, y1);
   // cout << "New point: " << (*p1).x << " " << (*p1).y << " " << (*p2).x << " " << (*p2).y << endl;
   bool r3 = find_intersect(p1, p2, x1, y1, x1, y0);
   // cout << "New point: " << (*p1).x << " " << (*p1).y << " " << (*p2).x << " " << (*p2).y << endl;
   bool r4 = find_intersect(p1, p2, x1, y0, x0, y0);
   // cout << "New point: " << (*p1).x << " " << (*p1).y << " " << (*p2).x << " " << (*p2).y << endl;
   return (((*p1).x <= X0 || (*p1).x >= X1 || (*p1).y <= Y0 || (*p1).y >= Y1) && !r1 && !r2 && !r3 && !r4);
}

void clean_edges()
{
   vector<seg*>::iterator i;
   vector<seg*> new_output;

   // Merge parallel edges
   for (i = output.begin(); i != output.end(); i++) {
      point p0 = (*i)->start;
      point p1 = (*i)->end;
      double baseSlope = (p0.y - p1.y) / (p0.x - p1.x);
      bool merged = false;
      for (vector<seg*>::iterator j = new_output.begin(); j != new_output.end(); j++) {
         point currentp0 = (*j)->start;
         point currentp1 = (*j)->end;
         double currentSlope = (currentp0.y - currentp1.y) / (currentp0.x - currentp1.x);
         // start i = start j => merge
         if (compare_point(p0, currentp0) && std::abs(baseSlope-currentSlope)<1e-6){
            (*j)->start = p1;
            (*j)->end = currentp1;
            merged = true;
         } else if (compare_point(p1, currentp0) && std::abs(baseSlope-currentSlope)<1e-6){
            (*j)->start = p0;
            (*j)->end = currentp1;
            merged = true;
         } else if (compare_point(p1, currentp1) && std::abs(baseSlope-currentSlope)<1e-6){
            (*j)->start = p0;
            (*j)->end = currentp0;
            merged = true;
         } else if (compare_point(p0, currentp1) && std::abs(baseSlope-currentSlope)<1e-6){
            (*j)->start = p1;
            (*j)->end = currentp0;
            merged = true;
         } else {
            //do nothing
         }
      }
      if (!merged) {
         new_output.push_back(*i);
      }
      // cout << p0.x << " " << p0.y << " " << p1.x << " " << p1.y << endl;
   }

   // Trim edges
   for (i = new_output.begin(); i != new_output.end(); i++) {
      point* point1 = &((*i)->start);
      point* point2 = &((*i)->end);
      bool shouldDelete = intersect_segments(point1, point2, X0, Y0, X1, Y1);
      if (shouldDelete){
         new_output.erase(i);
         i--;
      } else {
         (*i)->start = *point1;
         (*i)->end = *point2;
      }
   }
   output = new_output;
}

void create_border() {
   std::vector<double> bordersX = { X0, X1 };
   std::vector<double> bordersY = { Y0, Y1 };
   std::vector<double> borderPoints = { };
   std::vector<double>::iterator i;
   for (i = bordersX.begin(); i != bordersX.end(); i++){
      priority_queue<double> pq;
      pq.push(-Y0);
      vector<seg*>::iterator j;
      for (j = output.begin(); j != output.end(); j++) {
         point p0 = (*j)->start;
         point p1 = (*j)->end;
         if (p0.x == *i) {
            pq.push(-p0.y);
         }
         if (p1.x == *i) {
            pq.push(-p1.y);
         }
      }
      pq.push(-Y1);
      while (!pq.empty()) {
         double startY = -pq.top();
         pq.pop();
         if (!pq.empty()){
            double endY = -pq.top();
            point start = { *i, startY };
            point end = { *i, endY };
            seg *borderSeg = new seg(start, end, true);
            output.push_back(borderSeg);
         } else {
            break;
         }
      }
   }

   for (i = bordersY.begin(); i != bordersY.end(); i++){
      priority_queue<double> pq;
      pq.push(-X0);
      vector<seg*>::iterator j;
      for (j = output.begin(); j != output.end(); j++) {
         point p0 = (*j)->start;
         point p1 = (*j)->end;
         if (p0.y == *i) {
            pq.push(-p0.x);
         }
         if (p1.y == *i) {
            pq.push(-p1.x);
         }
      }
      pq.push(-X1);
      while (!pq.empty()) {
         double startX = -pq.top();
         pq.pop();
         if (!pq.empty()){
            double endX = -pq.top();
            point start = { startX, *i };
            point end = { endX, *i };
            seg *borderSeg = new seg(start, end, true);
            output.push_back(borderSeg);
         } else {
            break;
         }
      }
   }
}

void initiate_room_wall() {
   for (int j = 0; j < numSite; j++){
      double curX = active_grid[j][0];
      double curY = active_grid[j][1];
      visited[{curX, curY}] = 0;
   }


   vector<seg*>::iterator i;
   for (i = output.begin(); i != output.end(); i++) {
      point p0 = (*i)->start;
      point p1 = (*i)->end;
      double midX = (p0.x + p1.x)/2;
      double midY = (p0.y + p1.y)/2;
      double dist1 = DBL_MAX;
      point point1;
      double dist2 = DBL_MAX;
      point point2;
      for (int j = 0; j < numSite; j++){
         double curX = active_grid[j][0];
         double curY = active_grid[j][1];
         double curDist = sqrt(pow(midX-curX, 2) + pow(midY-curY, 2));
         if (curDist < dist1) {
            dist2 = dist1;
            point2 = {point1.x, point1.y};
            dist1 = curDist;
            point1 = { curX, curY };
         } else if (curDist < dist2) {
            dist2 = curDist;
            point2 = { curX, curY };
         } else {
            // Nothing
         }
      }

      if (allRooms.find(point1) != allRooms.end()){
         allRooms[point1].push_back(*i);
      } else {
         allRooms[point1] = { *i };
      }
      if (!(((*i)->start.x == (*i)->end.x && (*i)->end.x == X0)
            || ((*i)->start.x == (*i)->end.x && (*i)->end.x == X1)
            || ((*i)->start.y == (*i)->end.y && (*i)->end.y == Y0)
            || ((*i)->start.y == (*i)->end.y && (*i)->end.y == Y1)
         ))
      {
         if (allRooms.find(point2) != allRooms.end()){
            allRooms[point2].push_back(*i);
         } else {
            allRooms[point2] = { *i };
         }
      }

      if (allWalls.find(*i) != allWalls.end()){
         allWalls[*i].push_back(point1);
      } else {
         allWalls[*i] = { point1 };
      }
      // Only add point 2 if Wall is not on the border
      if (!(((*i)->start.x == (*i)->end.x && (*i)->end.x == X0)
            || ((*i)->start.x == (*i)->end.x && (*i)->end.x == X1)
            || ((*i)->start.y == (*i)->end.y && (*i)->end.y == Y0)
            || ((*i)->start.y == (*i)->end.y && (*i)->end.y == Y1)
         ))
      {
         if (allWalls.find(*i) != allWalls.end()){
            allWalls[*i].push_back(point2);
         } else {
            allWalls[*i] = { point2 };
         }
      }
      includedWall[*i] = 1;
   }
}

bool hasNeighbors(point p) {
   std::vector<seg *> walls = allRooms[p];
   bool hasAnyAvailableNeighbor = false;
   vector<seg *>::iterator wall;
   for (wall = walls.begin(); wall != walls.end(); wall++){
      vector<point> currentWall = allWalls[*wall];
      vector<point>::iterator neighbor;
      for (neighbor = currentWall.begin(); neighbor != currentWall.end(); neighbor++){
         if ((*neighbor).x != p.x && (*neighbor).y != p.y){
            if (visited[{(*neighbor).x, (*neighbor).y}] == 0) {
               hasAnyAvailableNeighbor = true;
            }
         }
      }
   }
   return hasAnyAvailableNeighbor;
}

void maze_generation(){
   // cout << allRooms.size() << endl;
   // cout << allWalls.size() << endl;
   // cout << includedWall.size() << endl;
   // cout << visited.size() << endl;
   stack<point> sPoint;
   startCenter = ((double) rand() / (RAND_MAX)) * numSite;
   do {endCenter = ((double) rand() / (RAND_MAX)) * numSite;} while (endCenter == startCenter);
   point currentCenter = { active_grid[startCenter][0], active_grid[startCenter][1] };
   do {
      if (visited[{currentCenter.x, currentCenter.y}] == 0) {
         sPoint.push(currentCenter);
         visited[currentCenter] = 1;
      }
      double prevX = currentCenter.x;
      double prevY = currentCenter.y;
      int nextWall = ((double) rand() / (RAND_MAX)) * allRooms[currentCenter].size();
      seg * chosenWall = allRooms[currentCenter][nextWall];
      vector<point> wallPoints = allWalls[chosenWall];
      vector<point>::iterator neighbor;

      for (neighbor = wallPoints.begin(); neighbor != wallPoints.end(); neighbor++){
         // cout << (*neighbor).x << " " << (*neighbor).y << " " << currentCenter.x << " " << currentCenter.y << endl;
         if ((*neighbor).x != currentCenter.x && (*neighbor).y != currentCenter.y){
            currentCenter = { (*neighbor).x, (*neighbor).y };
            break;
         }
      }

      // if ((((chosenWall)->start.x == (chosenWall)->end.x && (chosenWall)->end.x == X0)
      //       || ((chosenWall)->start.x == (chosenWall)->end.x && (chosenWall)->end.x == X1)
      //       || ((chosenWall)->start.y == (chosenWall)->end.y && (chosenWall)->end.y == Y0)
      //       || ((chosenWall)->start.y == (chosenWall)->end.y && (chosenWall)->end.y == Y1)
      //    ))
      // {
      //    cout << "Start Debug" << endl;
      //    cout << "prev: " << prevX << " " << prevY << endl;
      //    cout << "cur: " << currentCenter.x << " " << currentCenter.y << endl;
      //    cout << "End Debug" << endl;
      // }

      if ((prevX != currentCenter.x || prevY != currentCenter.y)) {
         if (visited[{currentCenter.x, currentCenter.y}] == 0) {
            includedWall[chosenWall] = 0;

            // cout << "Visited " << visited[{currentCenter.x, currentCenter.y}] << endl;
            while (!sPoint.empty() && !hasNeighbors(currentCenter)){
               // cout << sPoint.top().x << " Inner " << sPoint.top().y << endl;
               currentCenter = sPoint.top();
               sPoint.pop();
            }
            // cout << "Queue size: " << sPoint.size() << endl;
         } else {
            currentCenter.x = prevX;
            currentCenter.y = prevY;
            continue;
         }
      }
      // cout << sPoint.top().x << " Outer " << sPoint.top().y << endl;
      // cout << "Status: " << !sPoint.empty() << endl;
   } while (!sPoint.empty());
   return;
}

void print_output()
{
   // Bounding box coordinates.
   // cout << "box: " << endl;
   // cout << X0 << " "<< X1 << " " << Y0 << " " << Y1 << endl;

   // Each output segment in four-column format.
   // cout << "segment: " << endl;
   cout << active_grid[startCenter][0] << " " << active_grid[startCenter][1] << endl;
   cout << active_grid[endCenter][0] << " " << active_grid[endCenter][1] << endl;
   vector<seg*>::iterator i;
   for (i = output.begin(); i != output.end(); i++) {
      if (includedWall[*i] == 1) {
         point p0 = (*i)->start;
         point p1 = (*i)->end;
         cout << p0.x << " " << p0.y << " " << p1.x << " " << p1.y << endl;
      }
   }
}