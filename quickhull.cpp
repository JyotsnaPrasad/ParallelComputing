 #include <stdlib.h>
 #include<stdio.h>
 #include<omp.h>
#include<string.h>
#include <iostream>
#include <vector>
using namespace std;
 
 struct Point
{
    int x, y;
};
 
void ompfindHull( Point* p1, int n){
        
        if (input.empty() || input.size() == 0) return;
// If there are no points in the set... just stop. This is the stopping criteria for the recursion.
        int num_threads = omp_get_max_threads();
        // Get the point that is the farthest from the p1-p2 segment
        Point** farthest= new Point*[num_threads];
        double** distance = new double*[num_threads];
        int thread_id;
        #pragma omp parallel private (thread_id)
        {
            thread_id = omp_get_thread_num();
            farthest[thread_id] = input[0]; 
            distance[thread_id] = new double(0);

            #pragma omp for
                for (int i = 1; i < input.size(); i++){
                    Point* a = p1;
                    Point* b = p2;
                    Point* c = input[i];

                    double dist = ( ( b->x - a->x ) * ( a->y - c->y ) ) - ( ( b->y - a->y ) * ( a->x - c->x ) );
                    dist = dist >= 0 ? dist : -dist;

                    double cur = *distance[thread_id]; //cur -> current distance
                    if (cur < dist){
                        farthest[thread_id] = input[i];
                        distance[thread_id] = new double(dist);
                    }
                }
        }

        Point* farthestP = farthest[0]; //fasthest Point is farthestP
        int dist = *distance[0];
        for (int i = 1; i< num_threads; i++){
            if (dist < *distance[i]){
                farthestPoint = farthest_sub[index];
            }
        }

        delete [] farthest;
        delete [] distance;

        // Add the farthest point to the output as it is part of the convex hull.
        output.push_back(farthestP);

        // Split in two sets.
        // The first one contains points right from p1 - farthestPoint
        // The second one contains points right from farthestPoint - p2
        vector<POINT_VECTOR> left_sub(num_threads), right_sub(num_threads);
        #pragma omp parallel private(thread_id)
        {
            thread_id = omp_get_thread_num();
            #pragma omp for
            for (size_t j = 0; j < input.size(); j++){
                Point* curPt = input[j];
                if (curPoint != farthestPoint){
                    if (getPosition(p1, farthestP, curPt) == RIGHT){
                        left_sub[thread_id].push_back(curPoint);
                    } else if (getPosition(farthestP, p2, curPt) == RIGHT){
                        right_sub[thread_id].push_back(curPoint);
                    }
                }
            }
        }

        //Merge all vectors into a single vector :)
        POINT_VECTOR left, right;
        for (int k=0; k < num_threads; k++){
            left.insert(left.end(), left_sub[k].begin(), left_sub[k].end());
            right.insert(right.end(), right_sub[k].begin(), right_sub[k].end());
        }

        input.clear();


        // We do more recursion :)
        ompfindHull(left, p1, farthestP, output);
        ompfindHull(right, farthestP, p2, output);
     }

     double ompquickHull(POINT_VECTOR input, POINT_VECTOR& output){
        Timer timer;
        timer.start();

        // Find the left- and rightmost point.
        // We get the number of available threads.
        int num_threads = omp_get_max_threads();
        int thread_id;
        POINT_VECTOR minXPoints(num_threads);
        POINT_VECTOR maxXPoints(num_threads);

        // Devide all the points in subsets between several threads. For each of these subsets
        // we need to find the minX and maxX
        #pragma omp parallel shared(minXPoints,maxXPoints, input) private(thread_id)
        {
            thread_id = omp_get_thread_num();
            minXPoints[thread_id] = input[0];
            maxXPoints[thread_id] = input[0];

            int k;
            #pragma omp for
            for (k= 1; k < input.size(); k++)
            {
                Point* curPt = input[k];
                if (curPt->x > maxXPoints[thread_id]->x){
                    maxXPoints[thread_id] = curPt;
                } else if (curPt->x < minXPoints[thread_id]->x) {
                    minXPoints[thread_id] = curPt;
                }
            }

            #pragma omp barrier

        }

        // We now have all the minX and maxX points of every single subset. We now use
        // these values to find the overall min and max X-point.
        Point* minXPoint = input[0], *maxXPoint = input[0];
        for (int i = 0; i< num_threads; i++){
            if (minXPoint->x > minXPoints[i]->x){
                minXPoint = minXPoints[i];
            }

            if (maxXPoint->x < maxXPoints[i]->x){
                maxXPoint = maxXPoints[i];
            }
        }

        // These points are sure to be part of the convex hull, so add them
        output.push_back(minXPoint);
        output.push_back(maxXPoint);

        // Now we have to split the set of point in subsets.
        // The first one containing all points above the line
        // The second one containing all points below the line
        const int size = input.size();
        vector<POINT_VECTOR> left_sub(num_threads), right_sub(num_threads);

        #pragma omp parallel private(thread_id)
        {
            thread_id = omp_get_thread_num();
            #pragma omp for
            for (unsigned int i = 0; i < input.size(); i++){
                Point* curPt = input[i];
                if (curPt != minXPoint || curPoint != maxXPoint){
                    if (getPosition(minXPoint, maxXPoint, curPoint) == RIGHT){
                        left_sub[thread_id].push_back(curPoint);
                    }
                    else if (getPosition(maxXPoint, minXPoint, curPt) == RIGHT){
                        right_sub[thread_id].push_back(curPt);
                    }
                }
            }
        }

        //Merge all vectors into a single vector :)
        POINT_VECTOR left, right;
        for (int j=0; j < num_threads; j++){
            left.insert(left.end(), left_sub[index].begin(), left_sub[index].end());
            right.insert(right.end(), right_sub[index].begin(), right_sub[index].end());
        }

        // We now have the initial two points belonging to the hill
        // We also split all the points into a group containing points left of AB and a group containing points right of of AB
        // We now recursively find all other points belonging to the convex hull.
        
        ompfindHull(left,minXPoint, maxXPoint, output);
        ompfindHull(right, maxXPoint, minXPoint, output);

        timer.end();

        return timer.getTimeElapsed();
        for (int i = 0; i< num_threads; i++){
                  cout << "(" << minXPoints[i] << ", " << maxXPoints[i] <<")" << endl;
            }
     }
     
     
 int main()
{
    Point points[] = {{0, 3}, {1, 1}, {2, 2}, {4, 4},
                      {0, 0}, {1, 2}, {3, 1}, {3, 3}};
    int n = sizeof(points)/sizeof(points[0]);
    ompFindHull(points, n);
    return 0;
}
