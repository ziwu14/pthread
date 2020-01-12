#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

using namespace std;
double calc_time(struct timespec start, struct timespec end);
void initMatrixFromFile(ifstream& file, int N, int** matrix);
void freeMatrix(int N, int** matrix);
//----------------------------Point-------------------------------------

class Point {
   private:
    int x;
    int y;
    int height;
    double absorbRate;
    double totalAbsorption;
    double remain;
    double trickleAmount;
    vector<pair<int, int>> neighbors;

   public:
    //construct and destruct
    Point(int x_,
          int y_,
          int height_,
          double absorbRate_,
          vector<pair<int, int>> neighbors_) : x(x_),
                                               y(y_),
                                               height(height_),
                                               absorbRate(absorbRate_),
                                               totalAbsorption(0),
                                               remain(0),
                                               trickleAmount(0),
                                               neighbors(neighbors_) {}
    //getters
    vector<pair<int, int>>& getNeighobors();
    double getTrickleAmount();
    double getTotalAbsorption();
    double getRemain();
    //Processing point during simulation
    void receive();
    void absorb();
    void calculate();
    void update(double trickleAmount);
    bool isNotFinished() { return remain; }

    //used to test content in the Point object
    friend ostream& operator<<(ostream& stream, const Point& data);
};

ostream& operator<<(ostream& stream, const Point& data) {
    stream << data.neighbors.size();
    return stream;
}

void Point::receive() {
    remain += 1;
}

/**
 * depends on whether remain > absorption rate
*/
void Point::absorb() {
    if (remain == 0)
        return;

    if (remain >= absorbRate) {
        totalAbsorption += absorbRate;
        remain -= absorbRate;
    } else {
        totalAbsorption += remain;
        remain = 0;
    }
}
/**
 * Note: the total trickle amount at each timestamp is 1
 * remain >= 1, trickle total amount 1 to neighbors
 * remain < 1, trickle all remain to neighbors
*/
void Point::calculate() {
    if (neighbors.size() == 0 || remain == 0) {
        trickleAmount = 0;
        return;
    }

    if (remain >= 1) {
        trickleAmount = (double)1 / neighbors.size();  //Note: rember casting 1 to double otherwise integer division
        remain -= 1;
    } else {
        trickleAmount = remain / neighbors.size();
        remain = 0;
    }
}

void Point::update(double trickleAmount) {
    remain += trickleAmount;
}

vector<pair<int, int>>& Point::getNeighobors() {
    return neighbors;
}

double Point::getTrickleAmount() {
    return trickleAmount;
}

double Point::getTotalAbsorption() {
    return totalAbsorption;
}

double Point::getRemain() {
    return remain;
}

//-------------------------------simulator----------------------------------

/**
 * contains matrix of Point objects
 * with configs:
 * -- matrix dimension
 * -- number of threads for parallelism
 * -- rain duration/time steps
 * -- absorption rate
*/
class Simulator {
   public:
    //construct and destruct
    Simulator(int N_, int timeSteps_, double absorbRate_, int** elevationMatrix);
    ~Simulator();
    // start the simulation and print result
    void run();
    void printResult();

   private:
    Point** points;
    int N;
    int timeSteps;
    double absorbRate;
    int totalTime;
    double elapsed_ns;

    //return vectors of coordinates indicating the set of lowest neighbors
    vector<pair<int, int>> findLowestNeighbors(int row, int col, int** matrix);
    void processDuringRain();
    void processAfterRain();
  // used to test container initialization
    void printPoints1();
    void printPoints2();
    void printNeighbors();

};

Simulator::Simulator(
    int N_,
    int timeSteps_,
    double absorbRate_,
    int** elevationMatrix) : N(N_),
                             timeSteps(timeSteps_),
                             absorbRate(absorbRate_),
                             totalTime(timeSteps_) {
    points = (Point**)malloc(N * sizeof(Point*));
    for (int row = 0; row < N; row++) {
        points[row] = (Point*)malloc(N * sizeof(Point));
        for (int col = 0; col < N; col++)
            points[row][col] = Point(row,
                                     col,
                                     elevationMatrix[row][col],
                                     absorbRate,
                                     findLowestNeighbors(row,
                                                         col,
                                                         elevationMatrix));
    }
}

vector<pair<int, int>> Simulator::findLowestNeighbors(
    int row,
    int col,
    int** matrix) {
    vector<pair<int, int>> ans;
    //neighbor at most matrix[row][col]-1
    int min = matrix[row][col] - 1;
    if (row > 0 && matrix[row - 1][col] <= min) {
        min = matrix[row - 1][col];
        ans.push_back({row - 1, col});
    }
    if (row < N - 1 && matrix[row + 1][col] <= min) {
        if (matrix[row + 1][col] < min) {
            min = matrix[row + 1][col];
            ans.clear();
        }
        ans.push_back({row + 1, col});
    }
    if (col > 0 && matrix[row][col - 1] <= min) {
        if (matrix[row][col - 1] < min) {
            min = matrix[row][col - 1];
            ans.clear();
        }
        ans.push_back({row, col - 1});
    }
    if (col < N - 1 && matrix[row][col + 1] <= min) {
        if (matrix[row][col + 1] < min) {
            min = matrix[row][col + 1];
            ans.clear();
        }
        ans.push_back({row, col + 1});
    }
    return ans;
}

Simulator::~Simulator() {
    for (int i = 0; i < N; i++)
        free(points[i]);
    free(points);
}

/**
 * two traversals to ensure deterministic result
 * 1st traversal:
 * -- receive 1 raindrop to remains for each point (if still raining)
 * -- if there are raindrops remaining on a point, absorb water into the point in absorbRate
 * -- calculate the # of raindrops that will trickle to the neighbors or no neighbors
 * 2nd traversal:
 * -- for each point, use the calculated # raindrops that will trickle to the neighbors
 *      to update the remains at each lowest neighobor
*/
void Simulator::run() {
    //printNeighbors();
    struct timespec start_time, end_time;
    clock_gettime(CLOCK_MONOTONIC, &start_time);
    processDuringRain();
    processAfterRain();
    clock_gettime(CLOCK_MONOTONIC, &end_time);
    elapsed_ns = calc_time(start_time, end_time);
}

void Simulator::processDuringRain() {
    for (int t = 0; t < timeSteps; t++) {
        //1st traversal
        for (int row = 0; row < N; row++) {
            for (int col = 0; col < N; col++) {
                Point* curr = &points[row][col];
                curr->receive();
                curr->absorb();
                curr->calculate();
            }
        }
       
        //2nd traversal
        for (int row = 0; row < N; row++) {
            for (int col = 0; col < N; col++) {
                Point* curr = &points[row][col];
                for (auto& coord : curr->getNeighobors())
                    points[coord.first][coord.second].update(curr->getTrickleAmount());
            }
        }
    }
}
void Simulator::processAfterRain() {
    while (1) {
        bool isNotFinished = false;
        //1st traversal
        for (int row = 0; row < N; row++) {
            for (int col = 0; col < N; col++) {
                Point* curr = &points[row][col];
                curr->absorb();
                curr->calculate();
                if (curr->isNotFinished())
                    isNotFinished = true;
            }
        }

        if (!isNotFinished) {
            totalTime++;
            break;
        }

        //2nd traversal
        for (int row = 0; row < N; row++) {
            for (int col = 0; col < N; col++) {
                Point* curr = &points[row][col];
                if (curr->getTrickleAmount() == 0)
                    continue;
                for (auto& coord : curr->getNeighobors())
                    points[coord.first][coord.second].update(curr->getTrickleAmount());
            }
        }

        totalTime++;
    }
}

void Simulator::printPoints1() {
    double sum = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            sum += (points[i][j].getRemain() + points[i][j].getTotalAbsorption() + points[i][j].getNeighobors().size() * points[i][j].getTrickleAmount());
    }
    cout << "total amount of rain drops after 1st traversal are " << sum << endl;
}

void Simulator::printPoints2() {
    double sum = 0;
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            sum += (points[i][j].getRemain() + points[i][j].getTotalAbsorption());
    }
    cout << "total amount of rain drops after 2nd traversal are " << sum << endl;
}

void Simulator::printNeighbors() {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++)
            cout << points[i][j].getTrickleAmount() << " ";
        cout << endl;
    }
}

void Simulator::printResult() {
    //cout << fixed;
    cout << "Rainfall simulation took " << totalTime << " time steps to complete." << endl
         << "Runtime = " << elapsed_ns / 1000000000 << " seconds" << endl
         << endl
         << "The following grid shows the number of raindrops absorbed at each point:" << endl;
    cout << noshowpoint;
    for (int row = 0; row < N; row++) {
        for (int col = 0; col < N; col++)
            cout << setw(8) << points[row][col].getTotalAbsorption();
        cout << endl;
    }
}

//----------------------non-class--------------------------------------

int main(int argc, char const* argv[]) {
    //invalid usage
    if (argc != 5) {
        cerr << "Usage: ./rainfall <# of time steps> "
                "<absorption rate> <elevation matrix dimension> <elevation file>"
             << endl;
        exit(1);
    }

    //invalid file
    ifstream file(argv[4]);
    if (!file) {
        cerr << "invalid file" << endl;
        exit(1);
    }

    //init
    int timeSteps = atoi(argv[1]);
    double absorbRate = stod(argv[2]);
    int N = atoi(argv[3]);

    int** matrix = (int**)malloc(N * sizeof(int*));
    for (int i = 0; i < N; i++) {
        matrix[i] = (int*)malloc(N * sizeof(int));
    }
    initMatrixFromFile(file, N, matrix);
    //run simulator
    Simulator simulator(N, timeSteps, absorbRate, matrix);
    freeMatrix(N, matrix);
    simulator.run();
    simulator.printResult();
    //clear
    return 0;
}

void initMatrixFromFile(ifstream& file, int N, int** matrix) {
    int num;
    for (int row = 0; row < N; row++) {
        for (int col = 0; col < N; col++) {
            file >> num;
            matrix[row][col] = num;
        }
    }
}

void freeMatrix(int N, int** matrix) {
    for (int i = 0; i < N; i++)
        free(matrix[i]);
    free(matrix);
}

double calc_time(struct timespec start, struct timespec end) {
    //get timestamp in nanosecond
    double start_sec = (double)start.tv_sec * 1000000000.0 + (double)start.tv_nsec;
    double end_sec = (double)end.tv_sec * 1000000000.0 + (double)end.tv_nsec;

    if (end_sec < start_sec) {
        return 0;
    } else {
        return end_sec - start_sec;
    }
}
