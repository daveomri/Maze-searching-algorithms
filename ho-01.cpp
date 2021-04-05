#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <list>
#include <set>
#include <algorithm>
#include <vector>
#include <ctime>
#include <cstdlib>
#include <functional>
#include <queue>
#include <climits>
#include <cmath>
//#include <stack> 
#include <thread>
//library for time, that is used for randomized shuffle
#include <chrono>
//library for replacing substrings in string
#include <regex>
//libraries for linux
#include <sys/ioctl.h>
#include <unistd.h>
using namespace std;

/**
 * Implementation of searching algorithms in maze
 * 
 * This class is bit overloaded, it has methods for individual searching algorithms and 
 * also methods for returning graphical representation of visited nodes and also path from 
 * start place S to end place E.
 */
class AreaSearch{
    private:
        struct MemoryBlock{
            /**
             * Structure for representing place in maze
             * 
             * x, y - represent positions in maze
             * status - is helpfull information that is used to mark place as visited, open, atd.., it depends on algorithm
             */
            MemoryBlock * prev = NULL;
            int x;
            int y;
            int status = 0;
            MemoryBlock(int x, int y): x(x), y(y) {}
            MemoryBlock(): x(0), y(0) {}
        };
        /// char used for stoping program during processing searching algorithm,
        /// 'f' for speed run, 'q' for jump to result
        char pauseMe;
        /// string where maze map is stored
        string map = "";
        /// two dimensional array for representing individual places in maze
        MemoryBlock ** realMap = NULL;
        /// dont know what this is
        int mapLength = 0;
        /// width of maze map
        int mapWidth = 0;
        /// heigth of maze map
        int mapHeigth = 0;
        /// expanded places in maze
        int expadnedNodes = 0;
        /// length of path from start to end
        int pathLength = 0;
        /// start position
        MemoryBlock *Start = NULL;
        /// end position
        MemoryBlock *End = NULL;

    public:
        ~AreaSearch();

        AreaSearch();

        void printMap();

        void printStats();

        bool useAlgorithm(string);

        bool setValuesFromFile(string);

    protected:
        void cleanAll();

        int indexOf(int, int );

        void cleanRealMap();

        void findPath(MemoryBlock *);

        static int randomInt (int);

        static double heuristicEuclidean(int, int, int, int);

        static int heuristicManhattan(int, int, int, int);

        static int heuristicDiagonal(int, int, int, int);

        int endDistance(int, int);

        void bfs();

        void dfs();

        void randomSearch();

        void dijkstraSearch();

        void greedySearch();

        void aStarSearch();
};

/**
 * Destructor calls cleanAll() method to free all allocated memory 
 */
AreaSearch::~AreaSearch(){
    cleanAll();
}

/**
 * Mehod for freeing memory and setting variables to default values
 */
void AreaSearch::cleanAll(){
    map = "";
    mapLength = 0;
    mapWidth = 0;
    mapHeigth = 0;
    expadnedNodes = 0;
    pathLength = 0;

    if (Start)
        delete Start;
    Start = NULL;
    if (End)
        delete End;
    End = NULL;
    for (int i = 0; i < mapHeigth; i++){
        if (realMap[i]){
            delete [] realMap[i];
            realMap[i] = NULL;
        }
    }
    if (realMap)
        delete [] realMap;
    realMap = NULL;
}

/**
 * Constructor, that does nothing
 */
AreaSearch::AreaSearch(){
}

/**
 * Method that reads data from file and using it to set class variables, map, length, width etc..
 * 
 * @param fileName name of input text file with maze and locations of start and end
 * @return true if no complications occurs, false otherwise
 */
bool AreaSearch::setValuesFromFile(string fileName){
    cleanAll();
    ///string for saving temp. input text, string for storing map - - - - - -
    string input = "";
    /// setting file as ifstream
    ifstream file (fileName);

    if (file.is_open() && file.good())
    {
        while ( getline (file,input) && input.at(0)=='X')
        {
            mapWidth = input.length();
            mapHeigth++;
            mapLength+=mapWidth;
            map+=input;
        }
        ///loading of start - - - - - - - - - - - - - - - - - - - - - - - - -
        stringstream sstart(input);
        int start_x, start_y;
        sstart >> input;
        sstart >> start_x;
        sstart >> input;
        sstart >> start_y;
        Start = new MemoryBlock(start_x, start_y);
        ///loading end- - - - - - - - - - - - - - - - - - - - - - - - - - - -
        getline(file, input);
        stringstream ssend(input);
        int end_x, end_y;            
        ssend >> input;
        ssend >> end_x;
        ssend >> input;
        ssend >> end_y;
        End = new MemoryBlock(end_x, end_y);
        ///creating real map- - - - - - - - - - - - - - - - - - - - - - - - -
        realMap = new MemoryBlock*[mapHeigth];
        for (int i = 0; i < mapHeigth; i++){
            realMap[i] = new MemoryBlock[mapWidth];
            for (int j = 0; j < mapWidth; j++){
                realMap[i][j].x = j;
                realMap[i][j].y = i;
            }
        }
        ///inserting start and end position - - - - - - - - - - - - - - - - -
        map[indexOf(start_y, start_x)] = 'S';
        map[indexOf(end_y, end_x)] = 'E';
        ///closing file and returning true as indicator, that everything went well
        file.close();
        return true;
    }
    return false;
}

/**
 * Method used to decide which algorithm does user want to use
 * 
 * @param alg name of searching algoritm
 * @return true if name is correct, false otherwise
 */
bool AreaSearch::useAlgorithm(string alg){
     transform(alg.begin(), alg.end(), alg.begin(), ::tolower);
     if (alg == "dfs")
        dfs();
    else if (alg == "bfs")
        bfs();
    else if(alg == "random")
        randomSearch();
    else if( alg == "dijkstra")
        dijkstraSearch();
    else if( alg == "greedy")
        greedySearch();
    else if (alg == "astar")
        aStarSearch();
    else
        return false;
    printMap();
    printStats();
    cleanRealMap();
    return true;
}

/**
 * Method prints maze map to terminal, it also changes color of some characters
 */
void AreaSearch::printMap(){
    struct winsize w;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &w);
    for(int i = 0; i < w.ws_row; i++)
        printf("\n");
    string output = "";
    for(int i = 0; i<mapLength; i+=mapWidth)
        output += map.substr(i,mapWidth) + "\n";
    /// reseting lost charasters
    output[indexOf(Start->y, Start->x)+Start->y] = 'S';
    output[indexOf(End->y, End->x)+End->y] = 'E';
    ///changing color of charasters in output string
    output = std::regex_replace(output, std::regex("S"), "\033[1;42mS\033[0m");
    output = std::regex_replace(output, std::regex("E"), "\033[1;41mE\033[0m");
    output = std::regex_replace(output, std::regex("o"), "\033[1;32mo\033[0m");
    output = std::regex_replace(output, std::regex("#"), "\033[1;35m#\033[0m");
    cout << "Maze\n" << 
            output << endl;;
}

/**
 * Method prints info about algorithm to terminal, it shows what individual character means, and addition info
 */
void AreaSearch::printStats(){
    cout << "---------------------------------------\n"<<
            "S Start\n" <<
            "E End\n"<<
            "# Opened node\n"<<
            "o Path\n"<<
            "X Wall\n"<<
            "space Fresh node\n"<< 
            "---------------------------------------\n"<<
            "Nodes expanded: " << expadnedNodes<<
            "\nPath length: "<< pathLength << endl;
}

/**
 * Method used to translate 2D position to 1D, where is character of 2D realMap array in map string
 * 
 * @param y vertical position
 * @param x horizontal position
 * @return positon in 1D array, 0<=output<mapLength
 */
int AreaSearch::indexOf(int y, int x){
    return mapWidth*y+x;
}

/**
 * Method used to clean status info in realMap and predecessors
 */
void AreaSearch::cleanRealMap(){
    for (int i = 0; i < mapHeigth; i++){
        for (int j = 0; j < mapWidth; j++){
            realMap[i][j].status = 0;
            realMap[i][j].prev = NULL;
        }
    }
}

/**
 * Method used to find path from current place to first one
 * this method draws path to map variable
 * 
 * @param endM represents end positon of path
 */
void AreaSearch::findPath(MemoryBlock * endM){
    while(endM){
        pathLength++;
        map[indexOf(endM->y, endM->x)] = 'o';//"\032[1;31mot\033[0m\n";
        endM = endM->prev;
    }
    pathLength--;
    map[indexOf(Start->y, Start->x)] = 'S';
    map[indexOf(End->y, End->x)] = 'E'; 
}

/**
 * Method is static, and returns distance from current place to end place using Euclidean heuristic
 * 
 * @param curY current vertical position
 * @param curX current horizontal position
 * @param endY end vertival position
 * @param endX end horizontal position
 * @return distance from curent location to end location
 */
double AreaSearch::heuristicEuclidean(int curY, int curX, int endY, int endX){
    return sqrt((curX - endX)*(curX - endX) + (curY - endY)*(curY - endY));
}

/**
 * Method is static, and returns distance from current place to end place using Manhattan heuristic
 * 
 * @param curY current vertical position
 * @param curX current horizontal position
 * @param endY end vertival position
 * @param endX end horizontal position
 * @return distance from curent location to end location
 */
int AreaSearch::heuristicManhattan(int curY, int curX, int endY, int endX){
    return abs(curX - endX) + abs(curY - endY);
}

/**
 * Method is static, and returns distance from current place to end place using diagonal heuristic
 * 
 * @param curY current vertical position
 * @param curX current horizontal position
 * @param endY end vertival position
 * @param endX end horizontal position
 * @return distance from curent location to end location
 */
int AreaSearch::heuristicDiagonal(int curY, int curX, int endY, int endX){
    return  max(abs(curX - endX), abs(curY - endY)) ;
}

//------------------------------------------------------------------------------------------------------------------------------------
/**
 * Mehod represents BFS algorithm, it uses class variables as input data, and return string filled with path and searched places
 * also returns informations about length of path and number of visited places in maze
 */
void AreaSearch::bfs(){
    MemoryBlock* tmpBlock = NULL, *topN, *downN, *leftN, *rightN;
    list<MemoryBlock*> quote;
    realMap[Start->y][Start->x].status = 1;
    quote.push_back(&realMap[Start->y][Start->x]);
    list<MemoryBlock*>::iterator iter;

    ///while there are unvisited places or end is not reached
    while(!quote.empty()){
        if (pauseMe != 'q'){
            printMap();
            printStats();
            if (pauseMe!='f')
                pauseMe = getchar();
        }         
        tmpBlock = quote.front();
        quote.pop_front();
        if (tmpBlock->y==End->y&&End->x==tmpBlock->x)
            break;
       // this_thread::sleep_for(chrono::milliseconds(500));
        leftN   = &realMap[tmpBlock->y][tmpBlock->x-1];
        rightN  = &realMap[tmpBlock->y][tmpBlock->x+1];
        topN    = &realMap[tmpBlock->y-1][tmpBlock->x];
        downN   = &realMap[tmpBlock->y+1][tmpBlock->x];

        for (auto dirN : {topN, leftN, downN, rightN}){
            if(dirN->y >= 0 && dirN->y < mapHeigth && dirN->x >= 0 && dirN->x < mapWidth && map[indexOf(dirN->y, dirN->x)]!='X'){
                if(realMap[dirN->y][dirN->x].status == 0){
                    realMap[dirN->y][dirN->x].status = 1;
                    realMap[dirN->y][dirN->x].prev = tmpBlock;
                    quote.push_back(&realMap[dirN->y][dirN->x]);
                    map[indexOf(dirN->y, dirN->x)] = '#';
                    expadnedNodes++;
                }
            }
        }
    }
    findPath(tmpBlock);
}
//------------------------------------------------------------------------------------------------------------------------------------
/**
 * Method representing DFS algorithm
 * it fills map variable with places where algorithm were and also it finds path from start to end it there is one
 */
void AreaSearch::dfs(){
    MemoryBlock * tmpBlock = NULL, *topN, *downN, *leftN, *rightN;
    stack<MemoryBlock*> stack;

    stack.push(&realMap[Start->y][Start->x]);
    realMap[Start->y][Start->x].status = 1;
    ///Search while there are unvisited places, or until you reach end
    while(!stack.empty()){
        if (pauseMe != 'q'){
            printMap();
            printStats();
            if (pauseMe!='f')
                pauseMe = getchar();
        }         
        tmpBlock = stack.top();
        stack.pop();

        if (tmpBlock->y==End->y&&End->x==tmpBlock->x)
            break;        
        /// representation of directions
        leftN   = &realMap[tmpBlock->y][tmpBlock->x-1];
        rightN  = &realMap[tmpBlock->y][tmpBlock->x+1];
        topN    = &realMap[tmpBlock->y-1][tmpBlock->x];
        downN   = &realMap[tmpBlock->y+1][tmpBlock->x];

         for (auto dirN : {topN, leftN, downN, rightN}){
            if(dirN->y >= 0 && dirN->y < mapHeigth && dirN->x >= 0 && dirN->x < mapWidth && map[indexOf(dirN->y, dirN->x)]!='X'){
                if(realMap[dirN->y][dirN->x].status == 0){
                    realMap[dirN->y][dirN->x].prev = tmpBlock;
                    realMap[dirN->y][dirN->x].status = 1;
                    map[indexOf(dirN->y, dirN->x)] = '#';
                    expadnedNodes++;
                    stack.push(&realMap[dirN->y][dirN->x]);
                }
            }
        }
    }
    findPath(tmpBlock);
}
//------------------------------------------------------------------------------------------------------------------------------------
/**
 * Method used by randomSearch
 * returns random number from given range
 * 
 * @param range range of number to return from
 * @return number from <0,...,range)
 */
int AreaSearch::randomInt(int range){
    return rand()%range;
}

/**
 * Method representing Random Search algorithm
 * it works like BFS, but instead of visiting places from oldest to new ones this algorithm uses random approach
 */
void AreaSearch::randomSearch(){
    srand(unsigned (time(0)));
    vector <MemoryBlock*> openNodes;
    MemoryBlock * tmpBlock = NULL, *topN, *downN, *leftN, *rightN;
    
    openNodes.emplace_back(&realMap[Start->y][Start->x]); 
    /// While there are unvisited places, or until it reachs end
    while(!openNodes.empty()){
        if (pauseMe != 'q'){
            printMap();
            printStats();
            if (pauseMe != 'f')
                pauseMe = getchar();
        }        
        random_shuffle(openNodes.begin(), openNodes.end(), randomInt);
        tmpBlock = openNodes.back();
        openNodes.pop_back();
        realMap[tmpBlock->y][tmpBlock->x].status = -1;

        if (tmpBlock->y == End->y && tmpBlock->x == End->x)
            break;

        map[indexOf(tmpBlock->y, tmpBlock->x)] = '#';
        /// representation of directions
        leftN   = &realMap[tmpBlock->y][tmpBlock->x-1];
        rightN  = &realMap[tmpBlock->y][tmpBlock->x+1];
        topN    = &realMap[tmpBlock->y-1][tmpBlock->x];
        downN   = &realMap[tmpBlock->y+1][tmpBlock->x];

        for (auto dirN : {topN, leftN, downN, rightN}){
            if(dirN->y >= 0 && dirN->y < mapHeigth && dirN->x >= 0 && dirN->x < mapWidth && map[indexOf(dirN->y, dirN->x)]!='X'){
                if(realMap[dirN->y][dirN->x].status == 0){
                    realMap[dirN->y][dirN->x].status = 1;
                    map[indexOf(dirN->y, dirN->x)] = '#';
                    expadnedNodes++;
                    realMap[dirN->y][dirN->x].prev = tmpBlock;
                    openNodes.emplace_back(&realMap[dirN->y][dirN->x]);
                }
            }
        }
    }
    findPath(tmpBlock);
}
//------------------------------------------------------------------------------------------------------------------------------------
/**
 * Method represents Dijkstra search algorithm
 * all edges has value 1, but it can be easily modified with lambda function
 */
void AreaSearch::dijkstraSearch(){
    int ** lenMap = new int*[mapHeigth];
    MemoryBlock * tmpBlock = NULL, *topN, *downN, *leftN, *rightN;

    for(int i = 0; i< mapHeigth; i++){
        lenMap[i] = new int[mapWidth];
        for (int j = 0; j < mapWidth; j++){
            lenMap[i][j] = INT_MAX;
        } 
    }
    ///comparator for priority queue, it compare values of positions in lenMap array
    auto comp = [&lenMap](MemoryBlock* one, MemoryBlock* two) { return (lenMap[one->y][one->x]) > (lenMap[two->y][two->x]); };
    priority_queue<MemoryBlock*, vector<MemoryBlock*>, decltype(comp)> priQ(comp);

    lenMap[Start->y][Start->x] = 0;
    realMap[Start->y][Start->x].status = 1;
    priQ.push(&realMap[Start->y][Start->x]);

    while(!priQ.empty()){
        if (pauseMe != 'q'){
            printMap();
            printStats();
            if (pauseMe!= 'f')
                pauseMe = getchar();
        }

        tmpBlock = priQ.top();
        priQ.pop();

        if (tmpBlock->y == End->y && tmpBlock->x == End->x)
            break;

        leftN   = &realMap[tmpBlock->y][tmpBlock->x-1];
        rightN  = &realMap[tmpBlock->y][tmpBlock->x+1];
        topN    = &realMap[tmpBlock->y-1][tmpBlock->x];
        downN   = &realMap[tmpBlock->y+1][tmpBlock->x];

        for (auto dirN : {topN, leftN, downN, rightN}){
            if(dirN->y >= 0 && dirN->y < mapHeigth && dirN->x >= 0 && dirN->x < mapWidth &&
                map[indexOf(dirN->y, dirN->x)]!='X'){
                if( realMap[dirN->y][dirN->x].status == 0 &&
                    lenMap[tmpBlock->y][tmpBlock->x] + 1 < lenMap[dirN->y][dirN->x]){
                    lenMap[dirN->y][dirN->x] = lenMap[tmpBlock->y][tmpBlock->x] + 1;
                    realMap[dirN->y][dirN->x].prev = tmpBlock;
                    priQ.push(&realMap[dirN->y][dirN->x]);
                    realMap[dirN->y][dirN->x].status = 1;
                    if(map[indexOf(dirN->y, dirN->x)] != '#')
                        expadnedNodes++;
                    map[indexOf(dirN->y, dirN->x)] = '#';
                }
            }
        }
    }
    findPath(tmpBlock);
}
//------------------------------------------------------------------------------------------------------------------------------------
/**
 * Method represents Greedy searching algorithm
 * it goes straight to end, it also uses euclidean heuristic
 */
void AreaSearch::greedySearch(){
    MemoryBlock * tmpBlock = NULL, *topN, *downN, *leftN, *rightN;
    const int endX = End->x, endY = End->y;
    auto comp = [endX, endY](MemoryBlock* one, MemoryBlock* two) {return (AreaSearch::heuristicEuclidean(one->y, one->x, endY, endX)) > (AreaSearch::heuristicEuclidean(two->y, two->x, endY, endX)); };
    priority_queue<MemoryBlock*, vector<MemoryBlock*>, decltype(comp)> priQ(comp);

    realMap[Start->y][Start->x].status = 1;
    priQ.emplace(&realMap[Start->y][Start->x]);

    while(!priQ.empty()){
        if (pauseMe != 'q'){
            printMap();
            printStats();
            if (pauseMe!= 'f')
                pauseMe = getchar();
        }

        tmpBlock = priQ.top();
        priQ.pop();

        if (tmpBlock->y == End->y && tmpBlock->x == End->x)
            break;

        leftN   = &realMap[tmpBlock->y][tmpBlock->x-1];
        rightN  = &realMap[tmpBlock->y][tmpBlock->x+1];
        topN    = &realMap[tmpBlock->y-1][tmpBlock->x];
        downN   = &realMap[tmpBlock->y+1][tmpBlock->x];

        for (auto dirN : {topN, leftN, downN, rightN}){
            if(tmpBlock->y-1 >= 0 && map[indexOf(dirN->y, dirN->x)]!='X'){
                if( dirN->status == 0 ){
                    map[indexOf(dirN->y, dirN->x)] = '#';
                    expadnedNodes++;
                    dirN->prev = tmpBlock;
                    priQ.emplace(dirN);
                    dirN->status = 1;
                }
            }
        }
    }
    findPath(tmpBlock);
}
//------------------------------------------------------------------------------------------------------------------------------------
/**
 * Method representing A* searching algorithm
 * it uses something from both Dijkstra and Greedy algorithm
 * it also uses Euclidean heuristic, with Manhatton and Diagonal it doesn't work properly
 */
void AreaSearch::aStarSearch(){
    int ** lenMap = new int*[mapHeigth];
    const int endX = End->x, endY = End->y;
    auto comp = [endX, endY, &lenMap](MemoryBlock* one, MemoryBlock* two, int disT = 0) {return (AreaSearch::heuristicEuclidean(one->y, one->x, endY, endX)+lenMap[one->y][one->x]) > (AreaSearch::heuristicEuclidean(two->y, two->x, endY, endX)+lenMap[two->y][two->x]+disT) ; };
    set<MemoryBlock*, decltype(comp)> openSet(comp);
    priority_queue<MemoryBlock*, vector<MemoryBlock*>, decltype(comp)> priQ(comp);
    MemoryBlock * tmpBlock = NULL;
    MemoryBlock * leftN, * rightN, * topN, * downN;

    for(int i = 0; i< mapHeigth; i++){
        lenMap[i] = new int[mapWidth];
        for (int j = 0; j < mapWidth; j++){
            lenMap[i][j] = INT_MAX;
        }
    }
    
    openSet.insert(&realMap[Start->y][Start->x]);
    priQ.emplace(&realMap[Start->y][Start->x]);
    lenMap[Start->y][Start->x] = 0;

    while(!priQ.empty()){
        if (pauseMe != 'q'){
            printMap();
            printStats();
            if (pauseMe != 'f')
                pauseMe = getchar();
        }

        tmpBlock = priQ.top();
        priQ.pop();
        if (tmpBlock->y == End->y && tmpBlock->x == End->x)
            break;

        leftN   = &realMap[tmpBlock->y][tmpBlock->x-1];
        rightN  = &realMap[tmpBlock->y][tmpBlock->x+1];
        topN    = &realMap[tmpBlock->y-1][tmpBlock->x];
        downN   = &realMap[tmpBlock->y+1][tmpBlock->x];

        for (auto dirN : {topN, leftN, downN, rightN}){
            if(dirN->status != -1 && dirN->y >= 0  && dirN->y < mapHeigth && dirN->x >= 0 && dirN->x < mapWidth && map[indexOf(dirN->y, dirN->x)]!='X'){
                if ( dirN->status == 0 || ( comp(tmpBlock, dirN, 1)) ){
                    dirN->prev = tmpBlock;
                    lenMap[dirN->y][dirN->x] = lenMap[tmpBlock->y][tmpBlock->x] + 1;
                    dirN->status = 1;
                    priQ.push(dirN);
                    dirN->status = 1;
                    if(map[indexOf(dirN->y, dirN->x)] != '#')
                        expadnedNodes++;
                    map[indexOf(dirN->y, dirN->x)] = '#';
                }
            }
        }
        tmpBlock->status = -1;
    }
    findPath(tmpBlock);
}

//------------------------------------------------------------------------------------------------------------------------------------
/**
 * On linux run program with this command
 * `./a.out [MAP FILE PATH] {bfs|dfs|random|dijkstra|greedy|astar}`
 * 
 * to run astar algorithm on map located in dataset/0.txt
 * `./a.out dataset/0.txt astar`
 */
int main(int argc, char **argv){
    if (argc != 3){
        cout << "Wrong input" << endl;
        return 1;
    }
    AreaSearch araSea;
    string input = "";
    if (!araSea.setValuesFromFile(argv[1])){
        cout << "File error" << endl;
        return 1;
    }
    if (!araSea.useAlgorithm(argv[2])){
        cout << "Wrong algorithm input" <<endl;
        return 1;
    }
    return 0;
}