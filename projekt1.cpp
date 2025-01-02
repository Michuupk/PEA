#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <cmath>
#include <chrono>

using namespace std;

// Global variables
string filename = "none.atsp";

string myText;
long long graphSize;
int bestKnown = -1;
long long maxTime = 30; // 30 minutes by default (minutes defined in program)
vector<vector<long long>> graph;
vector<vector<long long>> graphAdjacency;
vector<long long> bestpath;

// Settings variables
int selection = 0; // 1 - Random TSP, 2 - Nearest Neighbour, 3 - Brute Force
long long howManyRandomPaths = 0;
void loadSettings()
{
    ifstream settingsFile("settings1.txt");
    if (!settingsFile.is_open())
    {
        cerr << "Error: Could not open the file!" << endl;
        return;
    }

    while (getline(settingsFile, myText))
    {
        if (myText.find("NAME_OF") != string::npos)
        {
            istringstream iss(myText);
            string label;
            iss >> label >> filename;
        }
        if (myText.find("SELECTION") != string::npos)
        {
            istringstream iss(myText);
            string label;
            iss >> label >> selection;
        }
        if (myText.find("BEST_KNOWN") != string::npos)
        {
            istringstream iss(myText);
            string label;
            iss >> label >> bestKnown;
            cout << "Best known: " << bestKnown << endl;
        }
        if (myText.find("TIME_LIMIT") != string::npos)
        {
            istringstream iss(myText);
            string label;
            iss >> label >> maxTime;
            break;
        }
    }
    settingsFile.close();
}

void LoadData(long long &graphSize) // Checking the size of graph
{

    ifstream graphFile(filename);
    if (!graphFile.is_open())
    {
        cerr << "Error: Could not open the file!" << endl;
        return;
    }

    while (getline(graphFile, myText))
    {
        if (myText.find("BEST_KNOWN") != string::npos)
        {
            istringstream iss(myText);
            string label;

            iss >> label >> bestKnown;
            cout << "Best known: " << bestKnown << endl;
            break;
        }

        if (myText.find("DIMENSION") != string::npos)
        {
            istringstream iss(myText);
            // cout<<myText<<endl;
            string label;

            iss >> label >> graphSize;
            // cout<< label<<endl;
            // cout << graphSize<<endl;
            // graphFile.close();
            break;
        }
    }
}

void LoadGraph(vector<vector<long long>> &graph, vector<vector<long long>> &graphAdjacency) // Loading graph to memory
{
    ifstream graphFile(filename);
    if (!graphFile.is_open())
    {
        cerr << "Error: Could not open the file!" << endl;
        return;
    }

    while (getline(graphFile, myText))
    {
        if (myText.find("EDGE_WEIGHT_SECTION") != string::npos)
        {

            for (long long i = 0; i < graphSize; i++)
            {
                vector<long long> row;

                for (long long j = 0; j < graphSize; j++)
                {
                    long long value;
                    graphFile >> value;
                    row.push_back(value);
                }
                graph.push_back(row);
            }
        }
        else
        {
            continue;
        }
    }

    // for (int i = 0; i < graphSize; i++)     // Printing graph
    // {
    //     for (int j = 0; j < graphSize; j++)
    //     {
    //         cout << graph[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    for (long long i = 0; i < graphSize; i++) // Creating adjacency matrix
    {
        vector<long long> row;
        // cout << "Node " << i << ": ";
        for (long long j = 0; j < graphSize; j++)
        {
            if (i != j && graph[i][j] > 0)
            {
                row.push_back(j);
                // cout << j << " ";
            }
        }
        // cout << endl;
        graphAdjacency.push_back(row);
    }

    // for (int i = 0; i < graphSize; i++)     // Printing adjacency matrix
    // {
    //     cout << "Node " << i << ": ";
    //     for (int j = 0; j < graphAdjacency[i].size(); j++)
    //     {
    //         cout << graphAdjacency[i][j] << " ";
    //     }
    //     cout << endl;
    //}

    graphFile.close();
}

void randomTSP(vector<vector<long long>> &graph, long long &graphSize, long double &cost, vector<long long> &bestpath)
{
    cout << "Random TSP" << endl;

    auto startTimer = chrono::high_resolution_clock::now();

    long double smallestCost = numeric_limits<long double>::max();
    vector<vector<long long>> allPaths;
    vector<long long> path;
    long long current = 0;
    int bestPathID = 0;

    long long visited[graphSize];
    for (long long i = 0; i < graphSize; i++) // filling with zeros
    {
        visited[i] = 0;
    }
    for (long long j = 0; j < numeric_limits<int>::max(); j++) // creating multiple paths, not enough for optimal result, big enough for my PC to slow down
    {
        auto startPath = chrono::high_resolution_clock::now();
        if (chrono::duration_cast<chrono::seconds>(chrono::high_resolution_clock::now() - startTimer).count() > maxTime * 60)
        {
            cout << "Time limit exceeded "<<maxTime<<" minutes" << endl;
            break;
        }

        long long k = j;
        long long limitcheck = 0;
        for (long long i = 0; i < graphSize - 1; i++) // creating random path
        {
            if (i == 0)
            {
                current = rand() % graphSize; // random starting polong long
                visited[current] += 1;
                path.push_back(current);
            }
            limitcheck = 0;
            long long next = rand() % graphSize;
            while (visited[next] == j + 1 || next == current) //|| graph[current][next] > 0) // checking if node was visited or is the same as current
            {
                next = rand() % graphSize;
                limitcheck++;
                if (limitcheck > 1000) // if limit is reached, break the loop
                {
                    break;
                }
            }
            if (limitcheck > 1000) // path is invalid, clear it
            {
                path.clear();
                for (int z = 0; z < graphSize; z++)
                {
                    visited[z] = j + 1;
                }
                break;
            }
            visited[next] += 1;
            path.push_back(next);
            cost += graph[current][next];
            current = next;
        }
        if (limitcheck < 1000)
        {
            path.push_back(path[0]);
            cost += graph[current][path[0]];
            allPaths.push_back(path);

            for (int i = 0; i < (allPaths.size()) - 1; i++)
            {
                if (allPaths[i] == path)
                {
                    allPaths.pop_back();
                    j--;
                }
            }
            path.clear();
            if (cost < smallestCost)
            {
                smallestCost = cost;
                bestPathID = k;
            }
            cost = 0;
        }
    }
    if (allPaths.size() == 0)
    {
        cout << "No valid paths found" << endl;
        bestpath.clear();
        bestpath.push_back(0);
        return;
    }
    cost = smallestCost;
    if (cost == bestKnown && bestKnown != -1)
    {
        cout << "Found optimal =" << cost << endl;
    }

    for (long long i = 0; i < graphSize; i++)
    {
        bestpath.push_back(allPaths[bestPathID][i]);
    }
}

void printPath(const vector<long long> &path){
    for (auto node : path)
    {
        cout << node << " ";
    }
    cout << endl;
}

void nearestNeighbour(vector<vector<long long>> &graph, long long &graphSize, long double &cost, vector<long long> &bestpath)
{
    cout << "Nearest Neighbour" << endl;

    auto startTimer = chrono::high_resolution_clock::now();
    long double smallestCost = numeric_limits<long double>::max();

    for (long long startNode = 0; startNode < graphSize; startNode++)
    {
        if (chrono::duration_cast<chrono::seconds>(chrono::high_resolution_clock::now() - startTimer).count() > maxTime * 60)
        {
            cout << "Time limit exceeded "<<maxTime<<" minutes" << endl;
            break;
        }

        vector<long long> path;
        vector<long long> visited(graphSize, 0); // array for visited cities
        long long current = startNode;
        long double currentCost = 0;

        visited[current] = 1;
        path.push_back(current);

        for (long long i = 1; i < graphSize; i++) // creating path
        {
            long long next = -1;
            long long smallest = numeric_limits<long long>::max();
            for (long long j = 0; j < graphSize; j++) // finding nearest neighbour
            {
                if (graph[current][j] < smallest && graph[current][j] > 0 && visited[j] == 0)
                {
                    smallest = graph[current][j];
                    next = j;
                }
            }
            if (next == -1)
                break; // no valid next node found
            visited[next] = 1;
            path.push_back(next);
            currentCost += graph[current][next];
            current = next;
        }
        
        currentCost += graph[current][startNode]; // Returning to the start point
        
        if (currentCost < smallestCost && path.size() == graphSize)
        {
            smallestCost = currentCost;
            bestpath = path;
        }
        visited.clear();
    }

    cost = smallestCost;

    if (cost == bestKnown && bestKnown != -1)
    {
        cout << "Found optimal =" << cost << endl;
    }
}

void dfs(vector<vector<long long>> &graph, vector<long long> &visited, vector<long long> &path, long long currentNode, long long &graphSize, long double currentCost, long double &smallestCost, vector<long long> &bestpath, const chrono::time_point<chrono::high_resolution_clock> &startTimer, long double bestknown, bool &foundBestPath)
{
    if (foundBestPath)
        return;

    auto currentTime = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::seconds>(currentTime - startTimer).count();
    if (duration >= maxTime * 60)
    { // 30 minutes
        return;
    }

    // cout << ", Current Cost: " << currentCost << endl;

    if (path.size() == graphSize) // If we have visited all nodes, check the total cost (including return to start)
    {
        currentCost += graph[currentNode][path[0]]; // Returning to the start point
        if (currentCost < smallestCost)
        {
            smallestCost = currentCost;
            bestpath = path;
            // cout << "New Best Path Found: ";
            // for (auto node : bestpath) {
            //     cout << node << " ";
            // }
            // cout << ", Cost: " << smallestCost << endl;
            if (smallestCost == bestknown)
            {
                foundBestPath = true;
            }
        }
        return;
    }

    // DFS on all unvisited nodes
    for (long long nextNode = 0; nextNode < graphSize; nextNode++)
    {
        if (foundBestPath)
            return;
        if (!visited[nextNode] && graph[currentNode][nextNode] > 0)
        {
            visited[nextNode] = 1;
            path.push_back(nextNode);
            dfs(graph, visited, path, nextNode, graphSize, currentCost + graph[currentNode][nextNode], smallestCost, bestpath, startTimer, bestknown, foundBestPath);
            // Backtracking
            visited[nextNode] = 0;
            path.pop_back();
        }
    }
}

void bruteForce(vector<vector<long long>> &graph, long long &graphSize, long double &cost, vector<long long> &bestpath, long double bestknown)
{
    cout << "Brute Force using DFS from every node" << endl;
    auto startTimer = chrono::high_resolution_clock::now();

    // Initialize variables
    long double smallestCost = numeric_limits<long double>::max();
    vector<long long> overallBestPath;
    bool foundBestPath = false;

    for (long long startNode = 0; startNode < graphSize; startNode++)
    {
        if (foundBestPath)
            break;
        vector<long long> visited(graphSize, 0);
        vector<long long> path;
        path.push_back(startNode);
        visited[startNode] = 1;
        dfs(graph, visited, path, startNode, graphSize, 0, smallestCost, overallBestPath, startTimer, bestknown, foundBestPath);
    }

    bestpath = overallBestPath;
    cost = smallestCost;

    if (graphSize < 16)
    {
        cout << "Overall Best Path: ";
        for (auto node : bestpath)
        {
            cout << node << " ";
        }
    }
    cout << ", Cost: " << cost << endl;
}

void validPath(vector<long long> &bestpath) // function checking if every city is visited only once
{
    if (bestpath.size() != graphSize)
    {
        cout << endl
             << "Path is invalid" << endl;
        return;
    }
    for (int i = 0; i < bestpath.size(); i++)
    {
        for (int j = i + 1; j < bestpath.size(); j++)
        {
            if (bestpath[i] == bestpath[j])
            {
                cout << endl
                     << "Path is invalid" << endl;
                return;
            }
        }
    }
    cout << endl
         << "Path is valid" << endl;
}

int main()
{
    srand(time(NULL)); // seed for random number generator

    loadSettings(); // loading settings from file
    LoadData(graphSize);
    // graphSize = 4; // for testing purposes
    cout << "Loaded graph has " << graphSize << " Nodes" << endl;

    LoadGraph(graph, graphAdjacency); // loading graph to memory
    cout << "Graph loaded to memory" << endl;

    long double cost = 0.0;
    chrono::time_point<std::chrono::high_resolution_clock> start_clock, end_clock; // variables for time measurement

    if (selection == 1) // selection of method depending on settings
    {
        start_clock = chrono::high_resolution_clock::now();
        randomTSP(graph, graphSize, cost, bestpath);
        end_clock = chrono::high_resolution_clock::now();
        chrono::duration<double> result = end_clock - start_clock;
        cout << "Time of execution: " << result.count() << "s" << endl;
        cout << "Lowest cost: " << cost << endl;
        cout << "Best path: ";
        for (long long i = 0; i < bestpath.size(); i++)
        {
            cout << bestpath[i] << "->";
        }
    }
    else if (selection == 2)
    {
        start_clock = chrono::high_resolution_clock::now();
        nearestNeighbour(graph, graphSize, cost, bestpath);
        end_clock = chrono::high_resolution_clock::now();
        chrono::duration<double> result = end_clock - start_clock;
        cout << "Time of execution: " << result.count() << "s" << endl;
        cout << "Lowest cost: " << cost << endl;
        cout << "Best path: ";
        for (long long i = 0; i < bestpath.size(); i++)
        {
            cout << bestpath[i] << "->";
        }
    }
    else if (selection == 3)
    {
        start_clock = chrono::high_resolution_clock::now();
        bruteForce(graph, graphSize, cost, bestpath, bestKnown);
        end_clock = chrono::high_resolution_clock::now();
        chrono::duration<double> result = end_clock - start_clock;
        cout << "Time of execution: " << result.count() << "s" << endl;
        cout << "Lowest cost: " << cost << endl;
        cout << "Best path: ";
        for (long long i = 0; i < bestpath.size(); i++)
        {
            cout << bestpath[i] << "->";
        }
    }

    validPath(bestpath);

    cout << endl
         << "End of program" << endl;
    cin.get();
    return 0;
}