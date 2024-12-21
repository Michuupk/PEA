#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>
#include <cmath>
#include <chrono>
#include <queue>
#include <tuple>
#include <algorithm>
#include <random>

using namespace std;

// Global variables
string filename = "none.atsp";

string myText;
long long graphSize;
vector<vector<long long>> graph;
vector<vector<long long>> graphAdjacency;
vector<long long> bestpath;

// Settings variables
int selection = 0; // 1 - Simulated Annealing, 2 - Genetic Algorithm
long long bestKnown = -1;
long double inittemp = 0;
long double endtemp = 0;
long double cooling = 0;
long long populationSize = 0;
long long generations = 0;
float mutationRate = 0;
int random_init = 0;

long long maxTime = 30; // 30 minutes by default (minutes defined in program)

void loadSettings()
{
    ifstream settingsFile("settings2.txt");
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
        
        if (selection == 1) // for AS
        {
            if (myText.find("RANDOM_INIT") != string::npos)
            {
                istringstream iss(myText);
                string label;
                iss >> label >> random_init;
                
            }
            if (myText.find("INIT_TEMP") != string::npos)
            {
                istringstream iss(myText);
                string label;
                iss >> label >> inittemp;
            }
            if (myText.find("MIN_TEMP") != string::npos)
            {
                istringstream iss(myText);
                string label;
                iss >> label >> endtemp;
            }
            if (myText.find("COOLING_RATE") != string::npos)
            {
                istringstream iss(myText);
                string label;
                iss >> label >> cooling;
            }
        }
        else if (selection == 2) // for Genetic Algorithm
        {
            if (myText.find("POPULATION_SIZE") != string::npos)
            {
                istringstream iss(myText);
                string label;
                iss >> label >> populationSize;
                cout << "Population size: " << populationSize << endl;
            }
            if (myText.find("GENERATIONS") != string::npos)
            {
                istringstream iss(myText);
                string label;
                iss >> label >> generations;
                cout << "Generations: " << generations << endl;
            }
            if (myText.find("MUTATION_RATE") != string::npos)
            {
                istringstream iss(myText);
                string label;
                iss >> label >> mutationRate;
                cout << "Mutation rate: " << mutationRate << endl;
            }
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
        if (myText.find("DIMENSION") != string::npos)
        {
            istringstream iss(myText);
            // cout<<myText<<endl;
            string label;

            iss >> label >> graphSize;
            // cout<< label<<endl;
            // cout << graphSize<<endl;
            // graphFile.close();
        }

        if (myText.find("BEST_KNOWN") != string::npos)
        {
            istringstream iss(myText);
            string label;

            iss >> label >> bestKnown;
            cout << "Best known: " << bestKnown << endl;
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

    for (long long i = 0; i < graphSize; i++) // Creating adjacency list
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

    // for (int i = 0; i < graphSize; i++)     // Printing adjacency list
    // {
    //     cout << "Node " << i << ": ";
    //     for (int j = 0; j < graphAdjacency[i].size(); j++)
    //     {
    //         cout << graphAdjacency[i][j] << " ";
    //     }
    //     cout << endl;
    // }

    graphFile.close();
}

void printPath(const vector<long long> path)
{
    for (auto node : path)
    {
        cout << node << " ";
    }
    cout << endl;
}

long long calculateCost(vector<long long> path, vector<vector<long long>> graph)
{
    long long cost = 0;
    for (int i = 0; i < path.size() - 1; i++)
    {
        if (graph[path[i]][path[i + 1]] <= 0)
        {
            return numeric_limits<long long>::max();
        }
        else
        {
            cost += graph[path[i]][path[i + 1]];
        }
    }
    cost += graph[path[path.size() - 1]][path[0]];
    return cost;
}


void nearestNeighbour(vector<vector<long long>> &graph, long long &graphSize, long long &cost, vector<long long> &bestpath)
{
    cout << "Nearest Neighbour" << endl;

    auto startTimer = chrono::high_resolution_clock::now();
    long double smallestCost = numeric_limits<long double>::max();

    for (long long startNode = 0; startNode < graphSize; startNode++)
    {
        if (chrono::duration_cast<chrono::seconds>(chrono::high_resolution_clock::now() - startTimer).count() > maxTime * 60)
        {
            cout << "Time limit exceeded (30 minutes)" << endl;
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

void generateRandomPath(long long &cost, vector<long long>&bestpath, vector<vector<long long>> &graph, long long &graphSize)
{
    bestpath.reserve(graphSize);
    for (long long i = 0; i < graphSize; i++)
    {
        bestpath.push_back(i);
    }

    random_device rd;
    mt19937 g(rd());
    shuffle(bestpath.begin(), bestpath.end(), g);

    cost = calculateCost(bestpath, graph);
}

vector<long long> generateNeighbour(vector<long long> &path)
{
    vector<long long> newPath = path;
    long long first = rand() % path.size();
    long long second = rand() % path.size();
    while(second == first)
    {
        second = rand() % path.size();
    } 
    swap(newPath[first], newPath[second]);
    return newPath;
}

void SimulatedAnnealing(vector<vector<long long>> &graph, long long &graphSize, long long &cost, vector<long long> &bestpath, long double inittemp, long double &endtemp, long double &cooling)
{
    cout << "Simulated Annealing" << endl;

    auto startTimer = chrono::high_resolution_clock::now();

    while (inittemp > endtemp)
    {
        if (chrono::duration_cast<chrono::seconds>(chrono::high_resolution_clock::now() - startTimer).count() > maxTime * 60)
        {
            cout << "Time limit exceeded (30 minutes)" << endl;
            break;
        }
        vector<long long> newpath;
        long long newcost = 0;
        do
        {
            //iteracje
            newpath = generateNeighbour(bestpath);
            
            newcost = calculateCost(newpath, graph);
        } while (newcost == numeric_limits<long long>::max());
        long long delta = newcost - cost;
        double random = ((double)rand() / (RAND_MAX)); // random number between 0 and 1
        if (delta < 0 || random < exp(-delta / inittemp))
        {
            bestpath = newpath;
            cost = newcost;
        }
        inittemp = inittemp - (inittemp * cooling);
    }
}

void generatePopulation(vector<vector<long long>> &population, long long &populationSize, long long &graphSize)
{
    population.reserve(populationSize);
    for (long long i = 0; i < populationSize; i++)
    {
        vector<long long> path(graphSize);
        iota(path.begin(), path.end(), 0);

        random_device rd;
        mt19937 g(rd());
        shuffle(path.begin(), path.end(), g);
        population.push_back(move(path));
    }
}

void createoffspring(vector<long long> &parent1, vector<long long> &parent2)
{
    int start = 0;
    int end = 0;

    while (start == end && end - start > graphSize - 1)
    {
        start = rand() % (parent1.size());
        end = start + rand() % (parent1.size() - start);
    }
    vector<long long> insertion1;
    vector<long long> insertion2;
    for (int i = start; i < end; i++) // Creating insertion
    {
        insertion1.push_back(parent1[i]);
        insertion2.push_back(parent2[i]);
    }
    int insertpoint = parent1.size();
    insertpoint = rand() % (insertpoint); // Insertion point

    parent2.erase(remove_if(parent2.begin(), parent2.end(), [&insertion1](long long element) { // Removing values of insertion1 from parent2
                      return find(insertion1.begin(), insertion1.end(), element) != insertion1.end();
                  }),
                  parent2.end());

    parent1.erase(remove_if(parent1.begin(), parent1.end(), [&insertion2](long long element) { // Removing values of insertion2 from parent1
                      return find(insertion2.begin(), insertion2.end(), element) != insertion2.end();
                  }),
                  parent1.end());

    parent1.insert(parent1.begin() + insertpoint, insertion2.begin(), insertion2.end()); // Inserting insertion to parent
    parent2.insert(parent2.begin() + insertpoint, insertion1.begin(), insertion1.end()); // Inserting insertion to parent

    insertion1.clear();
    insertion2.clear();
}

void mutateoffspring(vector<vector<long long>> &offspring, long long &graphSize, float &mutationRate)
{
    for (int i = 0; i < offspring.size(); i++)
    {
        double random = ((double)rand() / (RAND_MAX)); // random number between 0 and 1
        if (random < mutationRate)
        {
            int first = rand() % offspring[i].size();
            if (first == offspring[i].size() - 1)
            {
                first--;
            }
            int second = first + 1;
            swap(offspring[i][first], offspring[i][second]);
        }
    }
}

void GeneticAlgorithm(vector<vector<long long>> &graph, long long &graphSize, long long &cost, vector<long long> &bestpath, long long &populationSize, long long &generations, float &mutationRate)
{
    cout << "Genetic Algorithm" << endl;

    vector<vector<long long>> population;
    generatePopulation(population, populationSize, graphSize);

    for (long long i = 0; i < generations; i++)
    {
        vector<long long> fitness;
        int pathCost = 0;
        for (auto path : population)
        {
            pathCost = calculateCost(path, graph);
            fitness.push_back(pathCost);
        }

        vector<vector<long long>> newPopulation; // survivors
        int midway = population.size() / 2;
        for (long long j = 0; j < midway; j++)
        {
            if (fitness[j] == numeric_limits<long long>::max() && fitness[j + midway] == numeric_limits<long long>::max()) // if both paths are invalid
            {
                continue;
            }
            else if (fitness[j] == fitness[j + midway])
            {
                newPopulation.push_back(population[j]);
                newPopulation.push_back(population[j + midway]);
            }
            else if (fitness[j] < fitness[j + midway])
            {
                newPopulation.push_back(population[j]);
            }
            else
            {
                newPopulation.push_back(population[j + midway]);
            }
        }
        fitness.clear();
        population.swap(newPopulation); // survivors become new population
        newPopulation.clear();
        vector<long long> parent1;
        vector<long long> parent2;
        vector<vector<long long>> offspring;
        while (offspring.size() < populationSize + 1) // creating offspring
        {
            long long random = rand() % (populationSize / 2);
            parent1 = population[random]; // random parent selection
            int stop = 0;
            do
            {
                random = rand() % population.size();
                parent2 = population[random]; // random parent selection
                stop++;
            } while (parent2 == parent1 && stop < 100);
            createoffspring(parent1, parent2);

            offspring.push_back(parent1);
            offspring.push_back(parent2);
            parent1.clear();
            parent2.clear();
        }
        mutateoffspring(offspring, graphSize, mutationRate); // mutation

        population.swap(offspring);
        population.resize(populationSize);

        offspring.clear();

        for (auto path : population)
        {
            pathCost = calculateCost(path, graph);
            if (pathCost < cost)
            {
                cost = pathCost;
                bestpath = path;
            }
        }
    }
}

void validPath(vector<long long> path) // function checking if every city is visited only once
{
    if (path.size() != graphSize)
    {
        cout << endl
             << "Path is invalid" << endl;
        return;
    }
    for (int i = 0; i < path.size(); i++)
    {
        for (int j = i + 1; j < path.size(); j++)
        {
            if (path[i] == path[j])
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

    loadSettings();      // loading settings from file
    LoadData(graphSize); // loading graph size
    cout << "Loaded graph has " << graphSize << " Nodes" << endl;

    LoadGraph(graph, graphAdjacency); // loading graph to memory
    cout << "Graph loaded to memory" << endl;

    long long cost = 0;
    long long topLimit = 0;


    if (random_init == 1 )
    {
        generateRandomPath(cost, bestpath, graph, graphSize);
    }
    else nearestNeighbour(graph, graphSize, cost, bestpath); // calling nearest neighbour method for top limit

    chrono::time_point<std::chrono::high_resolution_clock> start_clock, end_clock; // variables for time measurement

    if (selection == 1)
    {
        start_clock = chrono::high_resolution_clock::now();
        SimulatedAnnealing(graph, graphSize, cost, bestpath, inittemp, endtemp, cooling);
        end_clock = chrono::high_resolution_clock::now();
        chrono::duration<double> result = end_clock - start_clock;
        cout << "Time of execution: " << result.count() << "s" << endl;
        if (cost != numeric_limits<long long>::max())
        {
            cout << "Lowest cost: " << cost << endl;
            if (bestpath.size() < 20)
            {
                cout << "Best path: ";
                printPath(bestpath);
            }
        }
        else
        {
            cout << "No path found" << endl;
        }
    }
    if (selection == 2)
    {
        start_clock = chrono::high_resolution_clock::now();
        GeneticAlgorithm(graph, graphSize, cost, bestpath, populationSize, generations, mutationRate);
        end_clock = chrono::high_resolution_clock::now();
        chrono::duration<double> result = end_clock - start_clock;
        cout << "Time of execution: " << result.count() << "s" << endl;
        if (cost != numeric_limits<long long>::max())
        {
            cout << "Lowest cost: " << cost << endl;
            if (bestpath.size() < 20)
            {
                cout << "Best path: ";
                printPath(bestpath);
            }
        }
        else
        {
            cout << "No path found" << endl;
        }
    }

    validPath(bestpath); // checking if path is valid
    cout << endl
         << "End of program" << endl;
    //cin.get();
    return 0;
}