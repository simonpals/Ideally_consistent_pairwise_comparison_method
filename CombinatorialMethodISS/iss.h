#pragma once

#include <iostream>
#include <vector>
#include <algorithm>
#include <map>
#include <queue>

using namespace std;

struct Edge;

vector<int> stack;
vector<vector<Edge>> spaningTrees;
vector<Edge> edges;
map<int, bool> vertixCounter;

// data structure to store graph edges
//which represent m[i,j] where i - src, j - dest
struct Edge {
    int src, dest;

    bool operator==(const Edge& edge) const
    {
        return (this->src == edge.src && this->dest == edge.dest);
    }

    bool operator!=(const Edge& edge) const
    {
        return !(*this == edge);
    }

    bool operator < (const Edge& edge) const
    {
        if (this->src == edge.src)
            return (this->dest < edge.dest);
        return this->src < edge.src;
    }
};

// class to represent a graph object
class Graph
{
public:
    int size = 0;
    // construct a vector of vectors to represent an adjacency list
    vector<vector<int>> adjList;

    // Graph Constructor
    Graph(vector<Edge> const &edges, int N)
    {
        size = N;
        // resize the vector to N elements of type vector<int>
        adjList.resize(N);

        // add edges to the undirected graph
        for (auto &edge : edges)
        {
            adjList[edge.src].push_back(edge.dest);
            adjList[edge.dest].push_back(edge.src);
        }
    }
};

//in m[i,j] we need i < j
Edge normalizeDirection(Edge edge)
{
    if (edge.src > edge.dest)
        return { edge.dest, edge.src };
    return edge;
}

//check stack content if there is new subset to choose
void fillEdgesFromStack(const Graph &graph, const vector<int> &stack)
{
    edges.clear();
    vertixCounter.clear();

    if (stack.size() >= graph.size)
    {
        for (int j = 1; j < stack.size(); j++)
        {
            Edge e = normalizeDirection({ stack[j - 1], stack[j] });
            if (std::find(edges.begin(), edges.end(), e) == edges.end())
            {
                vertixCounter[stack[j - 1]] = true;
                vertixCounter[stack[j]] = true;
                edges.push_back(e);
            }

            if (edges.size() == graph.size - 1)
                break;
        }
    }

    if (edges.size() == graph.size - 1 && vertixCounter.size() == graph.size)
    {
        sort(edges.begin(), edges.end());

        if (std::find(spaningTrees.begin(), spaningTrees.end(), edges) == spaningTrees.end())
            spaningTrees.push_back(edges);
    }
}


// Depth First Search (DFS) Recursive Implementation
// Function to perform DFS Traversal
void DFS(Graph const &graph, int v, int level)
{
    if (spaningTrees.size() == pow(graph.size, graph.size - 2))
        return;

    stack.push_back(v);

    if (level == (graph.size - 1) * 2) // graph.size)  
    {
        fillEdgesFromStack(graph, stack);
        stack.pop_back();
        return;
    }

    // do for every edge (v -> u)
    for (int u : graph.adjList[v])
    {
        DFS(graph, u, level + 1);
    }

    stack.pop_back();
}

//// Perform BFS on graph starting from vertex v
//void BFS(Graph const &graph, int v)
//{
//    // create a queue used to do BFS
//    vector<int> q;
//    map<int, bool> counter;
//
//    // push source vertex into the queue
//    q.push_back(v);
//    q.push_back(v);
//
//    // run till queue is not empty
//    while (!q.empty())
//    {
//        // pop front node from queue and print it
//        q.erase(q.begin());
//        v = q.front();
//        if(q.size() > 0) q.erase(q.begin());
//        int visitedVert = counter.size();
//        counter[v] = true;
//
//        if (visitedVert < graph.size)
//        {
//            // do for every edge (v -> u)
//            for (int u : graph.adjList[v])
//            {
//                q.push_back(v);
//                q.push_back(u);
//                fillEdgesFromStack(graph, q);
//            }
//        }
//        else
//        {
//            fillEdgesFromStack(graph, q);
//        }
//    }
//}

vector<vector<Edge>> informativelySignificantSubsets(int N)
{
    //N - Number of nodes in the graph

    vector<Edge> edges;
    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            edges.push_back({ i,j });
        }
    }

    //informatively significant subsets
    vector<vector<Edge>> iss;
    // create a graph from given edges
    Graph graph(edges, N);

    // Do DFS traversal from all undiscovered nodes to
    // cover all unconnected components of graph
    for (int i = 0; i < N; i++)
    {
        DFS(graph, i, 1);
    }

    //// Do BFS traversal to cover all unconnected components of graph
    //for (int i = 0; i < N; i++) {
    //    // start BFS traversal from vertex i
    //    BFS(graph, i);
    //}

    return spaningTrees;
}