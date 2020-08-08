
#include <iostream>
#include "iss.h"

//check pair of matrix element which can be used to define (and replace) searchedEdge matrix element in Ideally consistent PCM
bool findApropriatePair(Edge searchedEdge, int N, const vector<Edge> &subset, pair<Edge, Edge> &edgePair, bool &inverse)
{
    inverse = false;

    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            for (int k = 0; k < N; k++)
            {
                for (int m = k + 1; m < N; m++)
                {
                    bool inSameCol = (i == searchedEdge.src && k == searchedEdge.dest && j == m);
                    bool inSameRow = (j == searchedEdge.src && m == searchedEdge.dest && i == k);
                    inverse = (i == searchedEdge.src && m == searchedEdge.dest && j == k) || (j == searchedEdge.src && k == searchedEdge.dest && i == m);
                    
                    if (inSameRow || inSameCol || inverse)
                    {
                        Edge edge1 = { i,j };
                        Edge edge2 = { k,m };

                        if (edge1 != edge2 && edge1 != searchedEdge && searchedEdge != edge2)
                        {
                            if (std::find(subset.begin(), subset.end(), edge1) != subset.end() &&
                                std::find(subset.begin(), subset.end(), edge2) != subset.end())
                            {
                                if (inSameRow && !inverse)
                                    edgePair = make_pair(edge2, edge1);
                                else
                                    edgePair = make_pair(edge1, edge2);
                                return true;
                            }
                        }
                    }
                }
            }
        }
    }

    return false;
}

//function to calculate priority vector
vector<double> priorityVector(const vector<vector<double>> &pcm, int N, bool additive)
{
    vector<double> finalPriorityVec;
    vector<vector<double>> priorityVectors;

    //choose possible sets from completeS (one set size N-1)
    vector<vector<Edge>> iss = informativelySignificantSubsets(N);

    //completeS - all matrix elements in top triangle, above main diagonal (size N*(N-1)/2)
    vector<Edge> completeS;
    for (int i = 0; i < N; i++)
    {
        for (int j = i + 1; j < N; j++)
        {
            completeS.push_back({ i,j });
        }
    }

    for (auto &subset : iss)
    {
        std::vector<Edge>::iterator it;
        auto subsetCopy = subset;
        auto pcmCopy = pcm;
        vector<Edge> unknown(completeS.size());

        //find which points from Ideally consistent PCM should be defined (replaced)
        it = std::set_symmetric_difference(completeS.begin(), completeS.end(), subset.begin(), subset.end(), unknown.begin());
        unknown.resize(it - unknown.begin());

        if (unknown.size() != completeS.size() - subset.size())
        {
            continue;
        }

        int currEl = 0;
        //all points from Ideally consistent PCM should be defined - unknown.empty()
        while (!unknown.empty() && currEl < unknown.size())
        {
            Edge edge = unknown[currEl];

            bool inverse;
            pair<Edge, Edge> edgePair;
            //replace edge with pair of matrix element (from iss)
            if (findApropriatePair(edge, N, subsetCopy, edgePair, inverse))
            {
                double first = pcmCopy[edgePair.first.src][edgePair.first.dest];
                double second = pcmCopy[edgePair.second.src][edgePair.second.dest];
                second = (inverse ? (additive ? -second : 1. / second) : second);

                pcmCopy[edge.src][edge.dest] = (additive ? first - second : first / second);
                unknown.erase(unknown.begin() + currEl);
                subsetCopy.push_back(edge);
                currEl = 0;
            }
            else
            {
                currEl++;
            }
        }

        if (currEl >= unknown.size() && currEl>0)
            continue;

        double summ = 0;
        vector<double> priority;
        for (int i = 0; i < N; i++)
        {
            summ += pcmCopy[i][N - 1];
        }

        for (int i = 0; i < N; i++)
        {
            priority.push_back(pcmCopy[i][N - 1] / summ);
        }

        priorityVectors.push_back(priority);
    }

    //calculate average (mean) priority vector
    for (int i = 0; i < N; i++)
    {
        double summ = 0;
        for (auto &priority : priorityVectors)
        {
            summ += priority[i];
        }

        finalPriorityVec.push_back(summ / priorityVectors.size());
    }

    return finalPriorityVec;
}

int main()
{
    bool additive = false; //if false f(mi,mj) - multiplicative
    //size of PCM
    int N = 4;
    double diagonalElem = (additive ? 0. : 1.);

    //pairwise comparison matrix (PCM)
    vector<vector<double>> pcm = { {diagonalElem, 7, 6, 5}
                                  ,{0, diagonalElem, 1. / 2., 1}
                                  ,{0, 0, diagonalElem, 1. / 2.}
                                  ,{0, 0, 0, diagonalElem} };

    // PCM is antisymmetric
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < i; j++)
        {
            pcm[i][j] = (additive ? -pcm[j][i] : 1. / pcm[j][i]);
        }
    }

    //calculate priority vector
    vector<double> priority = priorityVector(pcm, N, additive);

    cout << "Priority vector: ";
    for (int i = 0; i < N; i++)
        cout << priority[i] << "\t";
}
