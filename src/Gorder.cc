#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <set>
#include <functional>
#include <climits>
#include <ctime>
#include <stdlib.h>
#include <cstdlib>
#include <cstdint>
#include <cstring>
#include <chrono>
#include <sys/time.h>
#include <algorithm>

#include "GorderGraph.h"
#include "GorderUtil.h"

#include "benchmark.h"
#include "builder.h"
#include "command_line.h"
#include "graph.h"
#include "pvector.h"


using namespace std;

const int INPUTNUM=1;

typedef float ScoreT;
const float kDamp = 0.85;

pvector<ScoreT> PageRankPull(const Graph &g, int max_iters,
                             double epsilon = 0) {
  const ScoreT init_score = 1.0f / g.num_nodes();
  const ScoreT base_score = (1.0f - kDamp) / g.num_nodes();
  pvector<ScoreT> scores(g.num_nodes(), init_score);
  pvector<ScoreT> outgoing_contrib(g.num_nodes());
  for (int iter=0; iter < max_iters; iter++) {
    double error = 0;
    #pragma omp parallel for
    for (NodeID n=0; n < g.num_nodes(); n++)
      outgoing_contrib[n] = scores[n] / g.out_degree(n);
    #pragma omp parallel for reduction(+ : error) schedule(dynamic, 64)
    for (NodeID u=0; u < g.num_nodes(); u++) {
      ScoreT incoming_total = 0;
      for (NodeID v : g.in_neigh(u))
        incoming_total += outgoing_contrib[v];
      ScoreT old_score = scores[u];
      scores[u] = base_score + kDamp * incoming_total;
      error += fabs(scores[u] - old_score);
    }
    //printf(" %2d    %lf\n", iter, error);
    if (error < epsilon)
      break;
  }
  return scores;
}


void PrintTopScores(const Graph &g, const pvector<ScoreT> &scores) {
  vector<pair<NodeID, ScoreT>> score_pairs(g.num_nodes());
  for (NodeID n=0; n < g.num_nodes(); n++) {
    score_pairs[n] = make_pair(n, scores[n]);
  }
  int k = 5;
  vector<pair<ScoreT, NodeID>> top_k = TopK(score_pairs, k);
  k = min(k, static_cast<int>(top_k.size()));
  for (auto kvp : top_k)
    cout << kvp.second << ":" << kvp.first << endl;
}


// Verifies by asserting a single serial iteration in push direction has
//   error < target_error
bool PRVerifier(const Graph &g, const pvector<ScoreT> &scores,
                        double target_error) {
  const ScoreT base_score = (1.0f - kDamp) / g.num_nodes();
  pvector<ScoreT> incomming_sums(g.num_nodes(), 0);
  double error = 0;
  for (NodeID u : g.vertices()) {
    ScoreT outgoing_contrib = scores[u] / g.out_degree(u);
    for (NodeID v : g.out_neigh(u))
      incomming_sums[v] += outgoing_contrib;
  }
  for (NodeID n : g.vertices()) {
    error += fabs(base_score + kDamp * incomming_sums[n] - scores[n]);
    incomming_sums[n] = 0;
  }
  PrintTime("Total Error", error);
  return error < target_error;
}

int main(int argc, char* argv[]){

  // Relabeling graph and running PR 
  CLPageRank cli(argc, argv, "pagerank", 1e-4, 20);
  if (!cli.ParseArgs())
    return -1;
  Builder b(cli);
  Graph g_pr = b.MakeGraph();

  auto PRBound = [&cli] (const Graph &g_pr) {
    return PageRankPull(g_pr, cli.max_iters(), cli.tolerance());
  };
  auto VerifierBound = [&cli] (const Graph &g_pr, const pvector<ScoreT> &scores) {
    return PRVerifier(g_pr, scores, cli.tolerance());
  };

    int W=5;
    int k = cli.relabel_top_k();
    bool undirected=cli.symmetrize();
    clock_t start, end;
    string filename = cli.filename();

    srand(time(0));
    GorderGraph g_tmp;
    string name;
    name=extractFilename(filename.c_str());
    g_tmp.setFilename(name);

    //start=clock();
    vector<int> removed;
    g_tmp.readGorderGraph(filename, undirected);
    removed = g_tmp.RemoveGreaterThanTopK(k);
    cout << "Vertices removed: " << removed.size() << endl; // should be zero if k <= 0

    cout << name << " readGorderGraph is complete." << endl;
    //end=clock();
    //cout << "Time Cost: " << (double)(end-start)/CLOCKS_PER_SEC << endl;

    start=clock();
    vector<int> order;
    g_tmp.GorderGreedy(order, removed, W);
    end=clock();
    cout << "ReOrdered Time Cost: " << (double)(end-start)/CLOCKS_PER_SEC << endl;
  
  g_pr = Builder::Relabel(g_pr, order);
  BenchmarkKernel(cli, g_pr, PRBound, PrintTopScores, VerifierBound);
  return 0;

}

