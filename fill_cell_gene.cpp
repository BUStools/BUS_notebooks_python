// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppProgress)]]
#include <Rcpp.h>
#include <progress.hpp>
#include <progress_bar.hpp>
#include <fstream>
#include <string>
#include <unordered_map>
#include <iterator>
#include <algorithm>
#include <unordered_set>
using namespace std;
using namespace Rcpp;
// Intersection of gene vectors
//[[Rcpp::export]]
vector<string> gene_intersect(vector<string> a, vector<string> b) {
  vector<string> out;
  sort(a.begin(), a.end());
  sort(b.begin(), b.end());
  set_intersection(a.begin(), a.end(), b.begin(), b.end(), back_inserter(out));
  return out;
}

//' Convert the output of kallisto into gene by cell matrix
//' 
//' This function takes the output file, usually \code{output.sort.txt}, from kallisto,
//' which has 4 columns: barcode, UMI, equivalence class, and counts, and generates
//' the row indices, column pointers, gene UMI count values, accepted barcodes, 
//' and genes (for dim names) to be used in R to construct a sparse matrix.
//' 
//' @param fn File name of the kallisto main output file.
//' @param genes A list with each element a string vector of genes that an equivalence
//' class maps to, generated earlier in R.
//' @param whitelist A string vector showing whiteliisted barcodes from 10x.
//' @param est_ncells Estimated number of cells; providing this argument will
//' speed up computation as it minimizes memory reallocation as vectors grow.
//' @param est_ngenes Estimated number of genes.
//' @param display_progress Whether to display progress.
//' 
//' @return A list with row indices, column pointers, gene UMI count values, 
//' barcodes (for column names), and gene IDs (for row names). You should construct
//' the sparse matrix in R.
//' 
//[[Rcpp::export]]
List fill_cell_gene(const char* fn, List genes, vector<string> whitelist,
                    int est_ncells, int est_ngenes,
                    bool display_progress = true) {
  ifstream infile(fn);
  string bc, umi, ec_str, cts;
  string pbar = "", pumi = "";
  int ec, n, i = 0; // i to keep track of # of iterations
  vector<string> gs, gl;
  unordered_map<string, unordered_map<string, double>> cell_gene;
  cell_gene.reserve(est_ncells);
  // Convert whitelist into unordered_set to speed up lookup
  unordered_set<string> wl(whitelist.begin(), whitelist.end());
  Rcout << "Reading data" << endl;
  while (infile >> bc >> umi >> ec_str >> cts) {
    if (i % 1000 == 0) {
      checkUserInterrupt();
    }
    // If barcode is not in whitelist, skip to the next barcode
    if (wl.find(bc) == wl.end()) {
      continue;
    }
    ec = stoi(ec_str);
    if (bc == pbar) {
      // Same barcode
      if (umi == pumi) {
        // Same umi, get intersection of gene list
        gl = as<vector<string>>(genes[ec]);
        gs = gene_intersect(gs, gl);
      } else {
        // New UMI, process the previous gene list
        n = gs.size();
        for (int j = 0; j < n; j++) {
          cell_gene[bc][gs[j]] += 1.0/(double)n;
        }
        pumi = umi;
        gs = as<vector<string>>(genes[ec]);
      }
    } else {
      // Previous gene list
      n = gs.size();
      for (int j = 0; j < n; j++) {
        cell_gene[pbar][gs[j]] += 1.0/(double)n;
      }
      pumi = umi; pbar = bc;
      gs = as<vector<string>>(genes[ec]);
    }
    // Some sense of progress
    if (display_progress) {
      if (i % 5000000 == 0 && i > 0) {
        Rcout << "Read " << i/1e6 << " million lines" << endl;
      }
    }
    i++;
  }
  // Remember the last gene
  n = gs.size();
  for (int j = 0; j < n; j++) {
    cell_gene[pbar][gs[j]] += 1.0/(double)n;
  }
  
  // Convert the unordered map into a sparse matrix
  // I'm using stl vectors here since they grow more nicely than arma::vec
  vector<string> barcodes, geneIDs;
  barcodes.reserve(est_ncells); geneIDs.reserve(est_ngenes);
  vector<double> values;
  values.reserve(i);
  vector<size_t> rowind, colptr;
  rowind.reserve(i); colptr.reserve(est_ncells + 1);
  unordered_map<string, size_t> rowind_map;
  unordered_map<string, double> g; // for individual genes
  i = 0; // Here keep track of number of iteration in case there's interruption
  size_t entry = 0, gene_row = 0; // keep track of how many entries
  string gn; // gene name
  Rcout << "Constructing sparse matrix" << endl;
  Progress p(cell_gene.size(), display_progress);
  // Consider parallelizing
  for (auto el : cell_gene) {
    if (i % 1000 == 0) {
      checkUserInterrupt();
    }
    barcodes.push_back(el.first);
    colptr.push_back(entry);
    // Construct the row index, iterate through each gene for this barcode
    g = el.second;
    for (auto k: g) {
      gn = k.first;
      if (rowind_map.find(gn) == rowind_map.end()) {
        // If not found
        rowind_map[gn] = gene_row;
        geneIDs.push_back(gn);
        gene_row++;
      }
      rowind.push_back(rowind_map[gn]);
      values.push_back(k.second); // The UMI count
      entry++;
    }
    p.increment();
    i++;
  }
  // Remember the last entry of colptr
  colptr.push_back(entry);
  return List::create(_["rowind"] = rowind, 
                      _["colptr"] = colptr,
                      _["values"] = values,
                      _["barcodes"] = barcodes,
                      _["genes"] = geneIDs);
}