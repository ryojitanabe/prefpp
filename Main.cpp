// Copyright (C) 2015-17 Andrzej Jaszkiewicz

#include "problem.h"
#include "quadtree.h"
#include "ttreeset.h"
#include "mfront.h"
#include "mfront2.h"
#include <sstream>
#include <chrono>

#include <iostream>
#include <vector>
#include <numeric>
#include <random>

TProblem Problem;
vector <TPoint*> NDS_solutions;

/**
 * Argsort(currently support ascending sort)
 * Reference https://gist.github.com/HViktorTsoi/58eabb4f7c5a303ced400bcfa816f6f5
 * @tparam T array element type
 * @param array input array
 * @return indices w.r.t sorted array
 */
template<typename T>
std::vector<size_t> argsort(const std::vector<T> &array) {
  std::vector<size_t> indices(array.size());
  std::iota(indices.begin(), indices.end(), 0);
  std::sort(indices.begin(), indices.end(),
	    [&array](int left, int right) -> bool {
	      // sort indices according to corresponding array element
	      return array[left] < array[right];
	    });

  return indices;
}


void setNDPoints(string res_file_path, int n_obj, vector< vector<double> > &point_2dvec) {
  char buffer[1024];
  
  Problem.setArtificialProblem(n_obj);
  int i = 1;

  snprintf(buffer, 1024, res_file_path.c_str());			
				
  TListSet <TPoint> allSolutions;

  allSolutions.Load(buffer);
  std::cout << "The size of the point set is " << allSolutions.size() << ".\n";  
   
  for (int l = 0; l < allSolutions.size(); l++) {
    TPoint* sol = allSolutions[l];
    int rpos = rand() % allSolutions.size();
    allSolutions[l] = allSolutions[rpos];
    allSolutions[rpos] = sol;
  }

  TTreeSet <TPoint> treeSet;

  for (int l = 0; l < allSolutions.size(); l++) {
    bool treeAdded = treeSet.Update(*((TPoint*)(allSolutions[l])));
  }

  std::cout << "The size of the nondominated point set is " << treeSet.numberOfSolutions() << ".\n";  
  
  treeSet.saveToList();
  treeSet.saveTo2dVec(point_2dvec);        
}

vector< vector<double> > distMatrix(vector< vector<double> > &archive) {  
  // Generate the distance matrix
  std::vector<std::vector<double>> dist_matrix(archive.size(), std::vector<double>(archive.size()));
  int n_obj = archive[0].size();

  std::vector<double> min_vec(n_obj, 1e+20);
  std::vector<double> max_vec(n_obj, -1e+20);
  std::vector<double> norm_vec(n_obj, 0);  

  for (int i = 0; i < n_obj; i++) {
    for (int j = 0; j < archive.size(); j++) {
      min_vec[i] = std::min(archive[j][i], min_vec[i]);
      max_vec[i] = std::max(archive[j][i], max_vec[i]);      
    }

    norm_vec[i] = (max_vec[i] - min_vec[i]);
  }
  
  for (int i = 0; i < archive.size(); i++) {
    for (int j = i; j < archive.size(); j++) {
      if (i == j) {
	dist_matrix[i][j] = 1e+20;
      }
      else {      
	double sum = 0;
	for (int k = 0; k < n_obj; k++) {
	  double tmp = (archive[i][k] - archive[j][k]) / norm_vec[k];
	  sum += tmp * tmp;
	}

	dist_matrix[i][j] = sqrt(sum);
	dist_matrix[j][i] = dist_matrix[i][j];
      }
    }    
  }

  return dist_matrix;
}

vector< vector<double> > iterativeDSS(vector< vector<double> > &archive, int pset_size) {
  std::mt19937 rnd_engine(1); 
  // Normalize the objective values of poitns in the archive.
  int n_obj = archive[0].size();
  int archive_size = archive.size();
  std::vector< std::vector<double> > psubset(pset_size, std::vector<double>(n_obj));
      
  // Randomly shuffle the order of IDs
  std::vector<int> shuffled_ids(archive.size());
  std::iota(shuffled_ids.begin(), shuffled_ids.end(), 0);
  std::shuffle(shuffled_ids.begin(), shuffled_ids.end(), rnd_engine);

  // If x[i] = 1, the i-th point in the set is selected.
  //  std::vector<int> x(archive.size(), 0);
  std::vector<int> sel_ids(pset_size+1);
  std::vector<int> tmp_ids(pset_size);
  
  for (int i = 0; i < pset_size; i++) {
    sel_ids[i] = shuffled_ids[i];
  }

  int iter = 0;
  int max_iter = 10000;
  int id = pset_size;
  int rid;
  bool isOverlapped;
  std::vector<std::vector<double>> dist_matrix = distMatrix(archive);    

  while (iter < max_iter) {   
    do {
      rid = shuffled_ids[id];
      isOverlapped = false;
      
      for (int i = 0; i < pset_size; i++) {
	if (sel_ids[i] == rid) {
	  isOverlapped = true;
	  break;
	}
      }

      id++;
      if (id >= archive.size()) {
	id = 0;
	std::shuffle(shuffled_ids.begin(), shuffled_ids.end(), rnd_engine);
      }
    } while(isOverlapped);

    sel_ids[pset_size] = rid;
    
    int worst_id = pset_size;
    double worst_ulevel = -1;
    
    for (int i = 0; i < pset_size+1; i++) {
      int pid = 0;
      for (int j = 0; j < pset_size+1; j++) {      
      	if (i != j) {
	    tmp_ids[pid] = sel_ids[j];
	    pid++;
      	}
      }
	
      double ulevel = 1e+20;
      for (int i = 0; i < pset_size; i++) {
      	for (int j = 0; j < pset_size; j++) {
      	  if (i != j && dist_matrix[tmp_ids[i]][tmp_ids[j]] < ulevel) {
      	    ulevel = dist_matrix[tmp_ids[i]][tmp_ids[j]];
      	  }
      	}
      }

      if (ulevel > worst_ulevel) {
	worst_ulevel = ulevel;
	worst_id = i;
      }                        
    }

    if (worst_id < pset_size) {
      sel_ids[worst_id] = sel_ids[pset_size];      
    }    
      
    iter++;
  }

  for (int i = 0; i < pset_size; i++) {
    for (int j = 0; j < n_obj; j++) {
      psubset[i][j] = archive[sel_ids[i]][j];
    }
  }  

  return psubset;
}

void save2DVecToFile(string save_file_path, vector< vector<double> > &point_2dvec) {
  // Save to the file
  ofstream ofs(save_file_path.c_str());
  int dim_row = point_2dvec[0].size();
  int dim_col = point_2dvec.size();
  
  for (int i = 0; i < dim_col; i++) {
    for (int j = 0; j < dim_row - 1; j++) ofs << point_2dvec[i][j] << ",";
    ofs << point_2dvec[i][dim_row - 1] << endl;
  }
  
  ofs.close();
}

double euc_distance(std::vector<double> p, std::vector<double> q, int n_obj){
  double sum = 0;
  for (int i = 0; i < n_obj; i++) {
    sum += (p[i] - q[i]) * (p[i] - q[i]);    
  }
  
  return sqrt(sum);
}

int main(int argc, char **argv) {
  srand(13);
  
  int n_obj;
  int pset_size;
  double roi_radius;
  string res_file_path;
  string pp_file_path;  
  string pref_point_file_path;
 
  int index = 1;
  std::string flag;
  std::string value; 

  while (index < argc) { // read (flag,value) pairs from command line arguments
    flag = argv[index];
    value = argv[index+1];
    //    std::cout << flag << " " << value << std::endl;

    if (flag == "-res_file_path") {
      res_file_path = value;
    }
    else if (flag == "-pp_file_path") {
      pp_file_path = value;
    }
    else if (flag == "-pref_point_file_path") {
      pref_point_file_path = value;
    }       
    else if (flag == "-n_obj") {
      n_obj = std::stoi(value);
    }
    else if (flag == "-pset_size") {
      pset_size = std::stoi(value);
    }
    else if (flag == "-roi_radius") {
      roi_radius = std::stod(value);
    }    
    else {
      std::cerr << "Error. Unkown argument " << flag << std::endl;
      exit(1);
    }
    
    index += 2;
  }

  std::vector<double> reference_point(n_obj);
  
  double *weight_vec = (double *) malloc (sizeof(double) * n_obj);
  for (int i = 0; i < n_obj; i++) weight_vec[i] = 1.0 / n_obj;
  
  /* Load the file for the reference point */
  FILE *fp = fopen(pref_point_file_path.c_str(), "r");
  if (fp == NULL) {
    printf("Error! %s cannot be opened.", pref_point_file_path.c_str());
    exit(1);
  }

  for (int i = 0; i < n_obj; i++) {
    fscanf(fp, "%lf", &(reference_point[i]));
  }

  vector< vector<double> > point_2dvec;
  // Select nondominated solutions from the unbouned external archive saved in "res_file_path".
  setNDPoints(res_file_path, n_obj, point_2dvec);  

  // Case 1: The size of the nondominated point set is smaller than the specified size "pset_size"
  if (point_2dvec.size() <= pset_size) {
    std::cout << "The size of the nondominated point set is smaller than the required size (" << point_2dvec.size() << " <= " << pset_size << ")" << std::endl;
    save2DVecToFile(pp_file_path, point_2dvec);
    return 0;
  }
  // Case 2: The size of the nondominated point set is larger than the specified size "pset_size"
  else {
    std::cout << "The size of the nondominated point set is larger than the required size (" << point_2dvec.size() << " > " << pset_size << ")" << std::endl;

    // Select points in the ROI, which is based on the closest point x^c (point_2dvec[best_id]) to the reference point.
    double min_dist = 1e+20;
    double tmp_dist;
    int best_id = 0;
    for (int i = 0; i < point_2dvec.size(); i++) {
      tmp_dist = euc_distance(point_2dvec[i], reference_point, n_obj);

      if (tmp_dist < min_dist) {
	min_dist = tmp_dist;
	best_id = i;
      }
    }

    // Calculate the distance from the closest point (point_2dvec[best_id]) to each point (point_2dvec[i]).
    std::vector<double> dist2best_vec(point_2dvec.size());
    for (int i = 0; i < point_2dvec.size(); i++) {
      dist2best_vec[i] = euc_distance(point_2dvec[i], point_2dvec[best_id], n_obj);
    }

    // Add points in the approximated ROI to psubset.
    std::vector<size_t> sorted_ids = argsort(dist2best_vec);
    vector< vector<double> > psubset;
    int last_i = 0;
    int id;
    double d;
    
    for (int i = 0; i < point_2dvec.size(); i++) {
      id = sorted_ids[i];
      d = euc_distance(point_2dvec[id], point_2dvec[best_id], n_obj);
      if (d < roi_radius) {
	psubset.push_back(point_2dvec[id]);
	last_i = i;
      }
    }
    
    if (psubset.size() == pset_size) {
      save2DVecToFile(pp_file_path, psubset);
    }
    else if (psubset.size() < pset_size) {
      // Case 2.1: This is Case 1 in Figure 2(a) in the GECCO paper. If the number of points in the approximated ROI is less than the specified size "pset_size", the closest points to the best point x^c are iteratively added to psubset until psubset.size() < pset_size.
      cout << "Case 1: The number of points in the approximated ROI is less than the specified size (" << psubset.size() << " < " << pset_size << "). The remaining points are selected from the out of the approximated ROI." << endl;
      for (int i = last_i+1; i < point_2dvec.size() && psubset.size() < pset_size; i++) {
	id = sorted_ids[i];
	psubset.push_back(point_2dvec[id]);            
      }

      save2DVecToFile(pp_file_path, psubset);
    }
    else if (psubset.size() > pset_size) {
      // Case 2.2: If the number of selected points in the ROI is larger than the specified size "pset_size", distance-based subset selection is performed.
      cout << "Case 2: The number of points in the approximated ROI is greater than the specified size (" << psubset.size() << " > " << pset_size << "). The iterative distance-based subset selection is performed." << endl;      
      vector< vector<double> > psubset_ = iterativeDSS(psubset, pset_size);
      save2DVecToFile(pp_file_path, psubset_);      
    }            
  } 
    
  return 0;
}
