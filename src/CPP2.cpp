#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix calculate_distances(CharacterMatrix genetics, double NA_weight,  int min_common_loci) {
    int n = genetics.nrow();       
    int m = genetics.ncol();       
    int num_loci = m / 2;           
    NumericMatrix dist_matrix(n, n); 
    
    for (int i = 0; i < n; i++) {           
        for (int j = i; j < n; j++) {        
            double distance = 0;            
            
            int common_count = 0;
            for (int k = 0; k < num_loci; k++) {   
                int idx1 = k * 2;         
                int idx2 = k * 2 + 1;    
                
                if (genetics(i, idx1) != "x" && genetics(j, idx1) != "x" &&
                    genetics(i, idx1) == genetics(j, idx1) &&
                    genetics(i, idx2) != "x" && genetics(j, idx2) != "x" &&
                    genetics(i, idx2) == genetics(j, idx2)) {
                    common_count++;
                }
            }
            
            if (i != j && common_count >= min_common_loci) {
                for (int k = 0; k < num_loci; k++) { 
                    int idx1 = k * 2;         
                    int idx2 = k * 2 + 1;     
                    
                    if (genetics(i, idx1) == "x" || genetics(j, idx1) == "x" ||
                        genetics(i, idx2) == "x" || genetics(j, idx2) == "x") { 
                        distance += NA_weight; 
                    } else if (genetics(i, idx1) != genetics(j, idx1) || genetics(i, idx2) != genetics(j, idx2)) { 
                        distance += 1;      
                    }
                }
            } else {
                distance = (NA_weight + 6) * m; 
            }
            
            dist_matrix(i, j) = distance; 
            dist_matrix(j, i) = distance; 
        }
    }
    
    return dist_matrix; 
}