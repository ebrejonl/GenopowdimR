#include <Rcpp.h>
using namespace Rcpp;

// Function to calculate distances based on mismatches in microsatellite data
// [[Rcpp::export]]
NumericMatrix calculate_distances(CharacterMatrix genetics, double NA_weight,  int min_common_loci) {
    // Step 1: Initialize variables
    int n = genetics.nrow();       // Get the number of rows (individuals) in the genetics matrix
    int m = genetics.ncol();       // Get the number of columns (alleles) in the genetics matrix
    int num_loci = m / 2;          // Calculate the number of loci pairs (each locus consists of 2 columns)
    NumericMatrix dist_matrix(n, n); // Create an n x n matrix to store the distances between individuals
    
    // Step 2: Loop through each pair of individuals
    for (int i = 0; i < n; i++) {           // Outer loop: iterate over each individual
        for (int j = i; j < n; j++) {       // Inner loop: iterate over each individual again, starting from i to avoid redundancy
            double distance = 0;            // Initialize the distance between individual i and individual j to 0
            
            // Step 3: Check if individuals i and j have at least 8 loci pairs in common without "x"
            int common_count = 0;
            for (int k = 0; k < num_loci; k++) {   // Loop through each locus pair
                // Get indices of alleles for locus pair k
                int idx1 = k * 2;         // First allele of locus k (odd index)
                int idx2 = k * 2 + 1;     // Second allele of locus k (even index)
                
                // Check alleles of individuals i and j at locus pair k
                if (genetics(i, idx1) != "x" && genetics(j, idx1) != "x" &&
                    genetics(i, idx1) == genetics(j, idx1) &&
                    genetics(i, idx2) != "x" && genetics(j, idx2) != "x" &&
                    genetics(i, idx2) == genetics(j, idx2)) {
                    common_count++;
                }
            }
            
            // Step 4: Calculate the distance based on the conditions
            if (i != j && common_count >= min_common_loci) { // Only calculate distances if i is not equal to j and have at least X loci pairs in common without "x"
                for (int k = 0; k < num_loci; k++) { // Loop through each locus pair
                    // Get indices of alleles for locus pair k
                    int idx1 = k * 2;         // First allele of locus k (odd index)
                    int idx2 = k * 2 + 1;     // Second allele of locus k (even index)
                    
                    // Check alleles of individuals i and j at locus pair k
                    if (genetics(i, idx1) == "x" || genetics(j, idx1) == "x" ||
                        genetics(i, idx2) == "x" || genetics(j, idx2) == "x") { // If either individual has an "x" at either allele
                        distance += NA_weight; // Add the NA_weight to the distance
                    } else if (genetics(i, idx1) != genetics(j, idx1) || genetics(i, idx2) != genetics(j, idx2)) { // Otherwise, if there is a mismatch
                        distance += 1;       // Add 1 to the distance
                    }
                }
            } else {
                distance = (NA_weight + 6) * m; // Assign a high value if less than 8 loci pairs in common without "x"
            }
            
            dist_matrix(i, j) = distance; // Assign the calculated distance to the distance matrix
            dist_matrix(j, i) = distance; // Since the distance matrix is symmetric, assign the same value to the symmetric element
        }
    }
    
    return dist_matrix; // Return the completed distance matrix
}