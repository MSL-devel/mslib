
#include <RInside.h>                    // for the embedded R via RInside

using namespace std;

int main(int argc, char *argv[]) {

    // create an embedded R instance 
    RInside R;

    /* 
       rinside_sample0 from RInside package
	 Shows
	     how you can call R commands using a parseEvalQ(string);
     */
    cout << "\n\n ************   TEST 1 : Hello World ***************** \n\n\n";
    R["txt"] = "Hello, world!\n";	// assign a char* (string) to 'txt'
    R.parseEvalQ("cat(txt)");           // eval the init string, ignoring any returns




    /* 
       rinside_sample1 from RInside package
        Shows 
              how you can get data into R and 
              how you can get data back from R
    */
    cout << "\n\n ************   TEST 2 : Passing Data ***************** \n\n\n";
    SEXP ans;                                 // SEXP data type for getting data back from R

    // Create a NumericMatrix  (4 by 4)
    Rcpp::NumericMatrix M(4,4);
    for (int i=0; i<4; i++) {
        for (int j=0; j<4; j++) {
            M(i,j) = i*10+j; 
        }
    }

    R["M"] = M;                                 // assign C++ matrix M to R's 'M' var

    std::string evalstr = "\
        cat('Running ls()\n'); print(ls());                    \
        cat('Showing M\n'); print(M);                          \
        cat('Showing colSums()\n'); Z <- colSums(M); print(Z); \
        Z";                     // returns Z

    ans = R.parseEval(evalstr);                 // eval the init string -- Z is now in ans
                                                
    Rcpp::NumericVector v(ans);                 // convert SEXP ans to a vector of doubles
    for (int i=0; i< v.size(); i++) {           // show the result
      std::cout << "In C++ element " << i << " is " << v[i] << " corresponds to sum of column "<<i+1<<" of matrix M"<<std::endl;
    }


    cout << "\n\n ************* Done *************"<<endl;

    exit(0);
}

