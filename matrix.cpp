#include <TH2F.h>
#include <TApplication.h>
#include <TCanvas.h>
#include <TRootCanvas.h>
#include <TMath.h>
#include <TFile.h>
#include <TMinuit.h>
#include <TObject.h>
#include <TRandom.h>
#include <cstddef>
#include <iostream>
#include <TMatrixD.h>
#include <TMatrixT.h>
#include <TVectorD.h>
#include <TTree.h>
#include <iostream>
#include <TDecompSVD.h>
#include <TH1F.h>
#include <ostream>
#include <TRandom3.h>
#include <cmath>
using namespace std;

TRandom3 randomgen;

double calc_avg(vector<double> &data) {
	double sum = 0.0;
	for(double &value : data) {
		sum += value;
	}
	return sum/data.size();
}

double calc_covariance(vector<double> &data1, vector<double> &data2) {
	double avg1 = calc_avg(data1);
	double avg2 = calc_avg(data2);
	// cout << data1.size() <<endl;
	double sum = 0.0;
	for(size_t i=0; i<data1.size(); i++) {
		sum += data1[i]*data2[i];
	}
	double avg_product = sum/data1.size();
	return avg_product - (avg1*avg2);
}

TMatrixD gen_mat_cov(vector<vector<double>> &matrix) {
	size_t size = matrix.size();
	size_t dim = matrix[0].size();
	TMatrixD cov_matrix(dim, dim);
	for(size_t i=0; i<size; i++) {
		for(size_t j=0; j<dim; j++) {
			for(size_t k=1; k<j; k++) {
				cov_matrix[j][k] = calc_covariance(matrix[j], matrix[k]);
				cov_matrix[k][j] = cov_matrix[j][k];

			}
		}
	}
	return cov_matrix;
}

void vector_plot(vector<double> &vec) 
{
	int size = vec.size();
	TH1F *hist = new TH1F("mat_hist", "Matrix Plot", size, 0, size);
	
	for(int i=0; i<size; i++) {
		hist->Fill(i+1);
	}
	TCanvas* canvas = new TCanvas("canvas", "TMatrix Plot", 800, 600);
	hist->Draw("COLZ"); // COLZ option for a color plot
	hist->GetXaxis()->SetTitle("X-Axis");
	hist->GetYaxis()->SetTitle("Y-Axis");
	canvas->Draw();
}
 


void psuedo_plt(vector<vector<double>> &psuedo_vector) 
{
	double ymax = 26;
	int size = psuedo_vector.size();
	
	for(int i=0; i<25; i++) {
		TH1F *histogram = new TH1F(Form("hist%d", i), Form("Histogram %d", i), 25, 0, 25);
		for(int j=0; j<size; j++) {
			histogram->Fill(i+1,psuedo_vector[j].at(i));
			
		}
		if(i==0) {
			histogram->SetMarkerColor(kRed);
			histogram->SetMarkerSize(5.5);
			
			histogram->SetLineColor(kRed);
			histogram->SetMaximum(ymax);
			histogram->Draw();
			
		}
		else {
			histogram->SetMarkerSize(5.5);
			histogram->SetMaximum(ymax);
			histogram->Draw("same");
			
		}
	}
} 

void matrix_plot(TMatrixT<double> &matrix) 
{
	int n_rows = matrix.GetNrows();
	int n_cols = matrix.GetNcols();
	TH2F *hist = new TH2F("mat_hist", "Matrix Plot", n_cols, 0, n_cols, n_rows, 0, n_rows);
	for(int i=0; i<n_rows; i++) {
		for( int j=0; j<n_rows; j++) {
			hist->Fill(i+1, j+1, matrix(i,j));
		}
	}
	hist->SetMinimum(-1); // Replace 'minimum_value' with the desired minimum value
	hist->SetMaximum(1); // Replace 'maximum_value' with the desired maximum value
    // Save the histogram to a root file
    TFile outputFile("matrixHistogram.root", "RECREATE");
	TCanvas* canvas = new TCanvas("canvas", "TMatrix Plot", 800, 600);
	hist->Draw("COLZ"); // COLZ option for a color plot
	// Optional: Set axis labels
	hist->GetXaxis()->SetTitle("X-Axis");
	hist->GetYaxis()->SetTitle("Y-Axis");
	// Show the canvas
	canvas->Draw();
}

void matrix_stat(TMatrixT<double> &matrix) 
{
	cout << "Matrix Info" << endl;
	cout << "Rows: " << matrix.GetNrows() << endl;
	cout << "Cols: " << matrix.GetNcols() << endl;
}

int main(int argc, char **argv) 
{
	TApplication app("app", &argc, argv);

	// Variables for holding extracted data
	
	TFile root_file("file_total_80dm2_2tue_t24_input.root");
	TMatrixD* mat_frac = (TMatrixD*)root_file.Get("matrix_fractional_flux_Xs_G4_det");

	int nRows = mat_frac->GetNrows();
	int nCols = mat_frac->GetNcols();
	TMatrixD mat_cov(nRows,nCols);
	TMatrixD mat_rho(nRows,nCols);
	TMatrixD mat_cv(1,1092);
	TMatrixD* mat_transition = (TMatrixD*)root_file.Get("matrix_transition");
	// Generate the correlation matrix
	for(int i=0; i<nRows; i++) {
		for(int j=0; j<nCols; j++) {
			mat_rho(i, j) = (*mat_frac)(i, j) / (TMath::Sqrt((*mat_frac)(i,i)*(*mat_frac)(j,j)));
			
			
		}
	}
	// We will only want to generate COV and fluctuate values relating to first 25 bins
	// Cut this off to pull out only 25 values and change the shape
	// Generate the Cov matrix
	TTree *tree(dynamic_cast<TTree*>(root_file.Get("tree_spectrum")));
	vector<float>* cv_vec = nullptr;
	tree->SetBranchAddress("vec_energy_spectrum", &cv_vec);
	tree->GetEntry(0);
	for(int i=0;i<25; i++) {
		mat_cv[0][i] = cv_vec->at(i);
	}
	for(int i=0; i<25; i++) {
		for(int j=0; j<25; j++) {
			mat_cov(i,j) = (*mat_frac)(i,j)*mat_cv(0,i)*mat_cv(0,j);
			
			}
	}
	TMatrixD mat_transpose(TMatrixT<double>::kTransposed, *mat_transition);
	TMatrixD mat_collapse = mat_transpose*mat_cov*(*mat_transition);

	// This code is for plotting the original 
	// TH1F *hist_cv = new TH1F("histo_cv", "CV Matrix", 25, 0, 25);
	// for(int i=0; i<25; i++) {
	// 	hist_cv->Fill(i+1,mat_cv(0,i));
		
	// }
	// hist_cv->Sumw2(0);
	// hist_cv->Draw();
	
// // Validating mat_cov
// 	TMatrixD mat_rho_from_cov(nRows, nCols);
	
// 	for(int i=0; i<nRows; i++) {
// 		for(int j=0; j<nCols; j++) {
// 			mat_rho_from_cov(i, j) = mat_cov(i, j) / (TMath::Sqrt(mat_cov(i,i)*mat_cov(j,j)));
			
// 		}
// 	}

	// Validating mat_cov
	// TMatrixD mat_rho_from_cov(364, 364);
	
	// for(int i=0; i<364; i++) {
	// 	for(int j=0; j<364; j++) {
	// 		mat_rho_from_cov(i, j) = mat_collapse(i, j) / (TMath::Sqrt(mat_collapse(i,i)*mat_collapse(j,j)));
			
	// 	}
	// }
	// TH1F *frac_values = new TH1F("histo_frc", "Fractional unc", 25, 0, 25);
	// for(int i=0; i<nCols; i++) {
	// 	frac_values->Fill(i+1,(*mat_frac)(i,i));
		
	// 	// cout << TMath::Sqrt((*mat_frac)(i,i)) <<endl;
	// }
	// frac_values->Draw();
// SVD Decomposition
	TDecompSVD udv(mat_collapse);
	TMatrixD U = udv.GetU();
	TMatrixD U_t (TMatrixD::kTransposed, U);
	TMatrixD cov_transformed = U_t*mat_collapse*U;
	cout << "Matrix U" <<endl;
	TMatrixD S(364,1);
	for(int i=0; i<1; i++) {
		for(int j=0; j<25; j++) {
			S(j,i) = mat_cv(i,j);
		}
	}
	cout << "Matrix U" <<endl;
	matrix_stat(U);
	cout << "Matrix S"<<endl;
	matrix_stat(S);
	TMatrixD spectrum_transformed = U_t*S; 
	// Pull out the sigma^2 values and create a lambda array
	// Loop over to create S'1-S'N values
	vector<double> product;
	vector<double> product_sum;
	// Here is where we will only care about 25 values
	// A Vector of vectors of double values
	vector<vector<double>> spectrum;
	vector<double> temp;			
	int num_psuedo = 100;
	// Save the original spectrum in slot 0
	for(int i=0; i<364; i++) {
		temp.push_back(spectrum_transformed(i,0));
	}
	spectrum.push_back(temp);
	// Generate a series of psuedo experiments and save each spectrum in the spectrum slots
	for(int i=0; i<num_psuedo; i++) {
		vector<double> row;
		for(int j=0; j<364; j++) {
			row.push_back(randomgen.Gaus(spectrum_transformed(j,0), sqrt(cov_transformed(j,j))));
			// element_sum_array[j] += psuedo_array[j];
		}
		spectrum.push_back(row);
	}

	// Now I need to use U to get back into realspace
	for(int i=0; i<num_psuedo; i++) {
		// Copy the vector into a TVector
		TVectorD temp_tvec(364, spectrum[i].data());
		// Perform the multiplication
		TVectorD result_tvec = U*temp_tvec;
		
		// Convert result TVec back into a vector
		vector<double> result_vec(result_tvec.GetNrows());
		for(int i=0; i<result_tvec.GetNrows(); i++) {
			result_vec[i] = result_tvec[i];
		}
		// Now Tack the vectir back into spectrum
		spectrum[i] = result_vec;
		
		
	}

	// Calculate the covariance Matrix
	for(int i=0; i<25; i++) {
		for(int j=0; j<i; j++)
			for(int k=0; k<num_psuedo; k++) {
				spectrum[k][j]*spectrum[k][j];
				
			}
	}
	// This code is for plotting the original 

	TH1F *hist_cv =  new TH1F("histo_cv", "CV Matrix", 25, 0, 25);
	for(int i=0; i<25; i++) {
		hist_cv->Fill(i+1,spectrum[99].at(i));
		
	}
	hist_cv->Sumw2(0);
	hist_cv->Draw();

	// psuedo_plt(spectrum);

	
	
	// Now I need to generate the the covariance matrix
	// TMatrixD mat_cov_psuedo = gen_mat_cov(spectrum);
	// cout << "Cov Stats" <<endl;
	// matrix_stat(cov_transformed);
	// matrix_stat(mat_cov_psuedo);
	
	// canvas->Update();
	app.Run();		
	return 0;
}

