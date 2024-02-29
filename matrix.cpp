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
#include <TGraph.h>
#include <TF1.h>
using namespace std;

double osc_func(double x, double *param) {
	double value = 1.0- param[0]*pow(sin(1.267*param[1]*470.0/x), 2.0);
	
	return value;
}


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
	int size = matrix.size();
	int dim = matrix[0].size();
	TMatrixD cov_matrix(dim, dim);
	double sum_1 = 0.0;
	double sum_2 = 0.0;
	double prod = 0.0;
	cout << "DEBUG: gen_mat_cov" <<endl;
	cout << "dim=" <<dim <<endl;
	cout << "size=" << size << endl;
	
	for(int i=0; i<dim; i++) {
		for(int j=i; j<dim; j++) {
			for(int k=0; k<size; k++) {
				prod += matrix[k][i]*matrix[k][j];
				sum_1 += matrix[k][i];
				sum_2 += matrix[k][j];
			}
			cov_matrix[i][j] = prod/size - (sum_1/size)*(sum_2/size);
			cov_matrix[j][i] = cov_matrix[i][j];
			sum_1=0;
			sum_2=0;
			prod = 0;
			
						
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
	double ymax = 50;
	TCanvas* canvas = new TCanvas("canvas", "TMatrix Plot", 800, 600);
	for (int i=0; i<psuedo_vector.size(); i++) {
		TH1F *histogram = new TH1F(Form("hist%d", i), Form("Histogram %d", i), 25, 0, 25);
		for(int j= 0; j<25; j++) {
			histogram->Fill(j+1, psuedo_vector[i][j]);
			if(i==0) {
				histogram->Sumw2(0);
				histogram->SetMaximum(ymax);
				histogram->SetLineColor(kRed);
				histogram->Draw();
			}
			else {
				
				histogram->Sumw2(0);
				histogram->SetMarkerSize(5.5);
				histogram->Draw("same");
			}

		}
	}
	TH1F *histogram = new TH1F("hist", "Pseudoexperiments", 25, 0, 25);
	for(int i=0; i<25; i++) {
		histogram->Fill(i+1, psuedo_vector[0][i]);
	}
	histogram->Sumw2(0);
	histogram->SetMaximum(ymax);
	histogram->SetLineColor(kRed);
	histogram->Draw("same");


	canvas->Draw();
}

void matrix_plot(TMatrixT<double> &matrix) 
{
	int n_rows = 25;
	int n_cols = 25;
	TH2F *hist = new TH2F("mat_hist", "Matrix Plot", n_cols, 0, n_cols, n_rows, 0, n_rows);
	for(int i=0; i<n_rows; i++) {
		for( int j=0; j<n_rows; j++) {
			hist->SetBinContent(i+1, j+1, matrix(i,j));
		}
	}

	// hist->SetMinimum(-0.3); // Replace 'minimum_value' with the desired minimum value
	// hist->SetMaximum(1); // Replace 'maximum_value' with the desired maximum valu

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

void percent_diff(TMatrixT<double> &val, TMatrixT<double> &ref) 
{
	int dim = 25;
	TH2F *hist = new TH2F("h1", "Percent Difference", dim,0,dim,dim,0,dim);
	for(int i=0; i<dim; i++) {
		for( int j=0; j<dim; j++) {
			hist->SetBinContent(i+1, j+1, 100*((val(i,j)-ref(i,j))/ref(i,j)));
		}
		hist->Draw("colz");

	        hist->SetMinimum(-20); // Replace 'minimum_value' with the desired minimum value
	hist->SetMaximum(20); // Replace 'maximum_value' with the desired maximum valu
	}
}

	
int main(int argc, char **argv) 
{
	TApplication app("app", &argc, argv);
	// Use epoch time to seed since I am a linux nerd
	auto now = std::chrono::system_clock::now().time_since_epoch();
	unsigned int seed = static_cast<unsigned int>(now.count());
	TRandom3 randomgen;
	randomgen.SetSeed(seed);
	// Variables for holding extracted data
	string filename = "file_total_80dm2_2tue_t24_input.root";
		
	TFile root_file(filename.c_str());
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
	// cout << "Matrix U" <<endl;
	// matrix_stat(U);
	// cout << "Matrix S"<<endl;
	// matrix_stat(S);

	// matrix_plot(cov_transformed);
	
	TMatrixD spectrum_transformed = U_t*S; 
	// Pull out the sigma^2 values and create a lambda array
	// Loop over to create S'1-S'N values
	vector<double> product;
	vector<double> product_sum;
	// Here is where we will only care about 25 values
	// A Vector of vectors of double values
	vector<vector<double>> spectrum;
	vector<double> temp;			
	int num_psuedo = 10000;
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
	
	// Now I need to generate the the covariance matrix
	// TMatrixD mat_cov_psuedo = gen_mat_cov(spectrum);
	// Generate a corr matrix to test against the original mat_rho
	// TMatrixD mat_cor_psuedo(25,25);
	// for(int i=0; i<25; i++) {
	// 	for(int j=0; j<25; j++) {
	// 		mat_cor_psuedo(i, j) = mat_cov_psuedo(i, j) / (TMath::Sqrt(mat_cov_psuedo(i,i)*mat_cov_psuedo(j,j)));
			
			
	// 	}
	// }
	// psuedo_plt(spectrum);
	// Validate the Covariance Matrix
	// TH2D *cov_hist = new TH2D("hist", "Covariance Matrix", 25, 0, 25, 25, 0, 25);
	
	// for(int i=0; i<25; i++) {
	// 	for(int j=0; j<25; j++) {
	// 		cout << mat_cov_psuedo(i,j);
			
	// 	}
	// }
	// cov_hist->Draw("colz");
	// matrix_plot(mat_cov_psuedo);
	
	// matrix_plot(mat_cor_psuedo);
	// percent_diff(mat_cov_psuedo, mat_cov);
	// percent_diff(mat_cor_psuedo, mat_rho);
	
	
	
	// cout << "Cov Stats" <<endl;
	// matrix_stat(cov_transformed);
	// matrix_stat(mat_cov_psuedo);
	
	// canvas->Update();
	cout << "Inverting Cov";
	
// Invert the cov matrix
	
	TMatrixD cov = mat_cov.GetSub(0, 24, 0, 24);

	TMatrixD cov_invert = cov.Invert();
	cout << " Done." <<endl;
	
	cout << "Generating predicted array";
	
// Generate the Predicted Array
	double param[2]={0, 0};
	
	TMatrixD pred(1,25);
	for(int i=0; i<25; i++) {
		pred[0][i] = osc_func(i*10, param)*mat_cv[0][i];
		cout << pred[0][i]<<endl;
		
	}
	pred[0][0]=0;
	
	cout <<" Done." <<endl;
	
	double chisquare_array[num_psuedo];
	cout << "Populating chisquare arrat";
	cout <<"Cov Inverted" <<endl;
	matrix_stat(cov_invert);
	

	for(int i =0; i<num_psuedo; i++) {
		TMatrixD diff_mp(1,25);
		// Calculate M-P
		for(int j=0; j<25; j++) {
			diff_mp[0][j] = spectrum[i][j] - pred[0][j];
		}
		TMatrixD diff_mp_transposed (TMatrixD::kTransposed, diff_mp);
		// Calculate chisquare
		// I know this collapses into a single value because science.
		// Or linear algebra
		// Reasons
		chisquare_array[i] = (diff_mp*cov_invert*diff_mp_transposed)[0][0];
	}
	TCanvas *canvas = new TCanvas("canvas", "Comparison", 800, 600);
	TH1F *chi_hist = new TH1F("chi_hist", "Chi Square", 100, 0, 60);
	for(int i=0; i<num_psuedo; i++) {
		chi_hist->Fill(chisquare_array[i]);
		
		}
	chi_hist->Sumw2(0);
	// chi_hist->Draw();
	    // TF1 *fChi2 = new TF1("fChi2", "[0]*x^[1]/TMath::Gamma([1]/2)/pow(2, [1]/2)*exp(-x/2)", 30, 100);
	    // fChi2->SetParameters(1.0, 25.0);

	TF1 *chiSquarePDF = new TF1("chiSquarePDF", "ROOT::Math::chisquared_pdf(x, [0], 0)", 0, 60);
	chiSquarePDF->SetParameter(0, 25);
	chiSquarePDF->SetLineColor(kRed);
	
	chi_hist->Scale(1.0 / chi_hist->Integral()/chi_hist->GetBinWidth(1));
	chi_hist->Draw();
	// chiSquarePDF->Draw();	
	chiSquarePDF->Draw("same");		     
	app.Run();
	
	return 0;
}

