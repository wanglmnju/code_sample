/*test.cpp
Read the real space vector R1/R2 and then the HDelta Matrix (16*16 orbitals) <R1|HDelta|R2>. Do Fourier transform of HDelta i.e. Sum exp(i*K1*R1)*exp(-i*K2*R2)*<R1|HDelta|R2> for every R1-R2 pair.
Read two projector matrix K1N1J1 and K2N2J2 (16*16), and calculate K1N1J1*FFT of HDelta*K2N2J2
Input: HDelta file name
Output: K1N1J1*FFT of HDelta*K2N2J2 
*/
#include <iostream>
#include <complex>
#include <fstream>
using namespace std;

int main(int argc, char *argv[])
{
  complex<double> HDelta[773][16][16]; /*Hamiltonian Delta Matrix HDelta[k][i][j],k-the R1-R2 pair, ij - the orbitals */
  complex<double> T1[16][16],T2[16][16];//Result Matrix
  complex<double> FFTHDelta[16][16];//Fourier transform of HDelta
  complex<double> K1N1J1[16][16],K2N2J2[16][16],K1N1J1C[16][16];//Projector matrix
  
  //Input Files
  ifstream HDeltaIn("HDelta.in",ios::in);//Open HDelta file
  fstream HDeltaOut("HDelta.out",ios::out);//Open HDelta.out file for comparing
  ifstream K1N1J1File("K1N1J1.dat",ios::in);//Open projector matrix 1 file
  ifstream K2N2J2File("K2N2J2.dat",ios::in);//Open projector matrix 2 file

  //Output Files
  fstream Result("Result.dat",ios::out);//Resulting file
  fstream Logfile("log.dat",ios::out);//Log File

  double R1[3],R2[3];//Real space vector R1 and R2
  double K1[3],K2[3];//Reciprocal vector K1 and K2
  complex<double> exponent1,exponent2,exponent,cmpx,cmpx1,cmpx2;
  double exp1,exp2;

  //Initializing reciprocal vector K1 and K2
  cout <<"Input K1 vector:"<<endl;
  cin >> K1[0] >> K1[1] >> K1[2];
  cout <<"Input K2 vector:"<<endl;
  cin >> K2[0] >> K2[1] >> K2[2];

  Logfile <<"K1 "<< K1[0] << " " <<K1[1] << " " <<K1[2]<< endl;
  Logfile <<"K2 "<< K2[0] << " " <<K2[1] << " " <<K2[2]<< endl;

  //Initializing results matrix and FFT matrix of HDelta to 0
  for (int i=0;i<16;i++)
  {
    for (int j=0;j<16;j++)
	{
	  FFTHDelta[i][j]=complex<double>(0,0);
	}
  }

/*Fourier transoform of <R1|HDelta|R2> */
  for (int k=0;k<773;k++)
  { //Read HDelta loops and every step for per R1-R2 pair
    for (int i=0;i<3;i++)//Read Real space vector R1 from HDelta file
    {
	  HDeltaIn >> R1[i]; 
	  HDeltaOut << R1[i] << " ";
    }
    HDeltaOut <<"  ";
    for (int i=0;i<3;i++)//Read R2
    {
		HDeltaIn >> R2[i];
		HDeltaOut << R2[i]<< " ";
    }
    HDeltaOut<<endl;

    //Calculate exp(i*K1*R1)*exp(-i*K2*R2)
    exp1= 2*3.1415926*(K1[0]*R1[0]+K1[1]*R1[1]+K1[2]*R1[2]);
    cmpx1=complex<double>(0,exp1);
    HDeltaOut << exp1 << cmpx1;
    exponent1=exp(cmpx1);
    exp2=-2*3.1415926*(K2[0]*R2[0]+K2[1]*R2[1]+K2[2]*R2[2]);
    cmpx2=complex<double>(0,exp2);
    HDeltaOut << exp2 << cmpx2;
    exponent2=exp(cmpx2);
    exponent=exponent1*exponent2;
    HDeltaOut << exponent1 << exponent2 << exponent << endl;

    for(int i=0;i<16;i++)
    {
      for(int j=0;j<16;j++)
      {
        HDeltaIn >> HDelta[k][i][j]; //Read the matrix element of <R1|HDelta|R2>
        cmpx = exponent * HDelta[k][i][j];
        HDeltaOut << cmpx;
        // exponent*HDelta Matrix Element and then sum
	FFTHDelta[i][j] = FFTHDelta[i][j]+cmpx;
      }
      HDeltaOut << endl;
    } 
  }//End of the Fourier transform and store it in FFTHDelta[i][j]

  //Output FFTHDelta to log file
  Logfile << "FFT HDelta:"<<endl;
  for (int i=0;i<16;i++)
  {
    for (int j=0;j<16;j++)
      Logfile << FFTHDelta[i][j];
    Logfile << endl;
  }
  Logfile <<endl;

  //Read the projector matrix and write them to log file
  for (int i=0;i<16;i++)
  {
    for (int j=0;j<16;j++)
    {
	K1N1J1File >> K1N1J1[i][j];
	K2N2J2File >> K2N2J2[i][j];
    }
  }
  Logfile << "K1N1J1:" << endl;
  for (int i=0;i<16;i++)
  {
	for (int j=0;j<16;j++)
	{
		Logfile << K1N1J1[i][j];
	}
	Logfile << endl;
  }
  Logfile << "Transposed conjugate of K1N1J1:"<<endl;
  for (int i=0;i<16;i++)
  {
    for (int j=0;j<16;j++)
    {
	K1N1J1C[i][j]=conj(K1N1J1[j][i]);  
	Logfile << K1N1J1C[i][j];
    }
	Logfile << endl;
  }
  Logfile <<endl<<"K2N2J2:"<<endl;
  for (int i=0;i<16;i++)
  {
    for (int j=0;j<16;j++)
    {
	  Logfile << K2N2J2[i][j];
    }
	Logfile << endl;
  }

  //Perform K1N1J1*FFTHDelta*K2N2J2
  //T1=FFTHDelta*K2N2J2 Matrix multiply
  Logfile << "T1=FFTHDelta*K2N2J2"<<endl;
  for (int i=0;i<16;i++)
  {
	  for (int j=0;j<16;j++)
	  {
		  T1[i][j]=0;
		  for (int k=0;k<16;k++)
		  {
			  T1[i][j]=T1[i][j]+FFTHDelta[i][k]*K2N2J2[k][j];
		  }
		  Logfile << T1[i][j];
	  }
	  Logfile << endl;
  }

  //T2=K1N1J1*T1 Matrix multiply
  Logfile << "T2=K1N1J1*T1"<<endl;
  for (int i=0;i<16;i++)
  {
	  for (int j=0;j<16;j++)
	  {
		  T2[i][j]=0;
		  for (int k=0;k<16;k++)
		  {
			  T2[i][j]=T2[i][j]+K1N1J1C[i][k]*T1[k][j];
		  }
		  Logfile << T2[i][j];
		  Result << T2[i][j];
	  }
	  Logfile << endl;
	  Result << endl;
  }
  
  HDeltaIn.close();
  HDeltaOut.close();
  Logfile.close();
  K1N1J1File.close();
  K2N2J2File.close();
  Result.close();

  return 0;
}
