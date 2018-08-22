#include<iostream>
#include<cmath>
#include<complex>
#include<fstream>       
#include "/armadillo-4.100.2/include/armadillo"
#include<time.h>
#include"FormFac.h"
 


using namespace std;
using namespace arma;




const double PI = 4.0*atan(1.0);



int main()
{

  clock_t time1, time2;

  time1 = clock();



  double a = 1.0; //  0.5* 6.261; //1.0; // this corresponds to "a" in the code and "c" by standard notation; it is along the chain.
  double aa =1.0;  //  7.445; //1.0; // this should be a; it is perpendicular to chain. 

  // Reciprocal lattice units. Refer to p. 2139 

  double ParU =  1.0; // (2.0*PI)/(2.0*a);  //1.0;  // (2.0*PI)/c;
  double PerpU = 1.0; //  (4.0*PI)/(aa * sqrt(3)); //1.0;

  double Jz = 13.8;  // 2.0*6.01;  // 13.8;  // 14.4;  // 13.8;  // 12.2711; // 13.8;   // 13.8; // 14.4 is to reconcile theoretical results and experimental data                                          // meV.
  double Jt = 1.9;  // 1.9;  // 1.9;       //  0.13*2.0*6.01;  // 1.9;// 2.08608;   //1.9; //  11.0*0.2    // 1.8355;   
  double B = -12.0; //  3.0; //0.5;  // 3.0                                            // T 
  double DM = 0.0; // Jt/7.0;// (1.9)/(7.0);  //  Jt/7.0                                        // meV
  double gfac =  2.00231930436153;  
  double mu =  5.7883818066e-5 * 1000;                       // meV T^-1

  // Next neareset neighbour interaction

  double Jn = 0.0;  // 2.0*0.44;// 0.9;
  double vareps = 0.0;  // 0.65; //  0.24;

  double Bz = 0.0; // -7.0;


  double T = 15.75;   // central value about which Gaussian convolution takes place.     
  double P = 0.0;  // Polarisation 
  double eta = 0.0001;  // 0.075;  // 0.1;  // 0.05;  // 0.0001;  // 0.3; //0.0001; //0.5;  // 0.001;  // PI/10.0; // 0.0001;  // PI/10.0;  //0.0001;  // PI/10.0; //0.0001;   // 0.3;   // 0.0001;   // PI/10.0;     // 0.0001;          //  PI/10.0;   // 0.1;  // 0.0001;                                          // This is just a



complex<double>ii = complex<double>(0.0,1.0);               // Defining the purely imaginary unit vector.                                                                                               


cx_mat H00;
cx_mat HNN; // Surface site Hamiltonian containing next nearest neighbour interaction

cx_mat d;
cx_mat dDag;



complex<double> g;
complex<double> gDag;

complex<double> t;
complex<double> tDag;


 // DM Normal

 
   complex<double> L;
  complex<double> LDag;






 mat id = eye<mat>(4,4);





 int  M = 20;  // 20;



 cx_mat alphaM(4*M,4);
 alphaM.zeros(); 


 
 cx_mat betaM(4*M,4);
 betaM.zeros();


 cx_mat epsilonSM(4*M,4);
 epsilonSM.zeros();



 cx_mat epsilonM(4*M,4);
 epsilonM.zeros();




 // Energy



 int J = 300; 
 
 double EnM[J];




// Below, the subatrix method is exploited .submatrix(first row, first column, last row, last column). Using this in conjunction with the large matrices defined above, one can perform the recursive algortihm simimlar to that in Mathemaitca. It will be something like the FIRST row will be in the set {0, 4, 8, 12, 16..} and the LAST row will be in the set {3, 7, 11, 15, ...}; the first column and last column will stay the same, 0 and 1 respectively. However, because the indexing in the recursive relation start at 1, the indices in submat must be as follows: submat(4*i,0, (4*i)+3 ,3).In the recursion relation, we wish to relate one submatrix with the submatrix just before it. If we have submat(12, 0, 15, 3) we want to relate to submat(8, 0, 11 , 3). In general, we are relatingsubmat (4*i, 0 , 4*i + 3, 3) to submat( 4*(i - 1), 0 , 4*(i -1), 3  )  


 cx_mat fineSM(4*J,4); // This matrix is used to store all the final value of epsilonSM for a given energy value. J is included instead of M since we are dealing with given energies.

 cx_mat G00tab(4*J,4); // This matrix is stores the value of the Green function for a given energy value.


 // Now including matrices for different Q values


 // Q_{||} or Q_z

 int H = 300;

 double QQ[H];



 int F = 1;  

  double QP[F];






 
 mat SxxNum3D (H,J);
  
 mat SyyNum3D(H,J);

 mat SzzNum3D(H,J); 

 mat SyzzyNum3D (H,J);              // Antiymmetric 

 mat SyzzyNum3DS (H,J);             // Symmetric   

 mat SxyyxNum3D(H,J);

 mat SzxxzNum3D(H,J);               // Symmetric

 mat SzxxzNum3DAS(H,J);             // Antisymmetric

 mat SzxxzNum3DAS2(H,J);            // Antisymmetric

  mat SdotSPerp(H,F );            // Unpolarised part of cross section for Q_||, Q_perp plots

  mat SyzzyPerpNum3D(H,F);        // Polarised part of cross section for Q_||, Q_perp plots

 mat TOTAL(H,J);                    // For total cross section, specifically for an energy scan with Q_perp = 2.2 and Q_|| = 1.4. 

 mat SxyyxNum3D02(H,J);      // Antisymmetic in 0_2 sector. 

 mat SyzzyNum3D02(H,J);     // Antisymmetric in 0_2 sector. 

 mat SzxxzNum3D02(H,J);     // Symmetric in 0_2 sector.

  mat SyzzyPerpNum3DPlus(H,F);    // Second term in total cross section with P = 1.0.

  mat SyzzyPerpNum3DMinus(H,F);    // Second term in total cross section with P = -1.0.

 
  mat TOTALQ_PARQ_PERP(H,F);     // Total cross section for Q_||, Q_perp for fixed energy.  

 mat TOTALQ_PARQ_PERP_MINUS_PLUS(H,F);     // Total cross section for Q_||, Q_perp for fixed energy.  

 double  NatPhysScanPlus[J];   // I+  For comparison with Nat Phys Fig 3  22.02.2015

 double  NatPhysScanMinus[J];   //I-  For comparison with Nat Phys Fig 3  22.02.2015


 int N = H*J;

  mat LM(N,F);   // This matrix to  store all the values of Total Cross Section for different energies. 



 
 for( int u = 0; u < H; u++ )
{



          QQ[u] = (u*PI*2.5)/(H-1);
 
  
  g = (Jt/2.0)*(1.0 + cos(2.0 * QQ[u] * a ) - ii * sin( 2.0 * QQ[u] * a ) );

  gDag = (Jt/2.0)*(1.0 + cos(2.0 * QQ[u] * a ) + ii * sin( 2.0 * QQ[u] * a ));

  t = ((gfac * mu * B)/2.0) * (1.0 + cos( QQ[u] * a ) - ii * sin( QQ[u] * a ));

  tDag = ((gfac * mu * B)/2.0) * (1.0 + cos( QQ[u] * a ) + ii * sin( QQ[u] * a ));

  

  // DM normal

   L =  ((ii*  DM)/2.0) *( 1.0 - cos( QQ[u] * a) + ii * sin( QQ[u] * a) );

   LDag = ((-ii*DM)/2.0) *(1.0 - cos( QQ[u] * a) - ii* sin( QQ[u] * a )  ); 



 
     
 H00 << Jz + gfac*mu*Bz + 2.0*Jn << t+L << g << 0 << endr
     << tDag+ LDag << Jz + 2.0*Jn << t + L << g << endr
     << gDag << tDag + LDag << Jz + gfac*mu*Bz +2.0*Jn  << t + L << endr
     << 0 << gDag << tDag + LDag << Jz  + 2.0*Jn << endr;



 // Matrix representation of Hamiltonian for surface site



 HNN << Jz + gfac*mu*Bz + 1.0*Jn - vareps*Jn*cos(2.0*QQ[u]) << t+L << g << 0 << endr
     << tDag+ LDag << Jz + 2.0*Jn - vareps*Jn*cos(QQ[u])  << t + L << g << endr
     << gDag << tDag + LDag << Jz + gfac*mu*Bz +2.0*Jn  << t + L << endr
     << 0 << gDag << tDag + LDag << Jz  + 2.0*Jn << endr;




  d << 0 << 0 << 0 << 0 << endr
    << 0 << 0 << 0 << 0 << endr
    << g << 0 << 0 << 0 << endr
    << t + L  << g << 0  << 0 << endr; 



 dDag << 0 << 0 << gDag  << tDag + LDag << endr
      << 0 << 0 << 0 << gDag << endr
      << 0 << 0 << 0 << 0 << endr
      << 0 << 0 << 0 << 0 << endr;  




 for(int n = 0; n < J; n++)
   {
 

     

       EnM[n] = 7.0 + 13.5*((double)n/ J); 


 // First submatrices are initialised.   

 alphaM.submat(0,0,3,3) = d * inv(EnM[n]*id - H00 + ii * eta * id  ) * d;

 betaM.submat(0,0,3,3) = dDag * inv(EnM[n]*id - H00 + ii * eta * id  ) * dDag;

 epsilonSM.submat(0,0,3,3) =  HNN + d * inv(EnM[n]*id - H00 + ii * eta * id  ) * dDag;  //HNN

 epsilonM.submat(0,0,3,3) = H00 + d * inv(EnM[n]*id - H00 + ii * eta * id  ) * dDag + dDag * inv(EnM[n]*id - H00 + ii * eta * id  ) * d; 
  


 


 for(int i = 1; i < M; i++)
   {

    
    


     alphaM.submat(4*i, 0 , (4*i) + 3, 3 ) = alphaM.submat(4*(i - 1), 0 ,(4*(i - 1) + 3) , 3 ) * inv(EnM[n] * id - epsilonM.submat(4*(i - 1), 0 ,(4*(i - 1) + 3) , 3) +  ii * eta * id    ) * alphaM.submat(4*(i - 1), 0 ,(4*(i - 1) + 3) , 3);




     betaM.submat(4*i, 0 ,(4*i) + 3 ,3) = betaM.submat(4*(i - 1), 0, (4*(i - 1) + 3  ), 3 ) * inv(EnM[n] * id  - epsilonM.submat(4*(i - 1), 0, (4*(i - 1) + 3), 3 ) + ii * eta * id ) * betaM.submat(4*(i - 1), 0, (4* (i - 1) + 3), 3);


     

     epsilonM.submat(4*i, 0 ,(4*i) + 3, 3) = epsilonM.submat( 4*(i - 1), 0, (4*(i - 1 ) + 3 ), 3) + alphaM.submat(4*(i - 1 ), 0, (4*( i - 1) + 3), 3) * inv(EnM[n] * id - epsilonM.submat(4*(i - 1 ), 0, (4 *(i - 1) + 3 ), 3) + ii * eta * id) * betaM.submat( 4*(i - 1), 0, ( 4*(i - 1 ) + 3), 3) + betaM.submat( 4*(i - 1), 0, (4*(i - 1 ) + 3), 3 ) * inv(EnM[n] *id - epsilonM.submat( 4*(i - 1), 0, (4*(i - 1 ) + 3) , 3 ) + ii * eta * id  ) * alphaM.submat( 4*(i - 1), 0, (4*(i - 1 ) + 3), 3 ); 
  
		    		 

     epsilonSM.submat(4*i, 0, (4*i) + 3, 3) = epsilonSM.submat( 4*(i - 1), 0, (4*(i - 1 ) + 3) , 3 ) + alphaM.submat( 4*(i - 1), 0, (4*(i - 1 ) + 3), 3 )  * inv(EnM[n] * id - epsilonM.submat( 4*(i - 1), 0, (4*(i - 1 ) + 3), 3 ) + ii * eta * id ) * betaM.submat( 4*(i - 1), 0, (4*(i - 1 ) + 3),3 );

     

   }






 // use submatrix to fill fineSM

 // Dimension of epsilonSM is 4M x 4. After all M iterations have been ran, we take the last 4 x 4 matrix in that large matrix and assign it a place in fineSM. The last row of epsilonSM is 4M - 1. 

 fineSM.submat(4*n , 0, (4*n + 3), 3) = epsilonSM.submat(4*M - 4, 0, (4*M -1) , 3 ); 

 


 // Calculate the Green function
  
 G00tab.submat(4*n, 0, (4*n + 3) , 3) = inv(EnM[n] * id - fineSM.submat(4*n, 0, (4*n + 3) ,3 ) + ii * eta * id  );


 


 
    
// Sxx(Q,w)

 
 complex<double> al;
 complex<double> alDag;
 complex<double> gam;
 complex<double> gamDag;
 complex<double> bet;
 complex<double> betDag;



al = (cos(QQ[u]*a) + ii * sin(QQ[u]*a ) )*(1.0 - (Jt/Jz)*cos( QQ[u]* a ) );              

alDag = (cos(QQ[u]*a) - ii * sin(QQ[u]*a ) )*(1.0 - (Jt/Jz)*cos( QQ[u]* a ) );

gam = (cos(QQ[u]*a) + ii * sin(QQ[u]*a) )*(-Jt/(2.0*Jz))*(1.0 + (cos(2*QQ[u]*a) + ii * sin(2*QQ[u]*a)  ));
 
gamDag = (cos(QQ[u]*a) - ii * sin(QQ[u]*a ) )* (-Jt/(2.0*Jz))*(1.0 + (cos(2*QQ[u]*a) - ii * sin(2*QQ[u]*a) ) );

bet = ((-gfac * mu * B)/(2.0*Jz))*(cos(QQ[u]*a) + ii * sin( QQ[u]*a )    )* ( 1.0 + cos(QQ[u]* a) + ii * sin(QQ[u]*a )   );
							   
betDag = ((-gfac * mu * B)/(2.0*Jz))*(cos(QQ[u]*a) - ii * sin( QQ[u]*a )    )* ( 1.0 + cos(QQ[u]* a) - ii * sin(QQ[u]*a )   );



 SxxNum3D(u,n) =  ((-1.0)/(4.0*(PI)))* imag(  alDag*al* G00tab.submat(4*n , 0, (4*n + 3)  , 3)(0,0)   +  alDag*gam* G00tab.submat(4*n , 0, (4*n + 3) , 3)(0,2)   + alDag*bet * G00tab.submat(4*n , 0,( 4*n  + 3) , 3 )(0,1) + al*gamDag*G00tab.submat(4*n , 0,( 4*n +3) , 3 )(2,0) + al*betDag* G00tab.submat(4*n , 0,( 4*n + 3) , 3  )(1,0) );    





 

 
 


 // Syy(Q,w)



 complex<double> alph;
 complex<double> alphDag;
 complex<double> beta;
 complex<double> betaDag;
 complex<double> del; 
 complex<double> delDag; 
 


 alph = (ii/2.0)*(cos(QQ[u]*a) + ii*sin(QQ[u]*a))*(1.0 + (-Jt/(2.0*Jz))*(2.0*cos(QQ[u]*a)) );

 alphDag = (-ii/2.0)*(cos(QQ[u]*a) - ii*sin(QQ[u]*a))*(1.0 + (-Jt/(2.0*Jz))*(2.0*cos(QQ[u]*a)) );


 beta =  (ii/2.0)*(cos(QQ[u]*a) + ii*sin(QQ[u]*a))*(-Jt/(2.0*Jz))*(cos(2.0*QQ[u]*a) + ii*sin(2.0*QQ[u]*a)  + 1.0 );

 betaDag =  (-ii/2.0)*(cos(QQ[u]*a) - ii*sin(QQ[u]*a))*(-Jt/(2.0*Jz))*(cos(2.0*QQ[u]*a) - ii*sin(2.0*QQ[u]*a)  + 1.0 );



 del =  (ii/2.0)*(cos(QQ[u]*a) + ii*sin(QQ[u]*a))*(cos(QQ[u]*a) + ii*sin(QQ[u]*a)   - 1.0)*((gfac * mu * B)/(2.0*Jz) ) ;  

 delDag =  (-ii/2.0)*(cos(QQ[u]*a) - ii*sin(QQ[u]*a))*(cos(QQ[u]*a) - ii*sin(QQ[u]*a)   - 1.0)*((gfac * mu * B)/(2.0*Jz) ) ;  





 SyyNum3D(u,n) = ((-1.0)/PI)*imag(alphDag*alph* G00tab.submat(4*n , 0,( 4*n + 3) , 3  )(0,0)   +  alphDag*beta* G00tab.submat(4*n , 0,( 4*n + 3) , 3  )(0,2)  + alphDag*del* G00tab.submat(4*n , 0,( 4*n + 3) , 3  )(0,1)  + betaDag*alph* G00tab.submat(4*n , 0,( 4*n + 3) , 3  )(2,0)  + delDag*alph* G00tab.submat(4*n , 0,( 4*n + 3) , 3  )(1,0)         );








 

 // Syy for |0_2>






 complex<double> aaa;
 complex<double> bbb;
 complex<double> ccc;

 complex<double> aaaDag;
 complex<double> bbbDag;
 complex<double> cccDag;


 aaa = (-1.0)*( (ii/2.0)*exp(ii*QQ[u]) + exp(ii*QQ[u])*(ii/2.0)*(-Jt/(2.0*Jz))*(exp(-ii*QQ[u])) + exp (ii*QQ[u])    );

 aaaDag = (-1.0)*( (-ii/2.0)*exp(-ii*QQ[u]) + exp(-ii*QQ[u])*(-ii/2.0)*(-Jt/(2.0*Jz))*(exp(ii*QQ[u])) + exp (-ii*QQ[u])    );



 bbb = (-1.0)*( (ii/2.0)*exp(ii*QQ[u])*(-Jt/(2.0*Jz))*(1.0 + exp(ii*2.0*QQ[u]))  );

 bbbDag = (-1.0)*( (-ii/2.0)*exp(-ii*QQ[u])*(-Jt/(2.0*Jz))*(1.0  + exp(-ii*2.0*QQ[u]))  );


 ccc = (-1.0)*((ii/(2.0))*exp(ii*QQ[u])*((gfac * mu * B)/(2.0*Jz) )*(exp(ii*QQ[u])  -  1.0    )  );


 






















 // Szz(Q,w)



 complex<double> chi;
 complex<double> xixi;

 complex<double> chiDag;
 complex<double> xixiDag;


 chi = (+Jt/(2.0*Jz))*(1 - cos(QQ[u]*a) - ii*sin(QQ[u]*a))*(cos(QQ[u]*a) + ii*sin(QQ[u]*a)  );

 chiDag = (+Jt/(2.0*Jz))*(1 - cos(QQ[u]*a) + ii*sin(QQ[u]*a))*(cos(QQ[u]*a) - ii*sin(QQ[u]*a)  );


 xixi = ((+gfac * mu * B)/(2.0*Jz) )*( cos(QQ[u]*a) + ii*sin(QQ[u]*a) );

 xixiDag  = ((+gfac * mu * B)/(2.0*Jz) )*( cos(QQ[u]*a) - ii*sin(QQ[u]*a) );



 SzzNum3D(u,n) = (-1.0/PI)*imag( chiDag*chi* G00tab.submat(4*n , 0,( 4*n + 3) , 3  )(1,1)    +  chiDag*xixi* G00tab.submat(4*n , 0,( 4*n + 3) , 3  )(1,0)  +  xixiDag*chi* G00tab.submat(4*n , 0,( 4*n + 3) , 3  )(0,1)  +  xixiDag*xixi* G00tab.submat(4*n , 0,( 4*n + 3) , 3  )(0,0) );


  





  
 


 
 
 // Syz(Q,w) - Szy(Q,w) and Syz(Q,w) + Szy(Q,w)



 complex<double> ch;
 complex<double> xi;
 complex<double> zet;
 complex<double> phi;
 complex<double> thet;



 complex<double> chDag;
 complex<double> xiDag;
 complex<double> zetDag;
 complex<double> phiDag;
 complex<double> thetDag;


  

// Normal DM

 ch = ii*(0.5)*( cos(QQ[u]*a) + ii * sin(QQ[u]*a ))*( 1.0 + ((- Jt)/(2.0*Jz))*(cos(QQ[u]*a) + ii * sin(QQ[u]*a) + cos(QQ[u]*a) - ii*sin(QQ[u]*a ) ));    
 xi = ii*(0.5)*(cos(QQ[u]*a) + ii*sin(QQ[u]*a) )*((-Jt)/(2.0*Jz))*(cos(2.0*QQ[u]*a) + ii*sin(2.0*QQ[u]*a) + 1.0);
 zet = ii*(0.5)*(cos(QQ[u]*a) + ii*sin(QQ[u]*a) )*(cos(QQ[u]*a) + ii*sin(QQ[u])  - 1.0)*((gfac* mu * B)/(2.0*Jz));
 phi = (+Jt/(2.0*Jz))*(1.0 - cos(QQ[u]*a) - ii*sin(QQ[u]*a)  )*(cos(QQ[u]*a) + ii*sin(QQ[u]*a) );
 thet = (+(gfac * mu * B)/(2.0*Jz))*(cos(QQ[u]*a) + ii*sin(QQ[u]*a)     );

 chDag = -ii*(0.5)*( cos(QQ[u]*a) - ii * sin(QQ[u]*a ))*( 1.0 + ((- Jt)/(2.0*Jz))*(cos(QQ[u]*a) - ii * sin(QQ[u]*a) + cos(QQ[u]*a) + ii*sin(QQ[u]*a) ));    
 xiDag = -ii*(0.5)*(cos(QQ[u]*a) - ii*sin(QQ[u]*a) )*((-Jt)/(2.0*Jz))*(cos(2.0*QQ[u]*a) - ii*sin(2.0*QQ[u]*a) + 1.0 );
 zetDag = -ii*(0.5)*(cos(QQ[u]*a) - ii*sin(QQ[u]*a) )*(cos(QQ[u]*a) - ii*sin(QQ[u])  - 1.0)*( (gfac* mu * B)/(2.0*Jz)); 
 phiDag = ((+Jt)/(2.0*Jz))*(1.0 - cos(QQ[u]*a) + ii*sin(QQ[u]*a)  )*(cos(QQ[u]*a) - ii*sin(QQ[u]*a) );
 thetDag = (+(gfac * mu * B)/(2.0*Jz))*(cos(QQ[u]*a) - ii*sin(QQ[u]*a)     ); 




 //Now calculating cross term



  // Just the first order terms.

  SyzzyNum3D(u,n) = (1.0/PI)*imag( ii* ( chDag*phi*G00tab.submat(4*n, 0, (4*n + 3), 3   )(0,1) + chDag*thet*G00tab.submat(4*n, 0, (4*n + 3), 3   )(0,0)  )  -   ii*( phiDag*ch*G00tab.submat(4*n, 0, (4*n + 3), 3  )(1,0) +  thetDag*ch*G00tab.submat(4*n, 0, (4*n +3),  3 )(0,0)  ));



  
     


  

 // For this calculation, the energy is fixed and the values Q_{||} (QQ[u]) are evaluated. 


 for(int t = 0; t < F; t++)
   {



   
        QP[t] = ( PerpU * (t + 0.00001)/(F-1))*5.5; // 20.01.2015 


       // 22.02.2015. For Qperp Q par

       SdotSPerp(u,t) =        (  1.0 -  ( ( QP[t]*QP[t] )/(QQ[u]*QQ[u] + QP[t]*QP[t])   )     )*2.0*SxxNum3D(u,n) + 2.0* SyyNum3D(u,n) + (  1.0   -  ( ( QQ[u]*QQ[u] )/( QQ[u]*QQ[u] + QP[t]*QP[t]))   )*2.0*SzzNum3D(u,n);   // Factor of 2.0 is introduced to account for |0_1> and |0_2>; Szxxz is omitted as it has zero value when both |0_1> and |0_2> are taken into account. 





    

   }






                    



























  

   // For the case of plotting Q_{||} and omega. Fixed Q_perp and Q_||

    

 for(int t = 0; t < F; t++)
   {

    


    


              QP[t] = ( PerpU * (t + 0.00001)/(F-1))*5.5; // 20.01.2015  takes into account explicit lattice constants 

     

              SyzzyPerpNum3DPlus(u,t) = (+1.0)*(1.0 - ( ( QQ[u]*QQ[u] )/( QQ[u]*QQ[u] + QP[t]*QP[t] ) )) *2.0* SyzzyNum3D(u,n); // Sxy - Syx is omitted since in the |0_1> and |0_2> regime this term vanishes. Also, a factor of 2.0 is included to account for two topological sectors. 


       // 22.02.2015 Qpar Qperp

         SyzzyPerpNum3DMinus(u,t) = (-1.0)*(1.0 - ( ( QQ[u]*QQ[u] )/( QQ[u]*QQ[u] + QP[t]*QP[t] ) )) *2.0* SyzzyNum3D(u,n); // Sxy - Syx is omitted since in the |0_1> and |0_2> regime this term vanishes. Also, a factor of 2.0 is included to account for two topological sectors. 




   }
 

   



 









 
 







 

 

 


 

    



 

 // Inclusion of the form factor will depend mainly on the absolute value of Q_|| and Q_perp


     for(int t = 0; t < F; t++)
       {



	





	    TOTALQ_PARQ_PERP_MINUS_PLUS(u,t) = (formfac(1.5, 1.0, 2.5, QQ[u], QP[t]))*(formfac(1.5, 1.0, 2.5, QQ[u], QP[t]))* (  SdotSPerp(u,t) + SyzzyPerpNum3DPlus(u,t) +  SdotSPerp(u,t) + SyzzyPerpNum3DMinus(u,t) );   // Just  I+ + I- and form factor.


	 

       }



    
 
 

















  

 }




   }





 // Out I+ for comparison with Nat Phys Fig 3


 ofstream NPP;
 NPP.open("NPP.dat");



 for (int n = 0; n < J; n++)
   {

     // 25.02.2015

     //  NPP << EnM[n] << " " << NatPhysScanPlus[n] << endl;


   }




 // Out I- for comparison with Nat Phys Fig 3


 ofstream NPM;
 NPM.open("NPM.dat");



 for (int n = 0; n < J; n++)
   {

     // 25.02.2015

     // NPM << EnM[n] << " " << NatPhysScanMinus[n] << endl;


   }






 // Out (I+ - I-)/(I+ + I-)



 
 ofstream NPAM;
 NPAM.open("NPAM.dat");



 for (int n = 0; n < J; n++)
   {


     // 25.02.2015

     //  NPAM << EnM[n] << " " << (NatPhysScanPlus[n] - NatPhysScanMinus[n])/(NatPhysScanPlus[n] + NatPhysScanMinus[n]) << endl;


   }








 // Form Factor

 ofstream FF;
 FF.open("FF.dat");             



 for(int u = 0;  u < H; u++)
   {

      for(int t = 0; t < F; t++)

	{

	  FF <<  QQ[u] << " " << QP[t] <<  " "  <<  (formfac(1.5, 1.0, 2.5, QQ[u], QP[t]))*(formfac(1.5, 1.0, 2.5, QQ[u], QP[t])) << endl;;


       } 



   }







//Sxx(Q,w)

 // The final stage is to output the above data to a .dat file. The preferable way for Mathematica to read this set is with each element as a vector: x is QQ, y is EnM, z is SxxNum3D. This is simply achieved by outputting the information in this order, with a whitespace separating each data point.  







 // For |0_1>

 ofstream outfileSxx;                                // Making an instance of ofstream class, that is, an object.                                                                                            
 outfileSxx.open("green_Sxx_3D_B_perp_DM_numerical.dat");                        // This file should be initially empty.                                                                                        



 // For |0_2>

 ofstream outfileSxx02; 

 outfileSxx02.open("green_Sxx_02_3D_B_perp_DM_numerical.dat"); 





 // For |0_1>  and |0_2>

 ofstream outfileSxx01and02; 

 outfileSxx01and02.open("green_Sxx_02and01_3D_B_perp_DM_numerical.dat"); 

 // outfileSxx01and02.open("testHang.dat"); 






 ofstream outfileSxxLog;

 outfileSxxLog.open("green_Sxx_Log_3D_B_perp_DM_numerical.dat");                        // This file should be initially empty.                                                                                                          


 // Now you are interested in generating contour plots in Mathematica. I think the best way to go about this is to out three DIFFERENT lists. Then, in Mathematica, assign a name to each of the three lists. From there you can use ListContourPlot. Actually, wait; ListContourPlot takes in the heights at a given point in the Q omega plane. So all we need to do is out SxxNum3D(u,n) and run through u and n 

 // ofstream outfileContour;

 // outfileContour.open("green_3D_B_perp_DM_Contour_numerical.dat" );


 

 for(int u = 0; u < H; u++ )
   {



     for(int n = 0; n < J; n++)
       {



	
	  outfileSxx01and02 << QQ[u] << " " << EnM[n] << " " << 2.0* SxxNum3D(u,n) << endl;   // Total response (|0_1> and |0_2>)

	 



       }





   }





   















  
 

//Syy(Q,w)

 // The final stage is to output the above data to a .dat file. The preferable way for Mathematica to read this set is with each element as a vector: x is QQ, y is EnM, z is SxxNum3D. This is simply achieved by outputting the information in this order, with a whitespace separating each data point.  








 // For |0_1>

 ofstream outfileSyy;                                // Making an instance of ofstream class, that is, an object.                                                                                        

 outfileSyy.open("green_Syy_3D_B_perp_DM_numerical.dat");                        // This file should be initially empty.                                                                                                 



 // For |0_2>

 ofstream outfileSyy02;                                // Making an instance of ofstream class, that is, an object.                                                                                        

 outfileSyy02.open("green_Syy_02_3D_B_perp_DM_numerical.dat");                        // This file should be initially empty.                                                                                                 






 // For |0_1> and |0_2>

 ofstream outfileSyy01and02;                                // Making an instance of ofstream class, that is, an object.                                                                                        

 outfileSyy01and02.open("green_Syy_01and02_3D_B_perp_DM_numerical.dat");                        // This file should be initially empty.                                                                                                 







 ofstream outfileSyyLog;    

 outfileSyyLog.open("green_Syy_Log_3D_B_perp_DM_numerical.dat");                        // This file should be initially empty.                                                                                                 


 // Now you are interested in generating contour plots in Mathematica. I think the best way to go about this is to out three DIFFERENT lists. Then, in Mathematica, assign a name to each of the three lists. From there you can use ListContourPlot. Actually, wait; ListContourPlot takes in the heights at a given point in the Q omega plane. So all we need to do is out SxxNum3D(u,n) and run through u and n 

 // ofstream outfileContour;

 // outfileContour.open("green_3D_B_perp_DM_Contour_numerical.dat" );




 for(int u = 0; u < H; u++ )
   {



     for(int n = 0; n < J; n++)
       {


 
	 
	   outfileSyy01and02 << QQ[u] << " "  << EnM[n] << " "  << 2.0* SyyNum3D(u,n) << endl;   // Total response (|0_1> and |0_2>)





       }





   }



   
 





//Szz(Q,w)

 // The final stage is to output the above data to a .dat file. The preferable way for Mathematica to read this set is with each element as a vector: x is QQ, y is EnM, z is SxxNum3D. This is simply achieved by outputting the information in this order, with a whitespace separating each data point.  








 // For |0_1> 

 ofstream outfileSzz;                                // Making an instance of ofstream class, that is, an object.                                                                                        

 outfileSzz.open("green_Szz_3D_B_perp_DM_numerical.dat");                        // This file should be initially empty.                                                                                                 




 // For |0_2>

 ofstream outfileSzz02;                                // Making an instance of ofstream class, that is, an object.                                                                                        
 outfileSzz02.open("green_Szz_02_3D_B_perp_DM_numerical.dat");                        // This file should be initially empty.                                                                                                 



 // For |0_1> and |0_2>

 ofstream outfileSzz01and02;                                // Making an instance of ofstream class, that is, an object.                                                                                        
 outfileSzz01and02.open("green_Szz_01and02_3D_B_perp_DM_numerical.dat");                        // This file should be initially empty.                                                                                                 






 ofstream outfileSzzLog;                                // Making an instance of ofstream class, that is, an object.                                                                                            
 outfileSzzLog.open("green_Szz_Log_3D_B_perp_DM_numerical.dat");                        // This file should be initially empty.                                                                                                 


 // Now you are interested in generating contour plots in Mathematica. I think the best way to go about this is to out three DIFFERENT lists. Then, in Mathematica, assign a name to each of the three lists. From there you can use ListContourPlot. Actually, wait; ListContourPlot takes in the heights at a given point in the Q omega plane. So all we need to do is out SxxNum3D(u,n) and run through u and n 

 // ofstream outfileContour;

 // outfileContour.open("green_3D_B_perp_DM_Contour_numerical.dat" );




 for(int u = 0; u < H; u++ )
   {



     for(int n = 0; n < J; n++)
       {


 
 
	    outfileSzz01and02 << QQ[u] << " "  << EnM[n] << " "  << 2.0* SzzNum3D(u,n) << endl;  // Total response (|0_1> and |0_2>)

	 



       }





   }


 

 
 //Syzzy(Q,w)



 // For |0_1>


 ofstream outfileSyzzy;                                // Making an instance of ofstream class, that is, an object.                                                                                            
 outfileSyzzy.open("green_Syzzy_3D_B_perp_DM_numerical.dat");                        // This file should be initially empty.                                                                                                 




 // For |0_2>

 ofstream outfileSyzzy02;                                // Making an instance of ofstream class, that is, an object.                                                                                            
 outfileSyzzy02.open("green_Syzzy_02_3D_B_perp_DM_numerical.dat");                        // This file should be initially empty.                                                                                                 




 // For |0_1> and |0_2>

 ofstream outfileSyzzy01and02;                                // Making an instance of ofstream class, that is, an object.                                                                                            
 outfileSyzzy01and02.open("green_Syzzy_01and02_3D_B_perp_DM_numerical.dat");                        // This file should be initially empty.                                                                                                 

 







 

 ofstream outfileSyzzyLog;                                // Making an instance of ofstream class, that is, an object.                                                                                              
 outfileSyzzyLog.open("green_Syzzy_Log_3D_B_perp_DM_numerical.dat");                        // This file should be initially empty.                                                                                                 










 


 ofstream outfileSyzzyS;    // Symmetric


 ofstream outfileSyzzySLog;  


 outfileSyzzyS.open("green_SyzzyS_3D_B_perp_DM_numerical.dat");                        // This file should be initially empty.                                                                                                 

 outfileSyzzySLog.open("green_SyzzyS_Log_3D_B_perp_DM_numerical.dat");                        // This file should be initially empty.                                                                                                 







 
 for(int u = 0; u < H; u++ )
   {

 

     for(int n = 0; n < J; n++)
       {




	    outfileSyzzy01and02 << QQ[u] << " "  << EnM[n] << " "  << 2.0* SyzzyNum3D(u,n) << endl;  // Total response (|0_1> and |0_2>)



        }




   }


 

 

 


 // Outputting I+ - I-/(I+ + I-)


 ofstream outfileTOTALQ_PARQ_PERP_MINUS_PLUS;                                // Making an instance of ofstream class, that is, an object.                                                                                            

 outfileTOTALQ_PARQ_PERP_MINUS_PLUS.open("green_TOTALQ_PARQ_PERP_MINUS_PLUS_3D_B_perp_DM_numerical.dat");                        // This file should be initially empty.                                                                                                 





 for(int u = 0; u < H; u++ )
   {



     for(int t = 0; t < F; t++)
       {


	 outfileTOTALQ_PARQ_PERP_MINUS_PLUS << QQ[u] << " " << QP[t] << " " << TOTALQ_PARQ_PERP_MINUS_PLUS(u,t) << endl;
       


       }





   }







 time2 = clock();


 double ttime;

 ttime = ((double)(time2 -time1))/CLOCKS_PER_SEC;

 cout << ttime << " sec  "   << endl;

 return 0;


}
