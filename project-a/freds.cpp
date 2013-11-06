#include <cmath>
#include <fstream>

double M;
const double m=0.1;                     //Set Physical Quantities as global variables
const double g=9.8;
const double l=9.8;

double Ldouble[4][4];

using namespace std;

double DoublePen(double yrkdouble[4], double h)              //Single step of Runge Kutta 4 Method for 4x4 L Matrix
{
    double kdouble1[4]; // for all arrays [1] - theta; [2] - phi; [3] - omega; [4] - lamda
    double kdouble2[4];
    double kdouble3[4];
    double kdouble4[4];
    
    for( int j=0; j<=3; j++)            //For loops go through matrix multiplication
    {
        kdouble1[j]=(Ldouble[j][0]*yrkdouble[0]+Ldouble[j][1]*yrkdouble[1]+Ldouble[j][2]*yrkdouble[2]+Ldouble[j][3]*yrkdouble[3])*h;
    }
    
    for( int j=0; j<=3; j++)
    {
        kdouble2[j]=(Ldouble[j][0]*(yrkdouble[0]+0.5*kdouble1[0])+
                     Ldouble[j][1]*(yrkdouble[1]+0.5*kdouble1[1])+
                     Ldouble[j][2]*(yrkdouble[2]+0.5*kdouble1[2])+
                     Ldouble[j][3]*(yrkdouble[3]+0.5*kdouble1[3]))*h;
    }
    
    for( int j=0; j<=3; j++)
    {
        kdouble3[j]=(Ldouble[j][0]*(yrkdouble[0]+0.5*kdouble2[0])+
                     Ldouble[j][1]*(yrkdouble[1]+0.5*kdouble2[1])+
                     Ldouble[j][2]*(yrkdouble[2]+0.5*kdouble2[2])+
                     Ldouble[j][3]*(yrkdouble[3]+0.5*kdouble2[3]))*h;
    }
    
    for( int j=0; j<=3; j++)
    {
        kdouble4[j]=(Ldouble[j][0]*(yrkdouble[0]+kdouble3[0])+
                     Ldouble[j][1]*(yrkdouble[1]+kdouble3[1])+
                     Ldouble[j][2]*(yrkdouble[2]+kdouble3[2])+
                     Ldouble[j][3]*(yrkdouble[3]+kdouble3[3]))*h;
    }
    
    for( int j=0; j<=3; j++)
    {
        yrkdouble[j]=yrkdouble[j]+(1.0/6.0)*(kdouble1[j]+2*kdouble2[j]+2*kdouble3[j]+kdouble4[j]);
    }
    return 0;
}

double energy2mass(double y[4])         //Energy Method uses Small angle approximations.
{
    
    double KE = 0.5*l*l*(m*y[2]*y[2]+M*(y[2]*y[2]+y[3]*y[3]+2.0*y[2]*y[3]));
    
    double PE = 0.5*g*((m+M)*l*y[0]*y[0]+M*l*y[1]*y[1]);
    
    return KE+PE;
}

int main ()
{
    ofstream doublependulum("/Users/Fred/Dropbox/ProjectA/doublepenG0p1R100.txt");
    
    double dt[5]={0.05, 0.1, 0.2, 0.5, 1};
    
    double t,steps;
    double yrkdouble[4];
    
    M=100*m;
    
    double G=0.1;
    double R=M/m;
    
    Ldouble[0][0]=0.0; // 4x4 LMatrix Set up.
    Ldouble[0][1]=0.0;
	Ldouble[0][2]=1.0;
	Ldouble[0][3]=0.0;
    
    Ldouble[1][0]=0.0;
	Ldouble[1][1]=0.0;
	Ldouble[1][2]=0.0;
	Ldouble[1][3]=1.0;
    
    Ldouble[2][0]=-(R+1.0);
	Ldouble[2][1]=R;
	Ldouble[2][2]=-G;
	Ldouble[2][3]=0.0;
    
    Ldouble[3][0]=R+1.0;
	Ldouble[3][1]=-(R+1.0);
	Ldouble[3][2]=G*(1.0-(1.0/R));
	Ldouble[3][3]=-G/R;
    
    for (int n=0; n<=4; n++)
    {
        steps = 30.0/dt[n];                        //Enough steps to fill 60 seconds
        t=0.0;
        
        yrkdouble[0]=0.1;                       //Reset initial conditions
        yrkdouble[1]=0.0;
        yrkdouble[2]=0.0;
        yrkdouble[3]=0.0;
        
        for (int i=0; i<=steps; i++)
        {
            t+=dt[n];
            DoublePen(yrkdouble, dt[n]);
            doublependulum  << t << '\t' << yrkdouble[0]<< '\t' << yrkdouble[1]<< '\t'
                            << yrkdouble[2]<< '\t' << yrkdouble[3]<< '\t' << energy2mass(yrkdouble)<< endl;
        }
    }
    
    return 0;
}