#include "libFoxWo.hh"

using namespace std;
using namespace lcio;

int foxwo(vector <ReconstructedParticle *> & particles,double & h10,double & h20, double & h30, double & h40)
{     
    int entries=particles.size();
    double h0=0.0;
    double hd=0.0;
    
    double **particleArray;
    // init
    h10=0;
    h20=0;
    h30=0;
    h40=0;
    
    
    
    if (entries<2)
    {
        cerr << ("FOXWO :Not enough entries")<<endl;;
        h10=-1;
        h20=-1;
        h30=-1;
        h40=-1;
    }
    else 
    {
        particleArray =new double*[entries];
        for (int n=0; n<entries; n++)
        {    
            particleArray[n]=new double[4];
        }    
        
        
        const double * momentum; 
        double magsquared=0;
        double ptemp=0;
        double cthe=0;
        
        //copy all particles into an array and calculate h0 and hd
        for (int i=0; i< entries; i++)
        {            
            momentum = particles[i]->getMomentum();
            particleArray[i][0]=momentum[0];
            particleArray[i][1]=momentum[1];
            particleArray[i][2]=momentum[2];
            //cout << i << " "<<  momentum[0]<< " "<<  momentum[1]<< " "<<  momentum[2]<<endl;
            magsquared=momentum[0]*momentum[0]+momentum[1]*momentum[1]+momentum[2]*momentum[2];
            particleArray[i][3]=sqrt(magsquared);
            h0=h0+particleArray[i][3];
            hd=hd+magsquared;
            
        }
        
        h0=h0*h0; //square it
        // now do the momenta
        for (int i=0; i<entries; i++)
        {
            for (int j=i+1; j<entries; j++)
            {
                ptemp=particleArray[i][3]*particleArray[j][3];
                cthe=(particleArray[i][0]*particleArray[j][0] +
                particleArray[i][1]*particleArray[j][1] +
                particleArray[i][2]*particleArray[j][2])/ptemp;
                h10=h10+ptemp*cthe;
                h20=h20+ptemp*(1.5*cthe*cthe-0.5);
                h30=h30+ptemp*(2.5*cthe*cthe*cthe-1.5*cthe);
                h40=h40+ptemp*(4.375*cthe*cthe*cthe*cthe-3.75*cthe*cthe+0.375);
                
            }
        }
        
        //normalize the momenta
        h10=(hd+2*h10)/h0;
        h20=(hd+2*h20)/h0;
        h30=(hd+2*h30)/h0;
        h40=(hd+2*h40)/h0;
        //free memory
        for( int i = 0 ; i < entries ; i++ )
        {
            //first we delete each row
            delete [] particleArray[i] ;
        }
        //finally, we delete the array of pointers
        delete [] particleArray ;
    } //end of else
    
    
    
    return 0;
}

//a simple implementation avoiding vectors should make it simpler to port to Python and Julia
int CrossProduct(double **v,int i, int j, double * result)
{    
    result[0]=v[i][1]*v[j][2]-v[i][2]*v[j][1];
    result[1]=v[i][2]*v[j][0]-v[i][0]*v[j][2];
    result[2]=v[i][0]*v[j][1]-v[i][1]*v[j][0];
    result[3]=0;
    return 0;
}



int Thrust(vector <ReconstructedParticle *> & particles,double & thrust)
{    
    int entries=particles.size();  
    double **particleArray;
    
    double FullVector[4], PartVector[4], RefVector[4],SumVector[4],MaxVector[4];
    
    //initialize the vectors 
    for (int i=0; i<4; i++)
    {
        FullVector[i]=PartVector[i]=RefVector[i]=SumVector[i]=MaxVector[i]=0;
    }
    
    if (entries<2)
    {
        cerr << ("Thrust : Not enough entries")<<endl;
        thrust=-1;
    }
    else 
    {
        particleArray =new double*[entries];
        for (int n=0; n<entries; n++)
        {    
            particleArray[n]=new double[4];
        }    
        
        
        const double * momentum;
        double RefVectorMag=0; 
        double dotProduct=0;
        double magsquared=0;
        int c1=0; 
        int c2=0;
        //copy all particles into an array and start calcualting the thrust
        for (int i=0; i< entries; i++)
        {            
            momentum = particles[i]->getMomentum();
            particleArray[i][0]=momentum[0];
            particleArray[i][1]=momentum[1];
            particleArray[i][2]=momentum[2];
            //cout << i << " "<<  momentum[0]<< " "<<  momentum[1]<< " "<<  momentum[2]<<endl;
            magsquared= momentum[0]*momentum[0]+momentum[1]*momentum[1]+momentum[2]*momentum[2];
            particleArray[i][3]=sqrt(magsquared); //store the abs(p)
            
            SumVector[0]+=particleArray[i][0];
            SumVector[1]+=particleArray[i][1];
            SumVector[2]+=particleArray[i][2];
            SumVector[3]+=particleArray[i][3];
            
        }
        
        // now do the thrust
        
        for (int i=0; i<entries-1; i++)
        {
            for (int j=i+1; j<entries; j++)
            {
                CrossProduct(particleArray,i,j,RefVector); //store the cross product in the ref vector
                RefVectorMag=sqrt(RefVector[0]*RefVector[0]+RefVector[1]*RefVector[1]+RefVector[2]*RefVector[2]);
                
                for (int n=0; n<4; n++)
                {
                    RefVector[n]=RefVector[n]/RefVectorMag;
                    PartVector[n]=0; 
                }
   
                // Add all momenta with sign; two choices for each reference particle.
                for (int k = 0; k < entries; k++) 
                {
                    if (k != i && k != j) 
                    {
                       dotProduct=particleArray[k][0]*RefVector[0]+particleArray[k][1]*RefVector[1]+particleArray[k][2]*RefVector[2];
                        if (dotProduct > 0.) 
                        {
                            for (int n=0; n<4; n++)
                            {
                              PartVector[n]+= particleArray[k][n];
                            }
                        }    
                        else
                        {
                            for (int n=0; n<4; n++)
                            {
                              PartVector[n]-= particleArray[k][n];
                            } 
                        }    
                    }
                }
                
                c1=0; 
                c2=0;
                for (int l = 0; l < 4; l++) 
                {
                
                    if(l == 0) 
                    {
                        c1=c2=1;
                    }                    
                    else if (l == 1) 
                    {    
                        c1=1;
                        c2=-1;
                    }    
                    else if (l == 2) 
                    {
                        
                        c1=-1;
                        c2=1;
                    }
                    else             
                    {   
         
                        c1=-1;
                        c2=-1;                       
                    }
                    
                    for (int n=0; n<3; n++) // only copy three of the four as we overwrite the fourth element again ...
                    {
                        FullVector[n]= PartVector[n]+c1*particleArray[i][n]+c2*particleArray[j][n];
                    }
                                     
                    FullVector[3]=sqrt(FullVector[0]*FullVector[0]+FullVector[1]*FullVector[1]+FullVector[2]*FullVector[2]);
                    if (FullVector[3] > MaxVector[3]) 
                    {
                        for (int n=0; n<4; n++) // set the new vector
                        {
                            MaxVector[n]=FullVector[n];
                        }
                    }    
                }
            }
        } 
        thrust = MaxVector[3] / SumVector[3];
//empty the memory    
        for( int i = 0 ; i < entries ; i++ )
        {
    //first we delete each row
            delete [] particleArray[i] ;
        }
        //finally, we delete the array of pointers
        delete [] particleArray ;            
    } //end of else     
    

    return 0; 
}

