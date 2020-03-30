import numpy as np
def foxWo(particles):
    entries=particles.size()
    h0=0.0
    hd=0.0
    
    # init
    h10=0.0
    h20=0.0
    h30=0.0
    h40=0.0  
    
    if entries<2:
        print("FOXWO :Not enough entries")
        h10=-1.0
        h20=-1.0
        h30=-1.0
        h40=-1.0
    else:
        particleArray = np.zeros((entries, 4)) 
        
        # copy all particles into an array and calculate h0 and hd
        for i in range(entries):
            momentum = particles[i].getMomentum()
            particleArray[i][0]=momentum[0]
            particleArray[i][1]=momentum[1]
            particleArray[i][2]=momentum[2]
            #print(i, momentum[0], momentum[1], momentum[2])
            magsquared=momentum[0]*momentum[0]+momentum[1]*momentum[1]+momentum[2]*momentum[2]
            particleArray[i][3]=np.sqrt(magsquared)
            h0=h0+particleArray[i][3]
            hd=hd+magsquared
        #print(h0, hd, particleArray.shape)
        
        h0=h0*h0 # square it
        # now do the momenta
        for i in range(entries):
            for j in range(i+1, entries):
                ptemp=particleArray[i][3]*particleArray[j][3]
                cthe=(particleArray[i][0]*particleArray[j][0] +
                particleArray[i][1]*particleArray[j][1] +
                particleArray[i][2]*particleArray[j][2])/ptemp
                h10 += ptemp*cthe
                h20 += ptemp*(1.5*cthe*cthe-0.5)
                h30 += ptemp*(2.5*cthe*cthe*cthe-1.5*cthe)
                h40 += ptemp*(4.375*cthe*cthe*cthe*cthe-3.75*cthe*cthe+0.375)
        
        # normalize the momenta
        h10=(hd+2*h10)/h0
        h20=(hd+2*h20)/h0
        h30=(hd+2*h30)/h0
        h40=(hd+2*h40)/h0    
    return h10, h20, h30, h40


# a simple implementation avoiding vectors should make it simpler to port to Python and Julia
def CrossProduct(v, i, j, result): 
    result[0]=v[i][1]*v[j][2]-v[i][2]*v[j][1]
    result[1]=v[i][2]*v[j][0]-v[i][0]*v[j][2]
    result[2]=v[i][0]*v[j][1]-v[i][1]*v[j][0]
    result[3]=0
    return 0


def Thrust(particles):
    entries=particles.size()
    
    FullVector = np.zeros(4)
    PartVector = np.zeros(4)
    RefVector = np.zeros(4)
    SumVector = np.zeros(4)
    MaxVector = np.zeros(4)
    thrust = 0.0
        
    if entries<2:
        print("Thrust : Not enough entries")
        thrust=-1
    else:
        particleArray = np.zeros((entries, 4))
        RefVectorMag=0
        dotProduct=0
        magsquared=0
        c1=0
        c2=0
        # copy all particles into an array and start calcualting the thrust
        for i in range(entries):
            momentum = particles[i].getMomentum()
            particleArray[i][0]=momentum[0]
            particleArray[i][1]=momentum[1]
            particleArray[i][2]=momentum[2]
            # print(i, momentum[0], momentum[1], momentum[2])
            magsquared= momentum[0]*momentum[0]+momentum[1]*momentum[1]+momentum[2]*momentum[2]
            particleArray[i][3]=np.sqrt(magsquared) # store the abs(p)
            
            SumVector[0]+=particleArray[i][0]
            SumVector[1]+=particleArray[i][1]
            SumVector[2]+=particleArray[i][2]
            SumVector[3]+=particleArray[i][3]

        # now do the thrust        
        for i in range(entries):
            for j in range(i+1, entries):
                CrossProduct(particleArray,i,j,RefVector) # store the cross product in the ref vector
                RefVectorMag = np.sqrt(RefVector[0]*RefVector[0]+RefVector[1]*RefVector[1]+RefVector[2]*RefVector[2])

                RefVector /= RefVectorMag
                PartVector *= 0
   
                # Add all momenta with sign; two choices for each reference particle.
                for k in range(entries):
                    if k != i and k != j:
                        dotProduct = particleArray[k, 0]*RefVector[0] + particleArray[k, 1]*RefVector[1] + particleArray[k, 2]*RefVector[2]
                        if dotProduct > 0.:
                            PartVector += particleArray[k, :]
                        else:
                            PartVector -= particleArray[k, :]                
                c1=0
                c2=0
                for l in range(4):
                    if l == 0:
                        c1=1
                        c2=1                 
                    elif l == 1:
                        c1=1
                        c2=-1
                    elif l == 2:
                        c1=-1
                        c2=1
                    else:         
                        c1=-1
                        c2=-1
                    
                    FullVector[0:3] = PartVector[0:3] + c1*particleArray[i, 0:3] + c2*particleArray[j, 0:3]
                    FullVector[3]=np.sqrt(FullVector[0]*FullVector[0]+FullVector[1]*FullVector[1]+FullVector[2]*FullVector[2])
                    if FullVector[3] > MaxVector[3]:
                        MaxVector[:] = FullVector[:]
        thrust = MaxVector[3] / SumVector[3]
    return thrust
