#include <mpi.h>
#include <stdio.h>
#include <math.h>

 
int main(int argc,char *argv[])
{
    int myid, numprocs;
    MPI_Status status;
 
    MPI_Init(&argc,&argv);
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&myid);

 
    //----------------------------------------------
    //Variables' initialization
    int tmax=1, xmax=1;
    int sizeT=10, sizeX=8;
    double K=0.01, D=0.4, F=1000;
    double dt=(double)tmax/sizeT, dx=(double)xmax/sizeX;
    double dx2inv=1/(dx*dx);
    double *U_old, *U_new;
    int DIM=sizeX+1;
    U_old=new double [DIM];
    U_new=new double [DIM];
    for(int i=0;i<sizeX+1;i++){  
       U_old[i]=0;
       U_new[i]=0;
    }
    //Source:
    U_old[sizeX/2]=F;
    //----------------------------------------------
    //Computational pieces for processors and their boundaries
    int piece=sizeX/numprocs;
    int lbound, rbound;
    if(myid==0)lbound=myid*piece+1;
    else lbound=myid*piece+1;
    rbound=lbound+piece-1;
    if(myid==numprocs-1 && sizeX%2==0) rbound--;
    //----------------------------------------------
   
    for(int t=1;t<=sizeT;t++){
       //Messages' exchange-------------------------------------------------------------
       if(myid>0)MPI_Recv(&U_new[lbound-1],1,MPI_DOUBLE,myid-1,1,MPI_COMM_WORLD,&status);
       if(t>1 && myid<numprocs-1)MPI_Recv(&U_new[rbound+1],1,MPI_DOUBLE,myid+1,2,
                                          MPI_COMM_WORLD,&status);
       //-------------------------------------------------------------------------------
       for(int i=lbound;i<=rbound;i++){
           U_new[i]=(K*(U_old[i+1]-2*U_old[i]+U_old[i-1])*dx2inv+D*U_old[i])*dt+U_old[i];
       }
       for(int i=lbound-1;i<=rbound+1;i++){
           U_old[i]=U_new[i];
       }
	   
       //Messages' exchange--------------------------------------------------------------
       if(myid<numprocs-1)MPI_Ssend(&U_new[rbound],1,MPI_DOUBLE,myid+1,1,MPI_COMM_WORLD);
       if(t<sizeT && myid>0)MPI_Ssend(&U_new[lbound],1,MPI_DOUBLE,myid-1,2,
                                      MPI_COMM_WORLD);
       //--------------------------------------------------------------------------------
    }
	
    if(myid>0)MPI_Send(&U_new[lbound],piece,MPI_DOUBLE,0,3,MPI_COMM_WORLD);
    if(myid==0){
       for(int i=1;i<numprocs;i++){
           MPI_Recv(&U_old[i*piece+1],piece,MPI_DOUBLE,i,3,MPI_COMM_WORLD,&status);
       }
       //endwtime = MPI_Wtime();
       //----------------------------------------------
     /*  //Printing to file
       ofstream ofs("results.txt");
       if(!ofs){cout << "Error!"; exit;}
       for(int i=0;i<sizeX+1;i++){
           ofs << U_old[i] << " ";
       }
       ofs.close();
       //----------------------------------------------
       cout << "The end" << endl;
       cout << "Time: " << endwtime-startwtime << " sec" << endl;
       cin.get();*/
    }
    MPI_Finalize();
    return 0;
}