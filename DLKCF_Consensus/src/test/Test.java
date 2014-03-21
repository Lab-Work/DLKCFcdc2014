package test;


import org.jblas.DoubleMatrix;
import doubleMatrix.Concat;
import model.D;
import model.RoadModel;
import Jama.*;
import section.Section;
import trueSolution.TriangularTrue;
import trueSolution.TrueSolution;
import filters.Estimation;
import trueSolutionCTM.*;

@SuppressWarnings("unused")
public class Test {

	public static void main(String[] args) {
		int cellsec=28;
		int cellsec1=28;
		int overlapsec=10;
		int overlapsec1=10;
		int Nsecs = 7;
		int Nsecs1 = 7;
		int Shocksec =3;
		int Shocksec1 =3;
		int cells = (Nsecs-1)*(cellsec-overlapsec)+cellsec;
		int cells1 = (Nsecs1-1)*(cellsec1-overlapsec1)+cellsec1;
		int distMeasure=overlapsec-1;

	    double rhoLeft = 0.8;
	    double rhoRight = 0.2;
	    
	    TrueSolutionCTM trueSolutionCTM= new TriangularTrueCTM(rhoLeft,rhoRight, 1d, 1000d/((double)cells),cells,distMeasure,Nsecs,true);	    
	    TrueSolution[] trueSolution=new TrueSolution[Nsecs];
	    for(int i=0;i<Shocksec;i++){
	    	trueSolution[i] = new TriangularTrue(trueSolutionCTM,rhoLeft, rhoLeft,rhoLeft,rhoRight, 1d, 1000d/((double)cells), cellsec, overlapsec, i, Nsecs,1,true,false);
	    }
	    trueSolution[Shocksec] = new TriangularTrue(trueSolutionCTM,rhoLeft, rhoRight,rhoLeft,rhoRight, 1d, 1000d/((double)cells), cellsec, overlapsec, Shocksec, Nsecs,1,true,false);
	    for(int i=Shocksec+1;i<Nsecs;i++){
	    	trueSolution[i] = new TriangularTrue(trueSolutionCTM,rhoRight, rhoRight,rhoLeft,rhoRight, 1d, 1000d/((double)cells), cellsec, overlapsec, i,Nsecs,1,true,false);
	    }
	    if (Nsecs==1){
	    	 trueSolution[0].trueRight=trueSolution[0];
	    	 trueSolution[0].trueLeft=trueSolution[0];
	    }
	    else{
	    	trueSolution[0].trueRight=trueSolution[1];
		    trueSolution[Nsecs-1].trueLeft=trueSolution[Nsecs-2];
		    for(int i=1;i<Nsecs-1;i++){
		    	trueSolution[i].trueLeft=trueSolution[i-1];
		    	trueSolution[i].trueRight=trueSolution[i+1];
		    }
	    }	    
	    TrueSolution[] trueSolution1=new TrueSolution[Nsecs1];
	    for(int i=0;i<Shocksec1;i++){
	    	trueSolution1[i] = new TriangularTrue(trueSolutionCTM,rhoLeft, rhoLeft,rhoLeft,rhoRight, 1d, 1000d/((double)cells1), cellsec1, overlapsec1, i, Nsecs1,1,true,false);
	    }
	    trueSolution1[Shocksec1] = new TriangularTrue(trueSolutionCTM,rhoLeft, rhoRight,rhoLeft,rhoRight, 1d, 1000d/((double)cells1), cellsec1, overlapsec1, Shocksec1, Nsecs1,1,true,false);
	    for(int i=Shocksec1+1;i<Nsecs1;i++){
	    	trueSolution1[i] = new TriangularTrue(trueSolutionCTM,rhoRight, rhoRight,rhoLeft,rhoRight, 1d, 1000d/((double)cells1), cellsec1, overlapsec1, i,Nsecs1,1,true,false);
	    }
	    if (Nsecs1==1){
	    	 trueSolution1[0].trueRight=trueSolution1[0];
	    	 trueSolution1[0].trueLeft=trueSolution1[0];
	    }
	    else{
	    	trueSolution1[0].trueRight=trueSolution1[1];
		    trueSolution1[Nsecs1-1].trueLeft=trueSolution1[Nsecs1-2];
		    for(int i=1;i<Nsecs1-1;i++){
		    	trueSolution1[i].trueLeft=trueSolution1[i-1];
		    	trueSolution1[i].trueRight=trueSolution1[i+1];
		    }
	    }

	    Estimation.exportResult(trueSolutionCTM, trueSolution,trueSolution1, 2000,2, "VaryTrustSensor(M0.0025VS0.09_0.0009)_trueCTM1_AllTrue11_EvenTrust(fixed)_FCConsensus");
	}
	
}
