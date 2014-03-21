package trueSolutionCTM;

import org.jblas.DoubleMatrix;
import org.jfree.data.xy.XYSeries;

import doubleMatrix.GaussianGenerator;
import filters.Filter;
import section.*;
import trueSolution.TrueSolution;
import model.*;


//this TrueSolutionCTM is used to compute the true solution of the entire network
public abstract class TrueSolutionCTM {


	int locDis; //location of the discontinuity	
	public double rhoLeft; //left initial condition 
	public double rhoRight; //right initial condition
	
	public double rhoMax;
	public double speedMax;
	public double rhoCritical;
	
	public double dt;
	public double dx;
	public int cells;
	public int numSections;
			
	public DoubleMatrix trueStatesCTM;//true state of the entire network
	public DoubleMatrix trueStatesCTMPrior;//true state at the previous step
	public DoubleMatrix measurementCTM;
	public boolean measurePerturb;//if the system has low performance sensors	
	public DoubleMatrix measureVarDensity;
	public GaussianGenerator measureGeneratorDensity;
	public int distMeasurement;
	
	public int numUp;

	abstract public double speed(double density);
	abstract public void propagateCTM();
	abstract public void getMeasurementCTM();	
	abstract public double sending(double density);
	abstract public double receiving(double density);
	abstract public double computeFlux(double density1, double density2);
	
	public void initial(double _rhoLeft, double _rhoRight, double _dt, double _dx, int _cells, int _distMeasurement, int _numSections, boolean _measurePerturb) {

		rhoLeft = _rhoLeft;
		rhoRight = _rhoRight;
		
		dx = _dx;
		dt = _dt;
		cells=_cells;
		numSections=_numSections;		
		
		distMeasurement=_distMeasurement;		
		measurePerturb=_measurePerturb;
		
		trueStatesCTM=DoubleMatrix.zeros(cells,1);
		
		locDis = (int)cells/2; //location of the discontinuity

		for(int i=0;i<locDis;i++){
			trueStatesCTM.put(i,0,rhoLeft);
		}
		for(int i=locDis;i<cells;i++){
			trueStatesCTM.put(i,0,rhoRight);
		}
		
		for(int i=0;i<5;i++){
			trueStatesCTM.put(i,0,rhoRight);
		}

		for(int i=cells-5;i<cells;i++){
				trueStatesCTM.put(i,0,0.35);
		}//to generate the shock propagating upstream
		
		trueStatesCTMPrior=trueStatesCTM.dup();
		
		measurementCTM=DoubleMatrix.zeros(((cells-1)/distMeasurement)+1,1);
		
		measureVarDensity = DoubleMatrix.eye(measurementCTM.getRows()).mul(0.0009);
		if (measurePerturb){
            for (int i=1;i<5;i++){
            	measureVarDensity.put(i*3,i*3,0.09);//assign low performance sensors
            }
		}
	    measureGeneratorDensity = new GaussianGenerator(measureVarDensity);
		
		numUp = 0;
	}

	public double[] getDensityBoundaries() {
		double[] res = new double[2];
		res[0] = trueStatesCTM.get(0,0);
		res[1] = trueStatesCTM.get(cells-1,0);
		return res;
	}

	public double flux(double density) { 	
		return density*speed(density);
	}

	public void newBoundaries(double _rhoLeft, double _rhoRight) {
		numUp = 0;
		rhoLeft = _rhoLeft;
		rhoRight = _rhoRight;
        trueStatesCTM=DoubleMatrix.zeros(cells,1);
		
		locDis = (int)cells/2; 
		
		for(int i=0;i<locDis;i++){
			trueStatesCTM.put(i,0,rhoLeft);
		}
		for(int i=locDis;i<cells;i++){
			trueStatesCTM.put(i,0,rhoRight);
		}
	}
		
	public void update() {		
		numUp++;
		propagateCTM();	
	}
	
	public void getMeasurement() {			
		getMeasurementCTM();	
	}

}
