package trueSolution;

import org.jblas.DoubleMatrix;
import org.jfree.data.xy.XYSeries;

import doubleMatrix.GaussianGenerator;
import filters.Filter;
import section.*;
import model.*;
import trueSolutionCTM.*;


//This TrueSolution is used to obtain the true state in a section from the true state of the entire network 
public abstract class TrueSolution {

//please refer to "TrueSolution.java" in the project "DLKCF_Observable" for definitions of the following variables (if not detailed)
	int locDis;	
	public double rhoLeft; 
	public double rhoRight; 
	public double totalrhoLeft; 
	public double totalrhoRight; 
	protected double fluxPrimeLeft;
	protected double fluxPrimeRight;
	
	public double rhoMax;
	public double speedMax;
	public double rhoCritical;
	
	public double dt;
	public double dx;
	public int cells;
	public int cellsec;
	public int overlapsec;
	public int index;	
	public int numSections;
	
	public Section section;
	public TrueSolutionCTM trueCTM;
	

	public int perturb;//type of disagreement on the system dynamics used by different agents in estimation
	public boolean measurePerturb;//if low quality sensors exist	
    public boolean individual;//if inter-agent communication exists 
	
	public DoubleMatrix measureVarDensity;//the measurement error covariance matrix
	public DoubleMatrix measurementsDensity;//the measurement vector
	
	public DoubleMatrix trueStates;
	public DoubleMatrix trueStatesPrior;
	
	public int stepMeasurements=1;
	public int sizeMeasurements; 
	
	public int numUp;
	
	public TrueSolution trueLeft;
	public TrueSolution trueRight;

	abstract public double speed(double density);
	abstract public void propagate();	
	abstract public double sending(double density);
	abstract public double receiving(double density);
	abstract public double computeFlux(double density1, double density2);

	
	public void initial(TrueSolutionCTM _trueCTM, double _rhoLeft, double _rhoRight, double _totalrhoLeft, double _totalrhoRight, double _dt, double _dx, int _cellsec, int _overlapsec, int _index,int _numSecs,int _perturb, boolean _measurePerturb, boolean _individual) {
		trueCTM=_trueCTM;
		index=_index;
		rhoLeft = _rhoLeft;
		rhoRight = _rhoRight;
		totalrhoLeft=_totalrhoLeft;
		totalrhoRight=_totalrhoRight;
		dx = _dx;
		dt = _dt;
		perturb=_perturb;
		measurePerturb=_measurePerturb;
		individual=_individual;

		cellsec=_cellsec;
		overlapsec=_overlapsec;
		numSections=_numSecs;
		cells=(numSections-1)*(cellsec-overlapsec)+cellsec;
				
		sizeMeasurements=4;
		if (overlapsec==1){
			measureVarDensity = trueCTM.measureVarDensity.getRange(index*3, index*3+4, index*3, index*3+4).dup();
		}
		else{
			measureVarDensity = trueCTM.measureVarDensity.getRange(index*2, index*2+4, index*2, index*2+4).dup();
		}

		trueStates=DoubleMatrix.zeros(cellsec,1);
		
		locDis = (int)cellsec/2; //localization of the discontinuity
		
		trueStates=trueCTM.trueStatesCTM.getRange(index*(cellsec-overlapsec), index*(cellsec-overlapsec)+cellsec, 0, 1);//get the true state in the section from the entire network
	    trueStatesPrior=trueStates.dup();
		
		numUp = 0;
	}

	public double[] getDensityBoundaries() {
		double[] res = new double[2];
		res[0] = trueStates.get(0,0);
		res[1] = trueStates.get(cellsec-1,0);
		return res;
	}

    public double flux(double density) { 	
		return density*speed(density);
	}

	public void newBoundaries(double _rhoLeft, double _rhoRight) {
		numUp = 0;
		rhoLeft = _rhoLeft;
		rhoRight = _rhoRight;
        trueStates=DoubleMatrix.zeros(cellsec,1);
		
		locDis = (int)cells/2; 
		
		for(int i=0;i<locDis;i++){
			trueStates.put(i,0,rhoLeft);
		}
		for(int i=locDis;i<cellsec;i++){
			trueStates.put(i,0,rhoRight);
		}
	}
		
	
	public void update() {		
		numUp++;
		propagate();	
	}
	
    public void updatemeasurement() {		
		newMeasurements();		
	}
    	
	private void newMeasurements() {
		if (numUp - stepMeasurements*(numUp/stepMeasurements) == 0) getMeasurementsAll();
	}
	
	private void getMeasurementsAll() {
		measurementsDensity = DoubleMatrix.zeros(sizeMeasurements, 1);			
		if (overlapsec==1){	
			measurementsDensity=trueCTM.measurementCTM.getRange(index*3, index*3+4, 0, 1);
		}
		else{
			measurementsDensity=trueCTM.measurementCTM.getRange(index*2, index*2+4, 0, 1);
		}
			
	}

 	@SuppressWarnings({ "rawtypes", "unchecked" })

 	public Section setSections() {
 		double[] paramsEven=new double [4];
 		double[] paramsOdd=new double [4];
 		
 		//use different parameters in the system dynamics of the SMM
 		if (perturb==0){
 	 		paramsEven[0]=1;//maximum density
 	 		paramsEven[1]=0.225;//critical density
 	 		paramsEven[2]=1;//maximum speed
 	 		paramsEven[3]=-1;
 	 		paramsOdd=paramsEven;
 		}
 		else if (perturb==1){
 	 		paramsOdd[0]=1.1;
 	 		paramsOdd[1]=0.3;
 	 		paramsOdd[2]=0.9;
 	 		paramsOdd[3]=-1;
 	 		paramsEven[0]=0.9;
 	 		paramsEven[1]=0.2;
 	 		paramsEven[2]=1.2;
 	 		paramsEven[3]=-1;	 		
 		}
 		else if (perturb==2){
 	 		paramsEven[0]=1.1;
 	 		paramsEven[1]=0.3;
 	 		paramsEven[2]=0.9;
 	 		paramsEven[3]=-1;
 	 		paramsOdd=paramsEven;
 		}
 		else{
 			paramsOdd[0]=0.9;
 	 		paramsOdd[1]=0.2;
 	 		paramsOdd[2]=1.2;
 	 		paramsOdd[3]=-1;
 	 		paramsEven=paramsOdd;
 		}				
 		Section secs ;
 		
 		DoubleMatrix _measurement=DoubleMatrix.zeros(sizeMeasurements, 0);
 		//identify mode based on the measurement
 		if (overlapsec==1){
			 _measurement=trueCTM.measurementCTM.getRange(index*3, index*3+4, 0, 1);
		}
		else{
			 _measurement=trueCTM.measurementCTM.getRange(index*2, index*2+4, 0, 1);
		}			
 			if ((index)%2==0){			
 				if (_measurement.get(0,0)<=paramsEven[1] && _measurement.get(3,0)<=paramsEven[1]){
 	 				    secs=Section.createSection("FF",this);
 	 				}
 	 				else if (_measurement.get(0,0)>paramsEven[1] && _measurement.get(3,0)>paramsEven[1]){
 	 				    secs=Section.createSection("CC",this);
 	 				}
 	 				else if (_measurement.get(0,0)>paramsEven[1] && _measurement.get(3,0)<=paramsEven[1]){
 	 				    secs=Section.createSection("CF",this);	 					
 	 			}
 	 				else {
 	 				    secs=Section.createSection("FC",this);
 	 				}
 			}
 			else{
 				if (_measurement.get(0,0)<=paramsOdd[1] && _measurement.get(3,0)<=paramsOdd[1]){
 	 				    secs=Section.createSection("FF",this);
 	 				}
 	 				else if (_measurement.get(0,0)>paramsOdd[1] && _measurement.get(3,0)>paramsOdd[1]){
 	 				    secs=Section.createSection("CC",this);
 	 				}
 	 				else if (_measurement.get(0,0)>paramsOdd[1] && _measurement.get(3,0)<=paramsOdd[1]){
 	 				    secs=Section.createSection("CF",this); 	 					
 	 			}
 	 				else {
 	 				    secs=Section.createSection("FC",this);
 	 				}
 			}
		    secs.densitysec1=trueStates;
		    secs.index=index;
		    //different agents use different model parameters
		    if ((secs.index)%2==0){
		    	secs.setNewParameters(paramsEven);
		    }
		    else{
		    	secs.setNewParameters(paramsOdd);
		    }
			    
		    if (secs.getClass().getSimpleName().equals("CF")){
		    	secs.getwavefront();
 		    }
		    else if (secs.getClass().getSimpleName().equals("FC")){
		    	secs.getwavefront();
 		    }		    
		    if (secs.getClass().getSimpleName().equals("FF")||secs.getClass().getSimpleName().equals("CC")){
		    	secs.getModelA();
	 		    secs.getModelB1();
			    secs.getModelB2();
		    }
 		return secs;
 	}
 		
 	public RoadModel setRoadModels(){
 		Section secs = setSections();
 		RoadModel rms;	
 		rms=RoadModel.createRoadModel(this, secs, "D");	
 		return rms;
 	}
}
