package model;


import org.jblas.DoubleMatrix;
import org.jfree.data.xy.XYSeries;
import doubleMatrix.Concat;
import doubleMatrix.GaussianGenerator;
import trueSolution.TriangularTrue;
import trueSolution.TrueSolution;
import section.*;

@SuppressWarnings("unused")
public abstract class RoadModel {
	abstract public void initialThis();
	abstract public DoubleMatrix getDensityMean(DoubleMatrix mean);
	abstract public DoubleMatrix getDensityVar(DoubleMatrix var);
	abstract public DoubleMatrix propagate(DoubleMatrix _density);
	abstract public DoubleMatrix getMeasureVector();	
	
	//please refer to "RoadModel.java" in the project "DLKCF_Observable" for definitions of the following variables (if not detailed)
	public TrueSolution trueSolution;
	public int size;
	int cellsec;
	public int sizeMeasurements; 
	public int stepMeasurements;
	public DoubleMatrix modelVar;
	public DoubleMatrix measureVar;
	public DoubleMatrix measure;
	
	public Section section;	
	public String nameModel;
		
	public DoubleMatrix initialMean;
	public DoubleMatrix initialVar;	
	DoubleMatrix initialDensity;
	public GaussianGenerator initialGeneratorDensity;
	public DoubleMatrix modelVarDensity;//used to set modelVar&modelVar in Filter
	DoubleMatrix initialVarDensity;//used to set initialVar
	
	public void initial(TrueSolution _trueSolution, Section _section) {		
		trueSolution = _trueSolution;
		section = _section;
		cellsec = section.cellsec;	
		stepMeasurements = trueSolution.stepMeasurements;

		defineVarAndMean();
		initialThis(); 
		
 		nameModel = getClass().getSimpleName()+"_";
		nameModel+=Integer.toString((int) Math.floor(10*trueSolution.rhoLeft));
		nameModel+=Integer.toString((int) Math.floor(10*trueSolution.rhoRight));
		nameModel+="_";
		nameModel+=Integer.toString((int) Math.floor(section.index));
	}
		
	public double getSpaceStep() {
		return trueSolution.dx;
	}
	
	public static RoadModel createRoadModel(TrueSolution _trueSolution, Section _section, String nameEvolution ) {		
 		String name = "model."+nameEvolution;
 		try {			
			@SuppressWarnings("rawtypes")
			Class cl = Class.forName(name);			
			@SuppressWarnings("unchecked")
			RoadModel m = (RoadModel) cl.getConstructor().newInstance(); 
			m.initial(_trueSolution, _section);
			return m;
		}
		catch(Throwable t) {System.out.println("Road model creation has failed");t.printStackTrace();}
		return null;
	}
 	
 	public void defineVarAndMean() {
		double rhoLeft = section.densitysec1.get(0,0);
		double rhoRight = section.densitysec1.get(cellsec-1,0);		
		double wholeRhoLeft=trueSolution.totalrhoLeft;
		double wholeRhoRight=trueSolution.totalrhoRight;
		DoubleMatrix initialGuess=getInitial(wholeRhoLeft, wholeRhoRight);//linear interpolation
			
 		initialVarDensity = DoubleMatrix.eye(cellsec);
		double medium = 0.01*(rhoLeft + rhoRight);
		double a = 2/((double)cellsec);
		for (int i = 0; i<cellsec; i++) {
			initialVarDensity.put(i,i, 0.001*(section.index+1) + medium*Math.pow(1-Math.abs(i*a-1), 1));//covariance of the initial additive Gaussian noise
		}
		
		initialGeneratorDensity = new GaussianGenerator(initialVarDensity);

		modelVarDensity = DoubleMatrix.eye(cellsec).mmul(1);//just an initial set, will be updated in each step
		initialDensity = initialGuess.getRange(section.index*(cellsec-trueSolution.overlapsec),section.index*(cellsec-trueSolution.overlapsec)+cellsec ,0,1);	
		DoubleMatrix initialNoiseDensity = initialGeneratorDensity.sample();
		initialDensity=initialDensity.add(initialNoiseDensity);
 	}
 	
	public DoubleMatrix getInitial(double left, double right) {
		DoubleMatrix initial = new DoubleMatrix(trueSolution.cells);
		double slope = (right - left) / ((double)trueSolution.cells); 
		for(int i=0; i<trueSolution.cells; i++) {
			initial.put(i,left + i*slope);
		}
		return initial;
	}
	
	public DoubleMatrix getRiemman(double left, double right) {
		DoubleMatrix initial = new DoubleMatrix(cellsec);
		for(int i=0; i<((int) cellsec/2); i++) {
			initial.put(i,left);
		}
		for(int i=((int) cellsec/2); i<cellsec; i++) {
			initial.put(i, right);
		}
		return initial;
	}
	
	public void newBoundaries(double rhoLeft, double rhoRight) {
		trueSolution.newBoundaries(rhoLeft, rhoRight);
		section=trueSolution.setSections();
		initial(trueSolution,section);
	}

	public void updateTrueSolution() {
		trueSolution.update();
	}
	
	public void setSetion(Section _section){
		section=_section;
		
	}
	
	public void updateModelVar() {
		if (trueSolution.individual){
			if ((section.index)%2==0){
				modelVar=DoubleMatrix.eye(cellsec).mmul(0.0025);
			}
			else{
				modelVar=DoubleMatrix.eye(cellsec).mmul(0.0025);
			}
		}
		else{
			if ((section.index)%2==0){
				modelVar=DoubleMatrix.eye(cellsec).mmul(0.0025);
			}
			else{
				modelVar=DoubleMatrix.eye(cellsec).mmul(0.0025);
			}
		}
	}
	
	public void updateMeasurement(){
        measureVar=trueSolution.measureVarDensity;		
		sizeMeasurements = measureVar.getColumns();		
		measure = DoubleMatrix.zeros(sizeMeasurements,size);
		if (sizeMeasurements==2){
			measure.put(0,0,1);
			measure.put(1,cellsec-1,1);
		}
		else{
			if (trueSolution.overlapsec!=1){
				for (int i=0;i<sizeMeasurements;i++){				
					if(i<2){
						measure.put(i, i*(trueSolution.overlapsec-1),1);
					}
					else{
						measure.put(i, cellsec-1-(3-i)*(trueSolution.overlapsec-1),1);
					}												
				}	
			}
			else{
				int distMeasurements=(cellsec-1)/(sizeMeasurements-1);
				for (int i=0;i<sizeMeasurements;i++){
				measure.put(i,i*distMeasurements,1);
				}
			}
		}			
	}	
}
