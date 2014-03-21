package trueSolutionCTM;

import org.jblas.DoubleMatrix;

public class TriangularTrueCTM extends TrueSolutionCTM{
	
	public TriangularTrueCTM(double _rhoLeft, double _rhoRight, double _dt, double _dx, int _cells, int _distMeasurement,int _numSections, boolean _measurePerturb){
		rhoCritical = 0.225; 
		rhoMax = 1; 
		speedMax = 1; 
		initial(_rhoLeft, _rhoRight, _dt, _dx, _cells,_distMeasurement, _numSections,_measurePerturb);
	}
	
	public double speed(double density) {
		if (density<=rhoCritical) return speedMax;
		else return rhoCritical*speedMax*(rhoMax-density)/(density*(rhoMax-rhoCritical));
	}

	@Override
	public void propagateCTM() {
		// TODO Auto-generated method stub
		DoubleMatrix _densitynext=DoubleMatrix.zeros(trueStatesCTM.getRows(), 1);
		trueStatesCTMPrior=trueStatesCTM;
		_densitynext.put(trueStatesCTM.getRows()-1,0,trueStatesCTM.get(trueStatesCTM.getRows()-1,0));
		_densitynext.put(0,0,trueStatesCTM.get(0,0));
		for (int i=1;i<trueStatesCTM.getRows()-1;i++){
			_densitynext.put(i,0,trueStatesCTM.get(i,0)+(dt/dx)*(computeFlux(trueStatesCTM.get(i-1,0),trueStatesCTM.get(i,0))-computeFlux(trueStatesCTM.get(i,0),trueStatesCTM.get(i+1,0))));
		}
		trueStatesCTM=_densitynext;
	}
	
	public void getMeasurementCTM() {
		DoubleMatrix noiseDensity = measureGeneratorDensity.sample();
		for (int i=0;i<measurementCTM.getRows();i++){
			double pointDensity=trueStatesCTM.get(i*distMeasurement);
			pointDensity+=noiseDensity.get(i);
			if (pointDensity<0){
				pointDensity=0;
			}
			else if (pointDensity>1){
				pointDensity=1;
			}
			measurementCTM.put(i,0,pointDensity);
		}
	}
	
	public double sending(double density) {
		// TODO Auto-generated method stub
		if (density>=rhoCritical){
			return rhoCritical*speedMax;
		}
		else return density*speedMax;
	}

	@Override
	public double receiving(double density) {
		// TODO Auto-generated method stub
		double w=(rhoCritical*speedMax)/(rhoMax-rhoCritical);
		if (density<=rhoCritical){
			return rhoCritical*speedMax;
		}
		else return w*(rhoMax-density);
	}
	
	@Override
	public double computeFlux(double density1, double density2) {
	    if (sending(density1)<=receiving(density2)){
	    	return sending(density1);
	    }
	    else return receiving(density2);
	}

}
