package trueSolution;

import org.jblas.DoubleMatrix;

import trueSolutionCTM.TrueSolutionCTM;

public class TriangularTrue extends TrueSolution{
	
	public TriangularTrue(TrueSolutionCTM _trueCTM,double _rhoLeft, double _rhoRight,  double _totalrhoLeft, double _totalrhoRight, double _dt, double _dx, int _cellsec, int _overlapsec, int _index, int _numSecs,int _perturb, boolean _measurePerturb, boolean _individual){
		rhoCritical = 0.225; 
		rhoMax = 1; 
		speedMax = 1; 
		initial(_trueCTM,_rhoLeft, _rhoRight, _totalrhoLeft, _totalrhoRight,_dt, _dx, _cellsec, _overlapsec, _index,_numSecs, _perturb, _measurePerturb,_individual);
	}
	public double speed(double density) {
		if (density<=rhoCritical) return speedMax;
		else return rhoCritical*speedMax*(rhoMax-density)/(density*(rhoMax-rhoCritical));
	}
	@Override
	public void propagate() {
		// TODO Auto-generated method stub
	    trueStatesPrior=trueStates;
	    trueStates=trueCTM.trueStatesCTM.getRange(index*(cellsec-overlapsec), index*(cellsec-overlapsec)+cellsec, 0, 1);

		
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
