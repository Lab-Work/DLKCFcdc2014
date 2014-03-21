package filters;
import filters.KCF;
import section.Section;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import model.RoadModel;
import org.jblas.DoubleMatrix;
import org.jfree.data.xy.XYSeries;
import Jama.EigenvalueDecomposition;
import Jama.Matrix;
import doubleMatrix.Concat;
import doubleMatrix.GaussianGenerator;
import doubleMatrix.InverseMatrix;
import trueSolution.TrueSolution;
import trueSolutionCTM.TrueSolutionCTM;


public class Estimation {

	RoadModel roadModel;	
	Filter filter;
	boolean consensus;

	public Estimation(RoadModel _roadModel, String _nameFilter, boolean _consensus) {
		roadModel = _roadModel;		
		filter = Filter.createFilter(this);			
		consensus=_consensus;
	}
	
	public void nextStepWithoutUpdateTrue() {
		filter.nextStep();
	}
	public void nextStepNoDataWithoutUpdateTrue() {
		filter.nextStepNoData();
	}
		
	public void updateTrueSolution() {
		roadModel.updateTrueSolution();//
	}

	public DoubleMatrix propagate(DoubleMatrix _density) {		
      return roadModel.propagate(_density);
	}
	
	public DoubleMatrix getDensityMean() {
		return roadModel.getDensityMean(filter.mean);
	}
	
	public void setSection(Section _section){		
			roadModel.setSetion(_section);
			roadModel.updateModelVar();
	}
	
	
	

	
	public static void exportResult(TrueSolutionCTM trueSolutionCTM, TrueSolution[] trueSolution,TrueSolution[] trueSolution1,int limite, int bdry, String folder) {
		try {					
			
			int numSections = trueSolution[0].numSections;
			int numSections1 = trueSolution1[0].numSections;

			//please refer to "Estimation.java" in the project "DLKCF_Observable" for definitions of the following variables (if not detailed)
			//another estimation (with suffix "1" added to all its variable names) introduced to compare estimates with and without consensus 
			BufferedWriter[] writerSection = new BufferedWriter[numSections];
			BufferedWriter[] writerSection1 = new BufferedWriter[numSections1];

			BufferedWriter writerTrue;
			BufferedWriter writerAvr;
			BufferedWriter writerAvr1;

			BufferedWriter Gamma;
			BufferedWriter Error;
			BufferedWriter Error1;
			
			BufferedWriter[] Mode= new BufferedWriter[numSections];
			BufferedWriter[] Mode1= new BufferedWriter[numSections1];

			BufferedWriter writerDisagreement;
			BufferedWriter writerDisagreement1;
			

				System.out.print("Beginning of simulation ");
				
				new File("results/"+folder).mkdir();
			    
				RoadModel[] roadModels=new RoadModel[numSections];
				RoadModel[] roadModels1=new RoadModel[numSections1];			
				for (int i=0;i<numSections;i++){					
					trueSolution[i].section=trueSolution[i].setSections();					
				}
				for (int i=0;i<numSections1;i++){					
					trueSolution1[i].section=trueSolution1[i].setSections();					
				}		
				for (int i=0;i<numSections;i++){					
					roadModels[i]=trueSolution[i].setRoadModels();
				}			
				for (int i=0;i<numSections1;i++){					
					roadModels1[i]=trueSolution1[i].setRoadModels();
				}
						
				writerTrue = new BufferedWriter(new FileWriter(new File("results/"+folder+"/"+"trueState.csv")));
				writerAvr = new BufferedWriter(new FileWriter(new File("results/"+folder+"/"+"writerAvr.csv")));
				writerAvr1 = new BufferedWriter(new FileWriter(new File("results/"+folder+"/"+"writerAvr1.csv")));				
				Error = new BufferedWriter(new FileWriter(new File("results/"+folder+"/"+"Error.csv")));
				Error1 = new BufferedWriter(new FileWriter(new File("results/"+folder+"/"+"Error1.csv")));
				Gamma = new BufferedWriter(new FileWriter(new File("results/"+folder+"/"+"Gamma.csv")));
				writerDisagreement = new BufferedWriter(new FileWriter(new File("results/"+folder+"/"+"Disagreement.csv")));
				writerDisagreement1 = new BufferedWriter(new FileWriter(new File("results/"+folder+"/"+"Disagreement1.csv")));

				Estimation[] estimations = new Estimation[numSections];
				Estimation[] estimations1 = new Estimation[numSections1];

				for (int i = 0; i<numSections; i++) {
					estimations[i] = new Estimation(roadModels[i], "KCF",true);		
					String S = "results/"+folder+"/"+roadModels[i].nameModel;
					writerSection[i] = new BufferedWriter(new FileWriter(new File(S+".csv")));			
					Mode[i]=new BufferedWriter(new FileWriter(new File(S+"Mode.csv")));
				} 
				
				for (int i = 0; i<numSections1; i++) {
					estimations1[i] = new Estimation(roadModels1[i], "KCF",false);
					String S = "results/"+folder+"/"+roadModels1[i].nameModel;	
					writerSection1[i] = new BufferedWriter(new FileWriter(new File(S+"_1.csv")));
					Mode1[i]=new BufferedWriter(new FileWriter(new File(S+"Mode1.csv")));
				}

				int overlap=roadModels[0].trueSolution.overlapsec;
				int cellsec=roadModels[0].section.cellsec;
				int overlap1=roadModels1[0].trueSolution.overlapsec;
				int cellsec1=roadModels1[0].section.cellsec;

				DoubleMatrix I1=DoubleMatrix.zeros(roadModels[0].section.cellsec,overlap);
				for (int i=0;i<overlap;i++){
					I1.put(i,i,1);
				}
				DoubleMatrix I2=DoubleMatrix.zeros(roadModels[0].section.cellsec,overlap);
				for (int i=roadModels[0].section.cellsec-overlap;i<roadModels[0].section.cellsec;i++){
					I2.put(i,i-(roadModels[0].section.cellsec-overlap),1);
				}
				DoubleMatrix I123=DoubleMatrix.concatVertically(I1.transpose(), I2.transpose());
				DoubleMatrix H_hat=I2.transpose();
				for (int i=0;i<numSections-2;i++){
					H_hat=Concat.Diagonal(H_hat, I123);
				}
				H_hat=Concat.Diagonal(H_hat, I1.transpose());
                
				DoubleMatrix L_tilde=DoubleMatrix.zeros((numSections-1)*2*overlap, numSections*cellsec);
				for(int i=0;i<numSections-1;i++){
					for(int j=0;j<overlap;j++){
						L_tilde.put(i*(2*overlap)+j,(i+1)*(cellsec)+j, 1);
					}
					for(int j=0;j<overlap;j++){
						L_tilde.put(i*(2*overlap)+overlap+j,(i+1)*(cellsec)-overlap+j, 1);
					}
					for(int j=0;j<overlap;j++){
						L_tilde.put(i*(2*overlap)+j,(i+1)*(cellsec)-overlap+j, -1);
					}
					for(int j=0;j<overlap;j++){
						L_tilde.put(i*(2*overlap)+overlap+j,(i+1)*(cellsec)+j, -1);
					}
				}
				
				DoubleMatrix[] S=new DoubleMatrix[numSections];
				
				trueSolutionCTM.getMeasurement();
			
				for (int k=0; k<limite; k++) {
					for (int i=0;i<trueSolutionCTM.trueStatesCTM.getRows();i++){
						writerTrue.write(trueSolutionCTM.trueStatesCTM.get(i,0)+",");
					}
					writerTrue.write("\n");
					
					DoubleMatrix[] mean = new DoubleMatrix[numSections];
					DoubleMatrix[] mean_1 = new DoubleMatrix[numSections1];
					
					double ErrorAvg=0;
					double ErrorAvg1=0;
					
					for (int i=0; i<numSections; i++) {
						mean [i] = estimations[i].getDensityMean();
						if (numSections==1){
							for (int j=0; j<mean[i].length ; j++) {
								writerAvr.write(mean[i].get(j)+",");							
							}
						}
						else{
							if (i==0){
								for (int j=0; j<mean[i].length-overlap ; j++) {
									writerAvr.write(mean[i].get(j)+",");							
								}
							}
							else if (i>0 && i<numSections-1){
								for (int j=0; j<overlap ; j++) {
									writerAvr.write((mean[i].get(j)+mean[i-1].get(cellsec-overlap+j))/2+",");	
								}
								for (int j=0; j<mean[i].length-2*overlap ; j++) {
									writerAvr.write(mean[i].get(overlap+j)+",");								
								}
							}
							else if (i==numSections-1){
								for (int j=0; j<overlap ; j++) {
									writerAvr.write((mean[i].get(j)+mean[i-1].get(cellsec-overlap+j))/2+",");	
								}
								for (int j=0; j<mean[i].length-overlap ; j++) {
									writerAvr.write(mean[i].get(overlap+j)+",");								
								}
							}
						}
						for (int j=0; j<mean[i].length ; j++) {
							writerSection[i].write(mean[i].get(j)+",");							
						}
						double errorSection=(((mean[i].sub(trueSolution[i].trueStates)).transpose().mmul((mean[i].sub(trueSolution[i].trueStates)))).get(0,0))/((double)cellsec);
						Error.write(errorSection+",");
						ErrorAvg=ErrorAvg+errorSection/((double)numSections);
						writerSection[i].write("\n");
					}
					
					Error.write(ErrorAvg+",");
					
					for (int i=0; i<numSections1; i++) {						
						mean_1 [i] = estimations1[i].getDensityMean();
						if (numSections1==1){
							for (int j=0; j<mean_1[i].length ; j++) {
								writerAvr1.write(mean_1[i].get(j)+",");							
							}
						}
						else{
							if (i==0){
								for (int j=0; j<mean_1[i].length-overlap1 ; j++) {
									writerAvr1.write(mean_1[i].get(j)+",");							
								}
							}
							else if (i>0 && i<numSections1-1){
								for (int j=0; j<overlap1 ; j++) {
									writerAvr1.write((mean_1[i].get(j)+mean_1[i-1].get(cellsec1-overlap1+j))/2+",");	
								}
								for (int j=0; j<mean_1[i].length-2*overlap1 ; j++) {
									writerAvr1.write(mean_1[i].get(overlap1+j)+",");								
								}
							}
							else if (i==numSections1-1){
								for (int j=0; j<overlap1 ; j++) {
									writerAvr1.write((mean_1[i].get(j)+mean_1[i-1].get(cellsec1-overlap1+j))/2+",");	
								}
								for (int j=0; j<mean_1[i].length-overlap1 ; j++) {
									writerAvr1.write(mean_1[i].get(overlap1+j)+",");								
								}
							}
						}										
						for (int j=0; j<mean_1[i].length ; j++) {
							writerSection1[i].write(mean_1[i].get(j)+",");						
						}						
						double errorSection1=(((mean_1[i].sub(trueSolution1[i].trueStates).transpose()).mmul((mean_1[i].sub(trueSolution1[i].trueStates)))).get(0,0))/((double)cellsec1);
						Error1.write(errorSection1+",");					
						ErrorAvg1=ErrorAvg1+errorSection1/((double)numSections1);
						writerSection1[i].write("\n");					
					}
					Error1.write(ErrorAvg1+",");
					writerAvr.write("\n");
					writerAvr1.write("\n");
					Error.write("\n");
					Error1.write("\n");					
					
					Section[] secs=new Section[numSections];
					Section[] secs1=new Section[numSections1];					
					for(int i=0;i<numSections;i++){
						secs[i]=trueSolution[i].setSections();//every step, update section 
					}
					for(int i=0;i<numSections1;i++){
						secs1[i]=trueSolution1[i].setSections();//every step, update section 
					}
					for(int i=0;i<numSections;i++){
						secs[i].Estimates=mean[i];
					}
					for(int i=0;i<numSections1;i++){
						secs1[i].Estimates=mean_1[i];
					}
					for(int i=0;i<numSections;i++){
						if(secs[i].getClass().getSimpleName().equals("FC")){
							secs[i].getwavefrontfromEstimation();
							if(secs[i].wavefront<0){
								secs[i].wavefront=0;
							}
							if(secs[i].wavefront<=cellsec-2){
							    if(trueSolution[i].flux(secs[i].Estimates.get(secs[i].wavefront,0))-trueSolution[i].flux(secs[i].Estimates.get(secs[i].wavefront+1,0))>0){			
								    secs[i].wavedirection=0;
							    }
						        else{
								    secs[i].wavedirection=1;
						        }
							}
							else if (secs[i].wavefront==cellsec-1){
								if(secs[i].index!=numSections-1){
									if(trueSolution[i].flux(secs[i].Estimates.get(secs[i].wavefront,0))-trueSolution[i+1].flux(secs[i].Estimates.get(overlap,0))>0){										
									    secs[i].wavedirection=0;									
								    }
								    else{					
										secs[i].wavedirection=1;									
								    }
								}
								else{
									secs[i].wavedirection=0;
								}								
							}
							secs[i].getModelA();
				 		    secs[i].getModelB1();
						    secs[i].getModelB2();						   						
						}
						else if (secs[i].getClass().getSimpleName().equals("CF")){
							secs[i].getwavefrontfromEstimation();
							secs[i].getModelA();
				 		    secs[i].getModelB1();
						    secs[i].getModelB2();
						}		
					}
					
					for(int i=0;i<numSections1;i++){
						if(secs1[i].getClass().getSimpleName().equals("FC")){
							secs1[i].getwavefrontfromEstimation();
							if(secs1[i].wavefront<0){
								secs1[i].wavefront=0;
							}
							if(secs1[i].wavefront<=cellsec1-2){
							    if(trueSolution1[i].flux(secs1[i].Estimates.get(secs1[i].wavefront,0))-trueSolution1[i].flux(secs1[i].Estimates.get(secs1[i].wavefront+1,0))>0){
							        secs1[i].wavedirection=0;		
						         }
						        else{
								    secs1[i].wavedirection=1;
						        }
							}
							else if (secs1[i].wavefront==cellsec1-1){
							    if(secs1[i].index!=numSections1-1){
								    if(trueSolution1[i].flux(secs1[i].Estimates.get(secs1[i].wavefront,0))-trueSolution1[i+1].flux(secs1[i].Estimates.get(overlap1,0))>0){										
										secs1[i].wavedirection=0;									
								    }
								    else{					
										secs1[i].wavedirection=1;									
								    }
								}
								else{
									secs1[i].wavedirection=0;
								}
							}
							secs1[i].getModelA();
				 		    secs1[i].getModelB1();
						    secs1[i].getModelB2();   						
						}
						else if (secs1[i].getClass().getSimpleName().equals("CF")){
							secs1[i].getwavefrontfromEstimation();
							secs1[i].getModelA();
				 		    secs1[i].getModelB1();
						    secs1[i].getModelB2();
						}		
					}					
					for (int i=0; i<numSections; i++) {
						estimations[i].setSection(secs[i]);
					}
					for (int i=0; i<numSections1; i++) {
						estimations1[i].setSection(secs1[i]);
					}
					for (int i=0; i<numSections; i++) {
						trueSolution[i].section=secs[i];
					}					
					for (int i=0; i<numSections1; i++) {
						trueSolution1[i].section=secs1[i];
					}
					
					for (int i=0; i<numSections; i++) {
						Mode[i].write(estimations[i].roadModel.section.getClass().getSimpleName()+",");
						Mode[i].write(estimations[i].roadModel.section.wavefront+",");
						Mode[i].write(estimations[i].roadModel.section._wavefront+",");
						Mode[i].write(estimations[i].roadModel.section._wavedirection+",");
						Mode[i].write(estimations[i].roadModel.section.wavedirection+",");
						Mode[i].write("\n");
					}
					for (int i=0; i<numSections1; i++) {
						Mode1[i].write(estimations1[i].roadModel.section.getClass().getSimpleName()+",");
						Mode1[i].write(estimations1[i].roadModel.section.wavefront+",");
						Mode1[i].write(estimations1[i].roadModel.section._wavefront+",");
						Mode1[i].write(estimations1[i].roadModel.section._wavedirection+",");
						Mode1[i].write(estimations1[i].roadModel.section.wavedirection+",");
						Mode1[i].write("\n");
					}

					
					for (int i=0; i<numSections; i++) {
						estimations[i].filter.getNewParametersFromModel();
					}
					
					for (int i=0; i<numSections1; i++) {
						estimations1[i].filter.getNewParametersFromModel();
					}					
									
					double disag=0;
					for (int i=0; i<numSections-1; i++) {
						disag= disag+(((mean[i].getRange(cellsec-overlap, cellsec,0,1).sub(mean[i+1].getRange(0, overlap,0,1))).transpose().mmul(mean[i].getRange(cellsec-overlap, cellsec,0,1).sub(mean[i+1].getRange(0, overlap,0,1)))).get(0,0))/((double)overlap);
					}
					disag=disag/((double)(numSections-1));
					writerDisagreement.write(disag+",");
					for (int i=0; i<numSections-1; i++) {
						writerDisagreement.write(((mean[i].getRange(cellsec-overlap, cellsec,0,1).sub(mean[i+1].getRange(0, overlap,0,1))).transpose().mmul(mean[i].getRange(cellsec-overlap, cellsec,0,1).sub(mean[i+1].getRange(0, overlap,0,1)))).get(0,0)+",");
					}
					writerDisagreement.write("\n");
					
					double disag1=0;
					for (int i=0; i<numSections1-1; i++) {
						disag1= disag1+(((mean_1[i].getRange(cellsec1-overlap1, cellsec1,0,1).sub(mean_1[i+1].getRange(0, overlap1,0,1))).transpose().mmul(mean_1[i].getRange(cellsec1-overlap1, cellsec1,0,1).sub(mean_1[i+1].getRange(0, overlap1,0,1)))).get(0,0))/((double)overlap1);
					}
					disag1=disag1/((double)(numSections-1));
					writerDisagreement1.write(disag1+",");
					for (int i=0; i<numSections1-1; i++) {
						writerDisagreement1.write(((mean_1[i].getRange(cellsec1-overlap1, cellsec1,0,1).sub(mean_1[i+1].getRange(0, overlap1,0,1))).transpose().mmul(mean_1[i].getRange(cellsec1-overlap1, cellsec1,0,1).sub(mean_1[i+1].getRange(0, overlap1,0,1)))).get(0,0)+",");
					}
					writerDisagreement1.write("\n");
										
					trueSolutionCTM.update();

					//bdry=0: assuming the cell density to the left (right) side of the upstream (downstream) cell is the same as the density of the upstream (downstream) cell
					//bdry=1: sinusoidal boundary condition with relatively large frequency
					//bdry=2: one-sided sinusoidal boundary condition with relatively low frequency, used in the paper
					int cells=trueSolutionCTM.trueStatesCTMPrior.getRows();
					if (bdry==0){
						double inflow=trueSolutionCTM.receiving(trueSolutionCTM.trueStatesCTMPrior.get(0,0));
						double outflow=trueSolutionCTM.sending(trueSolutionCTM.trueStatesCTMPrior.get(cells-1,0));
						double upOut=trueSolutionCTM.computeFlux(trueSolutionCTM.trueStatesCTMPrior.get(0,0),trueSolutionCTM.trueStatesCTMPrior.get(1,0));
						trueSolutionCTM.trueStatesCTM.put(0,0,trueSolutionCTM.trueStatesCTMPrior.get(0,0)+(trueSolutionCTM.dt/trueSolutionCTM.dx)*(inflow-upOut));

						double downIn=trueSolutionCTM.computeFlux(trueSolutionCTM.trueStatesCTMPrior.get(cells-2,0),trueSolutionCTM.trueStatesCTMPrior.get(cells-1,0));
						trueSolutionCTM.trueStatesCTM.put(cells-1,0,trueSolutionCTM.trueStatesCTMPrior.get(cells-1,0)+(trueSolutionCTM.dt/trueSolutionCTM.dx)*(downIn-outflow));
					}
					else if (bdry==1){
//						double inflow=0.1125+0.1125*Math.sin(k*(Math.PI/1000)+(Math.PI));
//						if (inflow>trueSolutionCTM.receiving(trueSolutionCTM.trueStatesCTMPrior.get(0,0))){
//							inflow=trueSolutionCTM.receiving(trueSolutionCTM.trueStatesCTMPrior.get(0,0));
//						}
//						
//						double outflow=trueSolutionCTM.sending(trueSolutionCTM.trueStatesCTMPrior.get(cells-1,0));
//						double outflow=0.1125+0.1125*Math.sin(k*(Math.PI/1000));
//						if (outflow>trueSolution[numSections-1].sending(trueSolution[numSections-1].trueStatesPrior.get(cellsec-1,0))){
//							outflow=trueSolution[numSections-1].sending(trueSolution[numSections-1].trueStatesPrior.get(cellsec-1,0));
//						}
						
//						double upOut=trueSolutionCTM.computeFlux(trueSolutionCTM.trueStatesCTMPrior.get(0,0),trueSolutionCTM.trueStatesCTMPrior.get(1,0));
//						trueSolutionCTM.trueStatesCTM.put(0,0,trueSolutionCTM.trueStatesCTMPrior.get(0,0)+(trueSolutionCTM.dt/trueSolutionCTM.dx)*(inflow-upOut));

//						double downIn=trueSolutionCTM.computeFlux(trueSolutionCTM.trueStatesCTMPrior.get(cells-2,0),trueSolutionCTM.trueStatesCTMPrior.get(cells-1,0));
//						trueSolutionCTM.trueStatesCTM.put(cells-1,0,trueSolutionCTM.trueStatesCTMPrior.get(cells-1,0)+(trueSolutionCTM.dt/trueSolutionCTM.dx)*(downIn-outflow));
					}
					else if (bdry==2){
						double inflow=0.1125+0.1125*Math.sin(k*(Math.PI/4000)+(Math.PI));
						if (inflow>trueSolutionCTM.receiving(trueSolutionCTM.trueStatesCTMPrior.get(0,0))){
							inflow=trueSolutionCTM.receiving(trueSolutionCTM.trueStatesCTMPrior.get(0,0));
						}
						
						double outflow=trueSolutionCTM.receiving(trueSolutionCTM.trueStatesCTMPrior.get(cells-1,0));
						
						double upOut=trueSolutionCTM.computeFlux(trueSolutionCTM.trueStatesCTMPrior.get(0,0),trueSolutionCTM.trueStatesCTMPrior.get(1,0));
						trueSolutionCTM.trueStatesCTM.put(0,0,trueSolutionCTM.trueStatesCTMPrior.get(0,0)+(trueSolutionCTM.dt/trueSolutionCTM.dx)*(inflow-upOut));

						double downIn=trueSolutionCTM.computeFlux(trueSolutionCTM.trueStatesCTMPrior.get(cells-2,0),trueSolutionCTM.trueStatesCTMPrior.get(cells-1,0));
						trueSolutionCTM.trueStatesCTM.put(cells-1,0,trueSolutionCTM.trueStatesCTMPrior.get(cells-1,0)+(trueSolutionCTM.dt/trueSolutionCTM.dx)*(downIn-outflow));
					}
					else{
						
					}
										
					trueSolutionCTM.getMeasurement();
					
					for(int i=0; i<numSections; i++){
						trueSolution[i].update();
					}
					
					for(int i=0; i<numSections; i++){
						trueSolution[i].updatemeasurement();
					}
					
					for(int i=0; i<numSections1; i++){
						trueSolution1[i].update();
					}
					
					for(int i=0; i<numSections1; i++){
						trueSolution1[i].updatemeasurement();
					}

					for (int i=0; i<numSections; i++) {
						estimations[i].nextStepWithoutUpdateTrue();
					}
					
					for (int i=0; i<numSections1; i++) {
						estimations1[i].nextStepWithoutUpdateTrue();
					}
					
					for (int i = 0; i<numSections; i++) {
						S[i]=(estimations[i].filter.measure.transpose()).mmul(InverseMatrix.invPoSym(estimations[i].filter.measureVar)).mmul(estimations[i].filter.measure);
					}
					
					if (estimations[0].consensus){					
					    DoubleMatrix[] G=new DoubleMatrix[numSections];
					    for (int i=0; i<numSections; i++) {
						    G[i]=estimations[i].roadModel.section.ModelA.mmul(estimations[i].filter.priorVar).mmul(estimations[i].roadModel.section.ModelA.transpose()).add(estimations[i].roadModel.modelVar).add(estimations[i].filter.f_var.mmul(S[i]).mmul(estimations[i].filter.f_var));
					    }
					
					    DoubleMatrix[] Lambda=new DoubleMatrix[numSections];
					    for (int i=0; i<numSections; i++) {
						    Lambda[i]=InverseMatrix.invPoSym(estimations[i].filter.priorVar).sub(estimations[i].roadModel.section.ModelA.transpose().mmul(InverseMatrix.invPoSym(G[i])).mmul(estimations[i].roadModel.section.ModelA));
					    }
					
					    DoubleMatrix BigG=G[0];
					    for (int i=1;i<numSections;i++){
						    BigG=Concat.Diagonal(BigG, G[i]);
					    }
					
					    DoubleMatrix BigLambda=Lambda[0];
					    for (int i=1;i<numSections;i++){
						    BigLambda=Concat.Diagonal(BigLambda, Lambda[i]);
					    }
								
					    DoubleMatrix A_hat=estimations[0].roadModel.section.ModelA;
					    for (int i=1;i<numSections;i++){
						    A_hat=Concat.Diagonal(A_hat, estimations[i].roadModel.section.ModelA);
					    }
					
					    DoubleMatrix Lower=A_hat.transpose().mmul(L_tilde.transpose()).mmul(H_hat).mmul(BigG).mmul(H_hat.transpose()).mmul(L_tilde).mmul(A_hat);
					    double[][] _BigLambda=new double[BigLambda.rows][BigLambda.columns];
					    double[][] _Lower=new double[Lower.rows][Lower.columns];
					    for (int i=0;i<BigLambda.rows;i++){
						    for (int j=0;j<BigLambda.columns;j++){
							    _BigLambda[i][j]=BigLambda.get(i, j);
						    }
					    }
					    for (int i=0;i<Lower.rows;i++){
						    for (int j=0;j<Lower.columns;j++){
							    _Lower[i][j]=Lower.get(i, j);
						    }
					    }
					    Matrix bigLambda=new Matrix(_BigLambda);
					    Matrix lower=new Matrix(_Lower);
				    	EigenvalueDecomposition bigLambdaEig=new EigenvalueDecomposition(bigLambda);
					    EigenvalueDecomposition lowerEig=new EigenvalueDecomposition(lower);
					    Matrix _bigLambda=bigLambdaEig.getD();
					    Matrix _lower=lowerEig.getD();
					
					    double temp1=_bigLambda.get(0, 0);
					    for(int i=1;i<_bigLambda.getColumnDimension();i++){
						    if(_bigLambda.get(i, i)<temp1){
							    temp1=_bigLambda.get(i, i);
						    }						
					    }
					
					    double temp2=_lower.get(0, 0);
					    for(int i=1;i<_lower.getColumnDimension();i++){
						    if(_lower.get(i, i)>temp2){
							    temp2=_lower.get(i, i);
						    }
						}
					    double _gamma=Math.sqrt(temp1/temp2);
					    Gamma.write(_gamma+",");
					    Gamma.write(temp1+",");
					    Gamma.write(temp2+",");
					    Gamma.write("\n");	
					    System.out.print(_gamma+"\n");
				
					    for (int i=0; i<numSections; i++) {
				            if(estimations[i].roadModel.section.getClass().getSimpleName().equals("FC")){
//					            estimations[i].filter.gamma=0;
				    	        estimations[i].filter.gamma=0.999999*_gamma;
				            }
				            else{
					            estimations[i].filter.gamma=0.999999*_gamma;
//			                    estimations[i].filter.gamma=0;
				            }						
					    }
					
					DoubleMatrix[] mean1=new DoubleMatrix[numSections];
//					if (estimations[1].roadModel.section.getClass().getSimpleName().equals("FC")){
//						mean1[0]=estimations[0].filter.mean;
//					}
//					else{
						mean1[0]=estimations[0].filter.Consensus2(estimations[1].filter);
//					}
					
					for (int i=1; i<numSections-1; i++) {
//						if (estimations[i-1].roadModel.section.getClass().getSimpleName().equals("FC")){
//							if (estimations[i+1].roadModel.section.getClass().getSimpleName().equals("FC")){
//								mean1[i]=estimations[i].filter.mean;
//							}
//							else{
//								 mean1[i]=estimations[i].filter.Consensus2(estimations[i+1].filter);
//							}
//						}
//						else{
//                          if (estimations[i+1].roadModel.section.getClass().getSimpleName().equals("FC")){
//                           	 mean1[i]=estimations[i].filter.Consensus1(estimations[i-1].filter);
//							}
//                          else{
                            	mean1[i]=estimations[i].filter.Consensus(estimations[i-1].filter,estimations[i+1].filter);
//                          }
//						}
					}
//                  if (estimations[numSections-2].roadModel.section.getClass().getSimpleName().equals("FC")){
//                   	mean1[numSections-1]=estimations[numSections-1].filter.mean;
//					}
//                  else{
                    	mean1[numSections-1]=estimations[numSections-1].filter.Consensus1(estimations[numSections-2].filter);
//                  }									
					for (int i=0; i<numSections; i++) {
						estimations[i].filter.mean=mean1[i];
					}
				}
			}
				
			for (int i = 0; i<numSections; i++) {
				writerSection[i].flush(); writerSection[i].close();
				Mode[i].flush(); Mode[i].close();
			}
			for (int i = 0; i<numSections1; i++) {
				writerSection1[i].flush(); writerSection1[i].close();
				Mode1[i].flush(); Mode1[i].close();
			}
			writerAvr.flush(); writerAvr.close();
			writerAvr1.flush(); writerAvr1.close();
			Error.flush(); Error.close();
			Error1.flush(); Error1.close();
            writerTrue.flush(); writerTrue.close();
			Gamma.flush(); Gamma.close();
		    writerDisagreement.flush(); writerDisagreement.close();
			writerDisagreement1.flush(); writerDisagreement1.close();
			
			System.out.println(" End");		
		}
		catch (Exception e) {e.printStackTrace();}		
	}
	
}
