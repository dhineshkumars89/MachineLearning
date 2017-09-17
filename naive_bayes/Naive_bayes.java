package assignment3;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;
import java.util.concurrent.TimeUnit;

import assignment3.Attribute;

public class Naive_bayes {

	public static void main(String[] args) {
//		long startTime = System.currentTimeMillis();
		
		try {
			BufferedReader trainingfile = new BufferedReader (new FileReader(args[0]));
			BufferedReader testfile = new BufferedReader (new FileReader(args[1]));
			String method = args[2];
			Integer extraoption;
			if(args.length<4)
				extraoption=0;
			else
				extraoption = Integer.parseInt(args[3]);
			if(method.isEmpty())
				System.out.println("Invalid method type");
			if(method.compareToIgnoreCase("histograms") == 0)
			{
				CreateHistogram(trainingfile,testfile,extraoption);
			}
			else if(method.compareToIgnoreCase("gaussians") == 0)
			{
				CreateGaussian(trainingfile,testfile);
			}
			else if(method.compareToIgnoreCase("mixtures") == 0)
			{
				CreateMixtures(trainingfile,testfile,extraoption);
			}
			else
				System.out.println("Invalid method type");
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

//		long endTime   = System.currentTimeMillis();
//		long totalTime = endTime - startTime;
//		 System.out.println(totalTime);
//		 System.out.println(TimeUnit.MILLISECONDS.toSeconds(totalTime));
	}

	private static void CreateMixtures(BufferedReader trainingfile, BufferedReader testfile, Integer extraoption) {
		// TODO Auto-generated method stub
		HashMap<Integer,ArrayList> hm = new HashMap();
		int col = 0; //no of columns
	    double total_records=0; //no of records
	    
	    ///Put all records of particular class type into arraylist  then into HashMap(classtype,arraylist)
	    try {
	    	String s = trainingfile.readLine();
			do{
				  String[] bit = s.trim().split("\\s+");
				  int classtype = Integer.valueOf(bit[bit.length-1]);
				  col = bit.length-1;
				  ArrayList lm;
				  
				  if(hm.containsKey(classtype))
				  {
					  lm = hm.get(classtype);
					  lm.add(Arrays.copyOfRange(bit,0,bit.length-1));
					  hm.put(classtype, lm);
				  }
				  else
				  {
					  lm = new ArrayList();
					  lm.add(Arrays.copyOfRange(bit,0,bit.length-1));
					  hm.put(classtype, lm);
				  }
				  s = trainingfile.readLine();
				  total_records++;
			}while( s != null );
		
			} catch (IOException e) {
			e.printStackTrace();
		}
		
	    // to store key and their respective attribute ALSO in hashmap
	    HashMap<Integer,HashMap> classhm= new HashMap();
	    HashMap<Integer,Double> classprob = new HashMap();
	    
	    Iterator entries = hm.entrySet().iterator();
		while (entries.hasNext()) {
			//classhm = new HashMap();
			Map.Entry entry = (Map.Entry) entries.next();
		    ArrayList arylst = (ArrayList) entry.getValue();
		    Integer key = (Integer) entry.getKey();
			int class_size = arylst.size();
			
			classprob.put(key, ((double)arylst.size()/total_records));
			
			double[][] classary = new double[class_size][col];  //record row and their columns
			double[] max = new double[col]; //max in columns
			double[] min = new double[col]; //min in columns
			
			
			//initialise max and min for class type (key)
			for(int i=0;i<col;i++)
			{
					max[i]=Double.MIN_VALUE;
					min[i]=Double.MAX_VALUE;
			}
			int k=0;
			Iterator its = arylst.iterator();
			while(its.hasNext())
			{
				String[] str = (String[]) its.next();
				for(int i=0;i<col;i++)
				{
					double val = Double.valueOf(str[i]);
					classary[k][i] = val;
					
					//store max n min value for attribute i
					if(val > max[i])
						max[i] = val;
					if((val < min[i]))
						min[i] = val;
				}
				k++;
			}
			
			
			//Number of gaussian
			int N = extraoption;
			double[][] p=new double[N][class_size];
			double sum = 0.0;
			
			//To store atr id and their values in hashmap
			HashMap<Integer,Attribute> attributehm = null; 
			
			for(int b=0;b<col;b++)
			{
				//intialization step
				if(attributehm == null)
					attributehm = new HashMap();
				else
					attributehm = classhm.get(key);
				Attribute atr = new Attribute(b,min[b],max[b],N);
				double[] W = atr.W;
				double[] mean=atr.mean;
				double[] sd=atr.sd;
				
				//EM algo loop
				for(int t=0;t<50;t++)
				{
					//Expectation Step - start
					
					//loop through records
					for(int a=0;a<class_size;a++)
					{
						sum = 0.0;
						for(int i=0;i<N;i++)
						{
							double first,last,result;
							if(sd[i] < Math.sqrt(0.0001))
								sd[i] = Math.sqrt(0.0001);
							first = 1.0/(sd[i] * Math.sqrt(2*Math.PI));
							last = Math.pow((classary[a][b]-mean[i]),2) / (2 * Math.pow(sd[i],2));
							result = first * Math.pow(Math.E, -last);
							//compute_gaussian(atr.mean[i],atr.sd[i],classary[a][b])
							p[i][a] = (W[i] * result);
							if(Double.isInfinite(p[i][a]) || (Double.isNaN(p[i][a])))
								p[i][a]=0.0;
								
							sum = sum + p[i][a];
							//System.out.println("itr:"+t+" classrow:"+a+" gaussian:"+i+" value::"+classary[a][b]+" "+first+"-"+last+"-"+result+"-"+p[i][a]+"-"+sum);
						}
						for(int i=0;i<N;i++)
						{
							p[i][a] = p[i][a]/sum;
							if(Double.isInfinite(p[i][a]) || (Double.isNaN(p[i][a])))
								p[i][a]=0.0;
						}
					}
					//Expectation Step - end
					
					//Maximization step - start
					double[] summ = new double[N];
					double[] msum = new double[N];
					double[] sdsum = new double[N];
					
					//loop through records
					for(int a=0;a<class_size;a++)
					{	
						sum = 0.0;
						for(int i=0;i<N;i++)
						{
							if(Double.isNaN(p[i][a]))
							{ 
								p[i][a]=0.0;
							}
							msum[i] = msum[i] + (p[i][a] * classary[a][b]);
							summ[i] = summ[i] + p[i][a];
						    sum = sum + summ[i];
						    sdsum[i] = sdsum[i] + (p[i][a] * (Math.pow((classary[a][b]-mean[i]),2)));
						    //System.out.println("itr:"+t+" classrow:"+a+" gaussian:"+i+" value::"+classary[a][b]+" "+msum[i]+"-"+summ[i]+"-"+sdsum[i]+"-"+p[i][a]+"-"+sum);
					 	}
					
						for(int i=0;i<N;i++)
						{
							mean[i] = msum[i]/summ[i];
							if(Double.isInfinite(mean[i]) || Double.isNaN(mean[i]))
								mean[i]=0.0;
							sd[i] = Math.sqrt(sdsum[i]/summ[i]);
							if((sd[i] < Math.sqrt(0.0001)) || Double.isInfinite(sd[i]) || Double.isNaN(sd[i]))
								sd[i] = Math.sqrt(0.0001);
							W[i] = summ[i]/sum;
							if(Double.isInfinite(W[i]) || Double.isNaN(W[i]))
								W[i]=0.0;
							//System.out.println(String.format("%.4f		%.4f	%.4f",W[i],sd[i],mean[i]));
						}
					}
					//Maximization step - stop
					
					///End of 50 loop  
				}
				//save the calulated values for mean,sd,weight
				atr.setMean(mean);
				atr.setSd(sd);
				atr.setW(W);
				attributehm.put(b,atr);
				classhm.put(key,attributehm);
			}
		}
		
		
		///Printing the Output
		Iterator entri = classhm.entrySet().iterator();
		while (entri.hasNext()) {
			Map.Entry entry = (Map.Entry) entri.next();
		    HashMap arylst = (HashMap) entry.getValue();
		    Integer key = (Integer) entry.getKey();
		    Iterator enty = arylst.entrySet().iterator();
		    while (enty.hasNext()) {
				Map.Entry entry1 = (Map.Entry) enty.next();
				Attribute arg = (Attribute) entry1.getValue();
			    Integer key1 = (Integer) entry1.getKey();
			    for(int r=0;r<arg.mean.length;r++)
			    {
			    	System.out.printf("Class %d, attribute %d, Gaussian %d, mean = %.2f, std = %.2f",key,key1,r,arg.mean[r],arg.sd[r]);
			    	System.out.println();
			    }
			}
		}
		
		
		
		//Classification part
		HashMap<Integer,Double> classconprob;// = new HashMap();
		int row=0;
		double conprobability;
		double numcorrectclassification=0.0;
	    try {
	    	String s = testfile.readLine();
			do{
				  String[] bit = s.trim().split("\\s+");
				  classconprob = new HashMap();
				  for(int j=0;j<bit.length-1;j++)
				  {
					double attribute = Double.valueOf(bit[j]);
					Iterator clssloop = classhm.entrySet().iterator();
					while (clssloop.hasNext()) {
						Map.Entry entry = (Map.Entry) clssloop.next();
					    HashMap arylst = (HashMap) entry.getValue();
					    Integer key = (Integer) entry.getKey();
					    Attribute val = (Attribute)arylst.get(j);
					    conprobability = 0.0;
					    for(int r=0;r<val.mean.length;r++)
					    {
					    	double first,last,result;
							first = 1.0/(val.sd[r] * Math.sqrt(2.0*Math.PI));
							last = Math.pow((attribute-val.mean[r]),2) / (2.0 * Math.pow(val.sd[r],2));
							result = first * Math.pow(Math.E, -last);
					    	conprobability = conprobability + (val.W[r] * result);
					    }
					    if(classconprob.containsKey(key))
					    {	
					    	Double intermin = (double)classconprob.get(key);
					    	intermin = conprobability * intermin.doubleValue();
					    	classconprob.put(key, intermin);
					    }
					    else
					    {
					    	classconprob.put(key, conprobability);
					    }
					    
					 }
				  }
		  
		  
				  //prediction for single record
				  double totalprob=0.0,max=0.0;
				  int maxc=0;
				  HashMap<Integer,Double> classval = new HashMap();
				  HashMap<Integer,Double> tieclass = new HashMap();
				  Iterator en = classprob.entrySet().iterator();
				  while (en.hasNext()) {
					    Map.Entry entry = (Map.Entry) en.next();
					    Double clsprob = (Double) entry.getValue();
					    Integer key = (Integer)entry.getKey();
					    
					    Double val = classconprob.get(key).doubleValue() * clsprob.doubleValue();
					    classval.put(key,val);
					    totalprob = totalprob + val.doubleValue();
				  }
			
				  Iterator it = classval.entrySet().iterator();
				  while (it.hasNext()) {
						Map.Entry entry = (Map.Entry) it.next();
						Integer key = (Integer)entry.getKey();
					    Double probb = (Double) entry.getValue();
					    double probabil = (probb.doubleValue()/totalprob);
					    if(max < probabil)
					    {
					    	max = probabil;
					    	maxc = key;
					    	tieclass.clear();
					    	tieclass.put(key, probabil);
					    }
					    else if(max == probabil)
					    {
					    	tieclass.put(key, probabil);
					    }
				  }
					//decision of max and maxc tie breaking for classes
					double accuracy = 0;
					if(tieclass.size()==1)
					{ 
						//do nothing 
						if(maxc == Integer.parseInt(bit[bit.length-1]))
						{
							accuracy = 1;
						}
					}
					else
					{
						int random = (int )(Math.random() * tieclass.size());
						Iterator iit = tieclass.entrySet().iterator();
						int ith = 0;
						while (it.hasNext()) {
							Map.Entry entry = (Map.Entry) iit.next();
							Integer key = (Integer)entry.getKey();
						    Double prob = (Double) entry.getValue();
						    if(ith == random)
						    {
						    	maxc=key;
						    	max=prob;
						    }
						    if(key == Integer.parseInt(bit[bit.length-1]))
							{
								accuracy = 1/tieclass.size();
							}
						}
						
					}
					System.out.printf("ID=%5d, predicted=%3d, probability = %.4f, true=%3d, accuracy=%4.2f\n", 
	                  row, maxc, max, Integer.parseInt(bit[bit.length-1]), accuracy);
					numcorrectclassification = numcorrectclassification + accuracy;
		  
		  
					  s = testfile.readLine();
					  row++;
			}while( s != null );
		
		} catch (IOException e) {
			e.printStackTrace();
		}
	    double classification_accuracy = (numcorrectclassification/(double)row);
	    System.out.printf("classification accuracy=%6.4f\n", classification_accuracy);
	}
	//mean,sd,weight
	private static double compute_gaussian(double m, double sd, double x) {
		// TODO Auto-generated method stub
		
		double first,last,result;
		first = 1.0/(sd * Math.sqrt(2*Math.PI));
		last = Math.pow((x-m),2) / (2 * Math.pow(sd,2));
		result = first * Math.pow(Math.E, -last);
		if(Double.isNaN(first) || Double.isNaN(last))
		{
//			System.out.println((Math.pow((x-m),2)+"-----"+Math.pow(sd,2))+"----"+x+"--"+m+"sd:"+sd);
		}
		return result;
	}

	private static void CreateGaussian(BufferedReader trainingfile, BufferedReader testfile) {

		//hashmap to store records as per class
		HashMap<Integer,ArrayList> hm = new HashMap();

		int col = 0; //no of columns
	    double total_records=0; //no of records
	    
	    //Store record of particular classtype in hashmap
	    try {
	    	String s = trainingfile.readLine();
			do{
			  String[] bit = s.trim().split("\\s+");
			  int classtype = Integer.valueOf(bit[bit.length-1]);
			  col = bit.length-1;
			  ArrayList lm;
			  
			  if(hm.containsKey(classtype))
			  {
				  lm = hm.get(classtype);
				  lm.add(Arrays.copyOfRange(bit,0,bit.length-1));
				  hm.put(classtype, lm);
			  }
			  else
			  {
				  lm = new ArrayList();
				  lm.add(Arrays.copyOfRange(bit,0,bit.length-1));
				  hm.put(classtype, lm);
			  }
			  s = trainingfile.readLine();
			  total_records++;
			}while( s != null );
		
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	    
	    
	    double[][] total = new double[hm.size()][col]; // matrix of total [no of class type][attributes]
	    HashMap<Integer,Integer> clas2indx = new HashMap(); //class to index mapping

		int k=0;
		Iterator entries = hm.entrySet().iterator();
		while (entries.hasNext()) {
		    Map.Entry entry = (Map.Entry) entries.next();
		    ArrayList arylst = (ArrayList) entry.getValue();  //arraylist of records of particular type
		    
		    clas2indx.put((Integer) entry.getKey(), k);
			Iterator its = arylst.iterator();
			while(its.hasNext())
			{
				String[] str = (String[]) its.next();
				for(int i=0;i<col;i++)
				{
					total[k][i] = total[k][i] + Double.valueOf(str[i]);
				}
			}
			if(!its.hasNext())
			{
				for(int i=0;i<col;i++)
				{
					total[k][i] = ((total[k][i])/((double)arylst.size()));
					//mean of each attribute for each class
				}
			}
			k++;
		}
		
		
		//Calculate variance
		k=0;
		
		HashMap<Integer,HashMap> classmap = new HashMap();
		HashMap<Integer,Double[]> attrmap; //stores attribute probabilities
		
		HashMap<Integer,Double> classprob = new HashMap();  ///stores class probability
		double[][] variance = new double[hm.size()][col]; // matrix of variance {no of class type][attributes]
		
		Iterator entry1 = hm.entrySet().iterator();
		while (entry1.hasNext()) {
		    Map.Entry entry = (Map.Entry) entry1.next();
		    ArrayList arylst = (ArrayList) entry.getValue();
		    Integer key = (Integer)entry.getKey();
		    
		    attrmap = new HashMap();
		    classprob.put(key, (double) ((double)arylst.size()/total_records));
		    
			Iterator its = arylst.iterator();
			while(its.hasNext())
			{
				String[] str = (String[]) its.next();
				for(int i=0;i<col;i++)
				{
					variance[k][i] = variance[k][i] + Math.pow((Double.valueOf(str[i])-total[clas2indx.get(key)][i]),2);
				}
			}
			if(!its.hasNext())
			{
				for(int i=0;i<col;i++)
				{
					Double[] values = new Double[2];
					variance[k][i] = ((variance[k][i])/((double)arylst.size()-1.0));
					if(variance[k][i] < 0.0001) { variance[k][i] = 0.0001;}
//					System.out.println(String.format("Class %d, attribute %d, mean = %.2f, std = %.2f",key,i+1,total[clas2indx.get(key)][i],Math.sqrt(variance[k][i])));
					System.out.printf("Class %d, attribute %d, mean = %.2f, std = %.2f",key,i,total[clas2indx.get(key)][i],Math.sqrt(variance[k][i]));
					System.out.println();
					values[0]=total[clas2indx.get(key)][i];
					values[1]=Math.sqrt(variance[k][i]);
					attrmap.put(i, values); //stores values in atr hashmap
					
				}
			}
			k++;
			
			classmap.put(key, attrmap); // stores atr hash map into their class map
		}
		
		//Class probability p(c)
		//computed while reading data
		
		//TRAINING 
		//Classification
		HashMap<Integer,ArrayList> trainhm = new HashMap();
		HashMap<Integer,Double> classconprob; // = new HashMap();

		int row=0;
		double numcorrectclassification=0.0;
	    try {
	    	String s = testfile.readLine();
			do{
			  String[] bit = s.trim().split("\\s+");
			  classconprob = new HashMap();  //stores class conditionally probability
		  
			  //loop thro attributes
			  for(int j=0;j<bit.length-1;j++)
			  {
					double attribute = Double.valueOf(bit[j]);
					Iterator en = classmap.entrySet().iterator();
					while (en.hasNext()) {
					    Map.Entry entry = (Map.Entry) en.next();
					    Integer key = (Integer)entry.getKey();
					    HashMap atr = (HashMap) entry.getValue();
					    Double[] values =  (Double[]) atr.get(j);
					    double mean=values[0],sd=values[1];
					    
					    //formula
					  	double first,last,result;
						first = 1.0/(sd * Math.sqrt(2.0*Math.PI));
						last = Math.pow((attribute-mean),2) / (2.0 * Math.pow(sd,2));
						result = first * Math.pow(Math.E, -last);
		
						if(classconprob.containsKey(key))
						{
							Double val = (Double) classconprob.get(key);
							val = val * result;
							classconprob.put(key, val);
						}
						else
						{
							classconprob.put(key, result);
						}
					    
					}
			  }
		  
			  ///process for each row
			  double totalprob=0.0,max=0.0;
			  int maxc=0;
			  HashMap<Integer,Double> classval = new HashMap();
			  HashMap<Integer,Double> tieclass = new HashMap();
			  Iterator en = classprob.entrySet().iterator();
			  while (en.hasNext()) {
				    Map.Entry entry = (Map.Entry) en.next();
				    Double clsprob = (Double) entry.getValue();
				    Integer key = (Integer)entry.getKey();
				    
				    Double val = classconprob.get(key).doubleValue() * clsprob.doubleValue();
				    classval.put(key,val);
				    totalprob = totalprob + val.doubleValue();
			  }
			
			  Iterator it = classval.entrySet().iterator();
			  while (it.hasNext()) {
					Map.Entry entry = (Map.Entry) it.next();
					Integer key = (Integer)entry.getKey();
				    Double probb = (Double) entry.getValue();
				    double probabil = (probb.doubleValue()/totalprob);
				    if(max < probabil)
				    {
				    	max = probabil;
				    	maxc = key;
				    	tieclass.clear();
				    	tieclass.put(key, probabil);
				    }
				    else if(max == probabil)
				    {
				    	tieclass.put(key, probabil);
				    }
			  }
				//decision of max and maxc tie breaking for classes
				double accuracy = 0;
				if(tieclass.size()==1)
				{ 
					//do nothing 
					if(maxc == Integer.parseInt(bit[bit.length-1]))
					{
						accuracy = 1;
					}
				}
				else
				{
					int random = (int )(Math.random() * tieclass.size());
					Iterator iit = tieclass.entrySet().iterator();
					int ith = 0;
					while (it.hasNext()) {
						Map.Entry entry = (Map.Entry) iit.next();
						Integer key = (Integer)entry.getKey();
					    Double prob = (Double) entry.getValue();
					    if(ith == random)
					    {
					    	maxc=key;
					    	max=prob;
					    }
					    if(key == Integer.parseInt(bit[bit.length-1]))
						{
							accuracy = 1/tieclass.size();
						}
					}
					
				}
				System.out.printf("ID=%5d, predicted=%3d, probability = %.4f, true=%3d, accuracy=%4.2f\n", 
		                  row, maxc, max, Integer.parseInt(bit[bit.length-1]), accuracy);
				numcorrectclassification = numcorrectclassification + accuracy;

				s = testfile.readLine();
				row++;
			}while( s != null );
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	    double classification_accuracy = (numcorrectclassification/(double)row);
	    System.out.printf("classification accuracy=%6.4f\n", classification_accuracy);
	}

	private static void CreateHistogram(BufferedReader trainingfile, BufferedReader testfile, Integer extraoption) {
		// TODO Auto-generated method stub
				
	    int col = 0;//number of columns
	    double total_records=0; //no of records
	    
	    //Read file lines and insert records of class type into their respective hashmap key (class type)
	    HashMap<Integer,ArrayList> hm = new HashMap();
	    try {
	    	String s = trainingfile.readLine();
			do{
			  String[] bit = s.trim().split("\\s+");
			  int classtype = Integer.valueOf(bit[bit.length-1]);
			  col = bit.length-1;
			  ArrayList lm;
			  
			  if(hm.containsKey(classtype))
			  {
				  lm = hm.get(classtype);
				  lm.add(Arrays.copyOfRange(bit,0,bit.length-1));
				  hm.put(classtype, lm);
			  }
			  else
			  {
				  lm = new ArrayList();
				  lm.add(Arrays.copyOfRange(bit,0,bit.length-1));
				  hm.put(classtype, lm);
			  }
			  
				s = trainingfile.readLine();
				total_records++;
			}while( s != null );
		
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	    HashMap<Integer,Double> classprob = new HashMap();
	    HashMap<Integer,HashMap> classhash = new HashMap();
	    //Iterate through class types and their records
	    Iterator entries = hm.entrySet().iterator();
		while (entries.hasNext()) {
			Map.Entry entry = (Map.Entry) entries.next();
			
			//key - class type, arylst - list of records of that type
		    ArrayList arylst = (ArrayList) entry.getValue();
		    Integer key = (Integer) entry.getKey();
		    //size of class records
			int class_size = arylst.size();
			
			classprob.put(key, ((double)arylst.size()/(total_records)));
			
			double[][] classary = new double[class_size][col];
			double[] max = new double[col];
			double[] min = new double[col];
			int k=0;
			
			for(int i=0;i<col;i++)
			{
					max[i]=Double.MIN_VALUE;
					min[i]=Double.MAX_VALUE;
			}
			Iterator its = arylst.iterator();
			while(its.hasNext())
			{
				String[] str = (String[]) its.next();
				for(int i=0;i<col;i++)
				{
					double val = Double.valueOf(str[i]);
					classary[k][i] = val;
					
					
					if(val > max[i])
						max[i] = val;
					if((val < min[i]))
						min[i] = val;
				}
				k++;
			}
			HashMap<Integer,HistoAttr> attrhash = new HashMap();
			HistoAttr hstatr;
			for(int j=0;j<col;j++)
			{
				double L = max[j];
				double S = min[j];
				int N = extraoption;
				hstatr = new HistoAttr(j,S,L,N);
				
				double G = (L-S)/(N-3);
				if((G<0.0001) || Double.isInfinite(G) || Double.isNaN(G))
					G=0.0001;
				int[] bin = new int[extraoption];
				
				for(int i=0;i<class_size;i++)
				{
					if((classary[i][j] < (S-(G/2.0)) )) // || (G==0.0)
					{
						bin[0] = bin[0] + 1;
					}
					else if((classary[i][j] >= (L+(G/2.0))))// && (N>3))
					{
						bin[N-1] = bin[N-1] + 1;
					}
					else if((classary[i][j] >= (S-(G/2.0))) && (classary[i][j] < (S+(G/2.0))))// && (N>=2))
					{
						bin[1] = bin[1] + 1;
					}
					else
					{
						for(int t=2;t<(N-1);t++)
						{
							if(  ( classary[i][j] >= (S+((t-2)*G)+(G/2.0)) ) && ( classary[i][j] < (S+((t-1)*G)+(G/2.0)) ) )
							{
								bin[t] = bin[t]+1;
							}
						}
					}
					
				}
				for(int t=0;t<N;t++)
				{
					//((double)(bin[t]/class_size))
 					double prob = ((bin[t])/(double)(class_size * G)); // 
					System.out.printf("Class %d, attribute %d, bin %d, P(bin | class) = %.2f",key,j,t,prob);
					System.out.println();
					hstatr.setidvBin(t, prob);
				}
				attrhash.put(j, hstatr);
			}
			classhash.put(key, attrhash);
		}
		
		
		//Classification part
		HashMap<Integer,Double> classconprob; // = new HashMap();
//	    double classprob[] = new double[hm.size()];
		int row=0;
		double numcorrectclassification=0.0;
	    try {
	    	String s = testfile.readLine();
			do{
			  String[] bit = s.trim().split("\\s+");
			  classconprob = new HashMap();
			  
			  for(int j=0;j<bit.length-1;j++)
			  {
				  double attribute = Double.valueOf(bit[j]);
				  
					Iterator en = classhash.entrySet().iterator();
					while (en.hasNext()) {
					    Map.Entry entry = (Map.Entry) en.next();
					    Integer key = (Integer)entry.getKey();
					    HashMap atr = (HashMap) entry.getValue();
					    HistoAttr histoatr =  (HistoAttr) atr.get(j);
					    double result = 0;
					    double S =histoatr.min;
					    double L =histoatr.max;
					    int N = histoatr.N;
					    double G = (L-S)/(N-3);
					    if((G<0.0001) || Double.isInfinite(G) || Double.isNaN(G))
							G=0.0001;
												
						for(int i=0;i<histoatr.N;i++)
						{
							if((attribute < (S-(G/2.0)))) // || (G==0.0)
							{
								result = histoatr.bin[0];
								break;
							}
							else if((attribute >= (S-(G/2.0))) && (attribute < (S+(G/2.0))))// && (N >= 2))
							{
								result = histoatr.bin[1];
								break;
							}
							else if((attribute >= (L+(G/2.0))))// && (N > 3))
							{
								result = histoatr.bin[N-1];
								break;
							}
							else
							{
								for(int t=2;t<(N-1);t++)
								{
									if(  ( attribute >= (S+((t-2)*G)+(G/2.0)))  && ( attribute < (S+((t-1)*G)+(G/2.0)) ) )
									{
										result = histoatr.bin[t];
										break;
									}
								}
							}
							
						} 
					   

						if(classconprob.containsKey(key))
						{
							Double val = (Double) classconprob.get(key);
							val = val * result;
							classconprob.put(key, val);
						}
						else
						{
							classconprob.put(key, result);
						}
					    
					}
				  
			  }
			  
			///process for each row
			  double totalprob=0.0,max=0.0;
			  int maxc=0;
			    HashMap<Integer,Double> classval = new HashMap();
			    HashMap<Integer,Double> tieclass = new HashMap();
			    Iterator en = classprob.entrySet().iterator();
				while (en.hasNext()) {
				    Map.Entry entry = (Map.Entry) en.next();
				    Double clsprob = (Double) entry.getValue();
				    Integer key = (Integer)entry.getKey();
				    
				    Double val = classconprob.get(key).doubleValue() * clsprob.doubleValue();
				    classval.put(key,val);
				    totalprob = totalprob + val.doubleValue();
				}
			  	Iterator it = classval.entrySet().iterator();
				while (it.hasNext()) {
					Map.Entry entry = (Map.Entry) it.next();
					Integer key = (Integer)entry.getKey();
				    Double probb = (Double) entry.getValue();
				    double probabil = (probb.doubleValue()/totalprob);
				    if(max < probabil)
				    {
				    	max = probabil;
				    	maxc = key;
				    	tieclass.clear();
				    	tieclass.put(key, probabil);
				    }
				    else if(max == probabil)
				    {
				    	tieclass.put(key, probabil);
				    }
				}
			  
			  	//decision of max and maxc tie breaking for classes
				double accuracy = 0;
				if(tieclass.size()==1)
				{ 
					//do nothing 
					if(maxc == Integer.parseInt(bit[bit.length-1]))
					{
						accuracy = 1;
					}
				}
				else
				{
					int random = (int )(Math.random() * tieclass.size());
					Iterator iit = tieclass.entrySet().iterator();
					int ith = 0;
					while (it.hasNext()) {
						Map.Entry entry = (Map.Entry) iit.next();
						Integer key = (Integer)entry.getKey();
					    Double prob = (Double) entry.getValue();
					    if(ith == random)
					    {
					    	maxc=key;
					    	max=prob;
					    }
					    if(key == Integer.parseInt(bit[bit.length-1]))
						{
							accuracy = 1/tieclass.size();
						}
					}
					
				}
			  System.out.printf("ID=%5d, predicted=%3d, probability = %.4f, true=%3d, accuracy=%4.2f\n", 
	                  row, maxc, max, Integer.parseInt(bit[bit.length-1]), accuracy);
			  numcorrectclassification = numcorrectclassification + accuracy;
			  s = testfile.readLine();
			  row++;
			}while( s != null );
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	    double classification_accuracy = (numcorrectclassification/(double)row);
	    System.out.printf("classification accuracy=%6.4f\n", classification_accuracy);
	}

}

