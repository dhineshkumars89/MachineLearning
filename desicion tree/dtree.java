
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Random;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

public class dtree {

	static double[] classes = null;
	static Integer[] attributes = null;
	TreeMap<Integer,node> hashnode = new TreeMap<Integer,node>();
	static Random rand = new Random();
	static double classification_acc = 0;
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String training_file = args[0];
		String test_file = args[1];
		String mode = args[2];
		int pruning_thr = Integer.parseInt(args[3]);
		
		
//		String training_file = "D:\\CSE 6363 Machine Learning\\Assignment\\Assignment 2\\pendigits_training.txt";
//		String test_file = "D:\\CSE 6363 Machine Learning\\Assignment\\Assignment 2\\pendigits_test.txt";
		
		dtree D = new dtree();
		double[][] train = D.getdata(training_file);
		double[][] test = D.getdata(test_file);
		
		D.classes = dtree.getClasses(train);
		D.attributes = dtree.getAttributes(train);
		
		if(mode.equalsIgnoreCase("optimized"))
		{
			node N = null;
			N = D.DTL(0,mode,train,D.attributes,distribution(train),pruning_thr,0);
			printDTL(D.hashnode);
			D.classification(test,N);
		}
		else if(mode.equalsIgnoreCase("randomized"))
		{
			node N = null;
			N = D.DTL(0,mode,train,D.attributes,distribution(train),pruning_thr,0);
			printDTL(D.hashnode);
			D.classification(test,N);
		}
		else if(mode.equalsIgnoreCase("forest3"))
		{
			node N[] = new node[3];
			HashMap hm = new HashMap();
			for(int i=0;i<3;i++)
			{
				dtree G = new dtree();
				N[i] = G.DTL(i,mode,train,D.attributes,distribution(train),pruning_thr,0);
				hm.put(i, G.hashnode);
				printDTL(G.hashnode);
			}
			D.classificationForest(test,N);
			
		}
		else if(mode.equalsIgnoreCase("forest15"))
		{
			node N[] = new node[15];
			HashMap hm = new HashMap();
			for(int i=0;i<15;i++)
			{
				dtree G = new dtree();
				N[i] = G.DTL(i,mode,train,D.attributes,distribution(train),pruning_thr,0);
				hm.put(i, G.hashnode);
				printDTL(G.hashnode);
			}
			D.classificationForest(test,N);
		}
		else
		{
			System.out.println("Invalid mode");
		}
		
		
		
		
	}

	public void classification(double[][] test, node n) {
		// TODO Auto-generated method stub
		double nrows = test.length;
		int cls_ind = test[0].length-1;
		classification_acc = 0;
		int cls_val=-1;
		for(int i=0;i<test.length;i++)
		{
			cls_val = (int)test[i][cls_ind];
			//System.out.println("checking at node "+n.id+" for atr "+n.feature_id+" threshold="+n.threshold);
			double atrval = test[i][(int)n.feature_id];
			if(atrval < n.threshold)
				this.classify(i,test[i],n.left,cls_val);
			else
				this.classify(i,test[i],n.right,cls_val);
		}
		System.out.printf("classification accuracy=%6.4f\n", (classification_acc/nrows));
	}
	
	public void classify(int row_id,double[] testrow,node n,int trueclass)
	{
		double atrval = testrow[(int)n.feature_id];
		while(n.threshold != -1)
		{
			//System.out.println("checking at node "+n.id+" for atr "+n.feature_id+" threshold="+n.threshold+" Value="+atrval);
			if(atrval < n.threshold)
			{
				n = n.left;
			}
			else
			{
				n = n.right;
			}
			if(n.feature_id>=0)
				atrval = testrow[(int)n.feature_id];
		}
		{
			double[] arry = n.prob;
			//System.out.println(arry);
			double max = -1;
			int class_ind = -1;
			for(int r=0;r<arry.length;r++)
			{
				if(max < arry[r])
				{
					max = arry[r];
					class_ind = r;
				}
			}
			if(class_ind>=0)
			{
				System.out.printf("ID=%5d, predicted=%3d, true=%3d, accuracy=%4.2f\n", row_id, (int)this.classes[class_ind], trueclass, max);
				if((int)this.classes[class_ind]==trueclass)
				{
					classification_acc= classification_acc + 1;
				}
			}
		}
			
	}

	public void classificationForest(double[][] test, node[] N) {
		// TODO Auto-generated method stub
		double nrows = test.length;
		int cls_ind = test[0].length-1;
		classification_acc = 0;
		int cls_val=-1;
		for(int i=0;i<test.length;i++)
		{
			cls_val = (int)test[i][cls_ind];
			double[][] probil = new double[N.length][classes.length];
			for(int r=0;r<N.length;r++)
			{
				probil[r] = this.classifyForest(i,test[i],N[r],cls_val);
			}
			
			int cls_len = probil[0].length;
			double[] prob = new double[cls_len];
			for(int w=0;w<probil[0].length;w++)
			{	
				for(int r=0;r<N.length;r++)
				{
					prob[w] += probil[r][w] ;
				}
				
			}
			
			
			double max = -1;
			int class_ind = -1;
			for(int w=0;w<prob.length;w++)
			{
				prob[w] = prob[w]/(double)3;
				if(max < prob[w])
				{
					max = prob[w];
					class_ind = w;
				}
			}
			if(class_ind>=0)
			{
				System.out.printf("ID=%5d, predicted=%3d, true=%3d, accuracy=%4.2f\n", i, (int)this.classes[class_ind], cls_val, max);
				if((int)this.classes[class_ind]==cls_val)
				{
					classification_acc= classification_acc + 1;
				}
			}
			

		}
		System.out.printf("classification accuracy=%6.4f\n", (classification_acc/nrows));
	}
	
	
	public double[] classifyForest(int row_id,double[] testrow,node n,int trueclass)
	{
		double atrval = testrow[(int)n.feature_id];
		while(n.threshold != -1)
		{
			if(atrval < n.threshold)
			{
				n = n.left;
			}
			else
			{
				n = n.right;
			}
			if(n.feature_id>=0)
				atrval = testrow[(int)n.feature_id];
		}
		
		return n.prob; 
		
	}
	public static void printDTL(TreeMap hashn) {
		Iterator iterator = hashn.keySet().iterator();

		while (iterator.hasNext()) {
		   Integer key = (Integer) iterator.next();
		   node nd = (node) hashn.get(key);

		   System.out.printf("tree=%2d, node=%3d, feature=%2d, thr=%6.2f, gain=%f \n", nd.tree, nd.id+1, (int)nd.feature_id, nd.threshold, nd.info_gain);
		   //System.out.println(Arrays.toString(nd.prob));
		}
		
			
		
	}

	public node DTL(int tree_num,String mode,double[][] train,Integer[] attributes,double[] deffault,int pruning_thr,int idx)
	{
		node N = new node();
		if(train.length < 50)
		{
	        N.id = idx;
	        N.tree = tree_num;
	        N.feature_id = -1;
	        N.threshold = -1;
	        N.info_gain = -1;
	        N.prob = deffault;
	        this.hashnode.put(idx, N);
		}
	    else if(check_single_class(train))
	    {
	        N.id = idx;
	        N.tree = tree_num;
	        int cls_ind = train[0].length-1;
	        int cls = (int) train[0][cls_ind];
	        double[] defi = new double[classes.length];
	        for(int h=0;h<classes.length;h++)
	        {
	        	if(classes[h] == (double)cls)
	        	{
	        		defi[h] = 1;
	        	}
	        	else
	        	{
	        		defi[h] = 0;
	        	}
	        }
	        N.feature_id = -1;
	        N.threshold = -1;
	        N.info_gain = -1;
	        N.prob = defi;
	        this.hashnode.put(idx, N);
		}
	    else
	    {
	        double[] ans = null;
	        if(mode.equalsIgnoreCase("optimized"))
	        {
	        	ans = chooseAttribute(train, attributes);
	        }
	        else
	        {	
	        	ans = chooseAttributeRandom(train, attributes);
	    	}
	        double best_attribute=ans[0];
	        double best_threshold=ans[1];
	        double max_gain=ans[2];
	        N.id = idx;
	        N.tree = tree_num;
	        N.feature_id = best_attribute;
	        N.threshold = best_threshold;
	        N.info_gain = max_gain;
	        
	        double[][] left_examples = getLessThan(train,(int)best_attribute,best_threshold);
		    double[][] right_examples = getGreaterThanEqualTo(train,(int)best_attribute,best_threshold); 
	        int left_ind = (2 * idx) + 1;
	        int right_ind = (2 * idx) + 2;
	        
	        //fprintf('tree=%2d, node=%3d, feature=%2d, thr=%6.2f, gain=%f \n',tree_num-1, idx, node.feature_id-1, node.threshold, node.info_gain);
	        N.left = this.DTL(tree_num,mode,left_examples, attributes, distribution(left_examples),pruning_thr,left_ind);
	        N.right = this.DTL(tree_num,mode,right_examples, attributes, distribution(right_examples),pruning_thr,right_ind);
	        this.hashnode.put(idx, N);
//	        al.add(left_ind, element);
//	        al.add(right_ind, element);
		}
		return N;
	}
	
	public double[][] getdata(String path)
	{
		double[][] data = null;
		try {
			File file = new File(path);
			FileReader fileReader = new FileReader(file);
			BufferedReader bufferedReader = new BufferedReader(fileReader);
			StringBuffer stringBuffer = new StringBuffer();
			String line;
			ArrayList<String> al = new ArrayList<String>();
			int row=0,column=0;
			while ((line = bufferedReader.readLine()) != null) {
				if(column == 0)
				{
					String ary[] = line.trim().split("\\s+");
					column = ary.length;
				}
				al.add(line);
				row++;
			}
			//fileReader.close();
			data = new double[row][column];
			int i = 0;
			//BufferedReader bufferedReaderr = new BufferedReader(fileReader);
			Iterator iterator = al.iterator();
			while (iterator.hasNext()) {
				line = (String) iterator.next();
				String[] ary = line.trim().split("\\s+");
				for(int j=0;j<column;j++)
				{
					data[i][j] = Double.parseDouble(ary[j]);
				}
				i++;
			}
			fileReader.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		return data;
	}
	
	public static double[] chooseAttribute(double[][] examples,Integer[] attributes)
	{
		double max_gain = -1;
		double best_attribute = -1;
		double best_threshold = -1;
		double threshold = -1,gain = -1;
		    for(int i=0;i<attributes.length;i++)
	    	{
		    	//A = attributes[i];
	    	
		        double[] attribute_values = selectColumn(examples, i);
		        Arrays.sort(attribute_values);
		        double L = attribute_values[0];
        		//min(attribute_values);
		        double M = attribute_values[attribute_values.length-1];
		        //max(attribute_values);
		        threshold = 0;
		        gain = 0;
		        for(int K = 1; K<=50; K++)
        		{
		             threshold = L + ((K*(M-L))/51);
		             gain = Information_Gain(examples, i, threshold);
		             if (gain > max_gain)
		             {
		                 max_gain = gain;
		                 best_attribute = i;
		                 best_threshold = threshold;
		             }
	    		}
		    }
		    double[] atrthrgain = new double[3];
		    atrthrgain[0]=best_attribute;
		    atrthrgain[1]=best_threshold;
		    atrthrgain[2]=max_gain;
		    
		    return atrthrgain;
	}
	
	public static double[] chooseAttributeRandom(double[][] examples,Integer[] attributes)
	{
		double max_gain = -1;
		double best_attribute = -1;
		double best_threshold = -1;
		double threshold = -1,gain = -1;
		int min = 0;
		int max = attributes.length - 1;
		int i = rand.nextInt(max - min + 1) + min;
        double[] attribute_values = selectColumn(examples, i);
        Arrays.sort(attribute_values);
        double L = attribute_values[0];
		//min(attribute_values);
        double M = attribute_values[attribute_values.length-1];
        //max(attribute_values);
        threshold = 0;
        gain = 0;
        for(int K = 1; K<=50; K++)
		{
             threshold = L + ((K*(M-L))/51);
             gain = Information_Gain(examples, i, threshold);
             if (gain > max_gain)
             {
                 max_gain = gain;
                 best_attribute = i;
                 best_threshold = threshold;
             }
		}
		   
		    double[] atrthrgain = new double[3];
		    atrthrgain[0]=best_attribute;
		    atrthrgain[1]=best_threshold;
		    atrthrgain[2]=max_gain;
		    
		    return atrthrgain;
	}
	
	//function distribution = DISTRIBUTION(examples,classes)
	public static double[] distribution(double[][] examples)
	{
		double rowss = examples.length;
		double[] distribution = new double[classes.length];
		double[][] record = null;
		for(int i=0;i<classes.length;i++)
		{
			record = getRecord(examples,classes[i]);
		    double k = record.length;
		    double val = k/rowss;
		    distribution[i] = val;
		}
		return distribution;
		
	}
	public static double[] selectColumn(double[][] examples, int i) {
		// TODO Auto-generated method stub
		ArrayList al = new ArrayList();
		for(int j=0;j<examples.length;j++)
		{
			al.add(examples[j][i]);
		}
		double[] record = new double[al.size()];
		//record = al.toArray();
		Iterator iterator = al.iterator();
		int j = 0;
		while (iterator.hasNext()) {
			record[j] = (double) iterator.next();
			j++;
		}
		return record;
	}

	public static double Information_Gain(double[][] examples, int A, double threshold)
	{
		double gain = -1;
		double intial_entropy = -1;
		intial_entropy = dtree.entropy(examples);
		double[][] left_examples = getLessThan(examples,A,threshold);

	    double[][] right_examples = getGreaterThanEqualTo(examples,A,threshold); 

	    double W1 = (double)left_examples.length/(double)examples.length;
	    double entrophy_left = dtree.entropy(left_examples);
	    double W2 = (double)right_examples.length/(double)examples.length;    
	    double entrophy_right = dtree.entropy(right_examples);
		gain = (intial_entropy - (W1*entrophy_left) - (W2* entrophy_right));
		return gain;
	}
	
		
	public static double[][] getLessThan(double[][] examples, int a, double threshold) {
		// TODO Auto-generated method stub
		int len = examples.length;
		ArrayList<double[]> al = new ArrayList();
		for(int i=0;i<len;i++)
		{
			if(examples[i][a] < threshold)
			{
				al.add(examples[i]);
			}
		}
		double[][] record = new double[al.size()][examples[0].length];
		//record = al.toArray();
		Iterator<double[]> iterator = al.iterator();
		int j = 0;
		while (iterator.hasNext()) {
			record[j] = iterator.next();
			j++;
		}
		return record;
	}
	
	public static double[][] getGreaterThanEqualTo(double[][] examples, int a, double threshold) {
		// TODO Auto-generated method stub
		int len = examples.length;
		ArrayList<double[]> al = new ArrayList();
		for(int i=0;i<len;i++)
		{
			if(examples[i][a] >= threshold)
			{
				al.add(examples[i]);
			}
		}
		double[][] record = new double[al.size()][examples[0].length];
		//record = al.toArray();
		Iterator<double[]> iterator = al.iterator();
		int j = 0;
		while (iterator.hasNext()) {
			record[j] = iterator.next();
			j++;
		}
		return record;
	}

	

	public static double entropy(double[][] examples)
	{
		    double tr = examples.length;
		    double entropy = 0;
		    double sum = 0;
		    double val = 0;
		    for(int i=0;i<classes.length;i++)
		    {
		    	val = 0;
		        double cls = classes[i];
		        double[][] record = getRecord(examples,cls);
		        double k = record.length;
		        if(k==0.0)
		        	val = 0;
		        else
		        	val = ((-1 * (k/tr)) * dtree.log2(k/tr));
		        
		        sum = sum + val;
		    }
		    entropy = sum;
		    return entropy;
	}
	
	public static double[][] getRecord(double[][] examples, double cls) {
		// TODO Auto-generated method stub
		ArrayList<double[]> al = new ArrayList();
		double[][] record = null;
		//System.out.println(examples.length);
		if(examples.length==0)
		{
			record =  new double[0][0];
			return record;
		}
		int cls_ind = examples[0].length - 1;
		for(int i=0;i<examples.length;i++)
		{
			if(examples[i][cls_ind] == cls)
			{
				al.add(examples[i]);
			}
		}
		record = new double[al.size()][examples[0].length];
		//record = al.toArray();
		Iterator<double[]> iterator = al.iterator();
		int j = 0;
		while (iterator.hasNext()) {
			record[j] = iterator.next();
			j++;
		}
		return record;
	}

	public static double[] getClasses(double[][] data)
	{
		double[] ary = null;
		int cls_ind = data[0].length - 1;
		Set<Double> set = new TreeSet<Double>();
		for(int i=0;i<data.length;i++)
		{
			set.add(data[i][cls_ind]);
		}
		ary = new double[set.size()];
		Iterator iterator = set.iterator();
		int j = 0;
		while (iterator.hasNext()) {
			ary[j] = (double) iterator.next();
			j++;
		}
		return ary;
	}
	
	public static Integer[] getAttributes(double[][] train) {
		// TODO Auto-generated method stub
		int atr_num = train[0].length-1;
		Integer[] atr = new Integer[atr_num];
		for(int i=0;i<atr_num;i++)
		{
			atr[i]=i;
		}
		return atr;
	}
	
	public static double log2(double n)
	{
	    return (Math.log(n) / Math.log(2));
	}
	
	public static boolean check_single_class(double[][] data)
	{
		Set<Double> set = new TreeSet<Double>();
		int cls_ind = data[0].length - 1;
		for(int i=0;i<data.length;i++)
		{
			set.add(data[i][cls_ind]);
		}
		if(set.size()==1)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	
//	public static int get_single_class(double[][] data)
//	{
//		Set<Double> set = new TreeSet<Double>();
//		int cls_ind = data[0].length - 1;
//		retrun
//	}
}
