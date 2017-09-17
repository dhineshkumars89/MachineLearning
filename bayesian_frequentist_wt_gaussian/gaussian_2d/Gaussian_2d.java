import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

public class Gaussian_2d {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		//"D:\CSE 6363 Machine Learning\Assignment\Assignment 2\yeast_training.txt"
long startTime = System.currentTimeMillis();
		
		try {
//			List<String> lines = Files.readAllLines(Paths.get(args[0]));
//			Iterator<String> it = lines.iterator();
			BufferedReader buffReader = new BufferedReader (new FileReader(args[0]));
			String s = buffReader.readLine();
		    HashMap<Integer,ArrayList> hm = new HashMap();
		    HashMap<Integer,Integer> clas2indx = new HashMap();
		    int col = 0;
			
			do {
			  String[] bit = s.trim().split("\\s+");
			  int classtype = Integer.valueOf(bit[bit.length-1]);
			  col = 2;
			  ArrayList lm;
			  
			  if(hm.containsKey(classtype))
			  {
				  lm = hm.get(classtype);
				  lm.add(Arrays.copyOfRange(bit,0,col));
				  hm.put(classtype, lm);
			  }
			  else
			  {
				  lm = new ArrayList();
				  lm.add(Arrays.copyOfRange(bit,0,col));
				  hm.put(classtype, lm);
			  }
			  s = buffReader.readLine();
			}while(s != null);
			
			double[][] total = new double[hm.size()][col];
			double[][] data = new double[hm.size()][col];
			double[][] sigma = new double[col][col];

			int k=0;
			Iterator entries = hm.entrySet().iterator();
			while (entries.hasNext()) {
			    Map.Entry entry = (Map.Entry) entries.next();
			    ArrayList arylst = (ArrayList) entry.getValue();
			    
			    clas2indx.put((Integer) entry.getKey(), k);
				Iterator its = arylst.iterator();
				while(its.hasNext())
				{
					String[] str = (String[]) its.next();
					for(int i=0;i<col;i++)
					{
						total[k][i] = total[k][i] + Double.parseDouble(str[i]);
						data[k][i] = Double.parseDouble(str[i]);
					}
				}
				if(!its.hasNext())
				{
					for(int i=0;i<col;i++)
					{
						total[k][i] = ((total[k][i])/(arylst.size()));
					}
				}
				k++;
			}
			//System.out.println(data);
			
			//Calculate variance
//			for(int r=0;r<=data.length;r++)
//			{
			int r=0;
			
			Iterator entrs = hm.entrySet().iterator();
			while (entrs.hasNext()) {
			    Map.Entry entry = (Map.Entry) entrs.next();
				
				ArrayList ary = (ArrayList) entry.getValue();
				
				
				Object key = entry.getKey();
//				if((Integer)key==4)
//				{
				for(int i=0;i<col;i++)
				{
					for(int j=0;j<col;j++)
					{
						Iterator its = ary.iterator();
						while(its.hasNext())
						{
							String[] str = (String[]) its.next();
							
							sigma[i][j]= sigma[i][j] + ( (Double.parseDouble(str[i])-total[clas2indx.get(key)][i]) * (Double.parseDouble(str[j])-total[clas2indx.get(key)][j]) );
							//System.out.println((Double.valueOf(str[i])-total[clas2indx.get(key)][i]) * (Double.valueOf(str[j])-total[clas2indx.get(key)][j]));
							//System.out.println("Class "+entry.getKey()+" sigma = "+ (Double.valueOf(str[i])-total[clas2indx.get(key)][i]) * (Double.valueOf(str[j])-total[clas2indx.get(key)][j]));
						}
						
						sigma[i][j] = (double)(sigma[i][j] / (ary.size()-1));
						
//						System.out.println(sigma[i][j]);
					}
				}
				System.out.println(String.format("Class %d, mean = [%.2f, %.2f], sigma = [%.2f, %.2f, %.2f, %.2f]", key,total[r][0],total[r][1],sigma[0][0],sigma[0][1],sigma[1][0],sigma[1][1]));
				sigma[0][0]=0;
				sigma[0][1]=0;
				sigma[1][0]=0;
				sigma[1][1]=0;
				r++;
//				}
			}
			
//			
//			for(int r=0;r<data.length;r++)
//			{
//			System.out.println(String.format("Class %d, mean = [%.2f, %.2f], sigma = [%.2f, %.2f, %.2f, %.2f]",k+1,total[k][i],total[k][i],variance[k][i]));
//			}
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		long endTime   = System.currentTimeMillis();
		 long totalTime = endTime - startTime;
		 //System.out.println(totalTime);
		 //System.out.println(TimeUnit.MILLISECONDS.toSeconds(totalTime));
	}

}
