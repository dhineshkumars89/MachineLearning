import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Reader;
//import java.nio.file.Files;
//import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.concurrent.TimeUnit;

public class Gaussian_1d {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		//"D:\CSE 6363 Machine Learning\Assignment\Assignment 2\yeast_training.txt"
		
		long startTime = System.currentTimeMillis();
		
		try {
			//List<String> lines = Files.readAllLines(Paths.get(args[0]));
			BufferedReader buffReader = new BufferedReader (new FileReader(args[0]));
			String s = buffReader.readLine();
            
            
//			Iterator<String> it = lines.iterator();
		    HashMap<Integer,ArrayList> hm = new HashMap();
		    HashMap<Integer,Integer> clas2indx = new HashMap();
		    int col = 0;
			
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
			  s = buffReader.readLine();
			}while( s != null );
			
			double[][] total = new double[hm.size()][col];

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
						total[k][i] = total[k][i] + Double.valueOf(str[i]);
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
			
			
			//Calculate variance
			k=0;
			Iterator entry1 = hm.entrySet().iterator();
			while (entry1.hasNext()) {
			    Map.Entry entry = (Map.Entry) entry1.next();
			    ArrayList arylst = (ArrayList) entry.getValue();
			    Object key = entry.getKey();
			    double[][] variance = new double[hm.size()][col];
			    
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
						variance[k][i] = ((variance[k][i])/(arylst.size()-1));
						System.out.println(String.format("Class %d, dimension %d, mean = %.2f, variance = %.2f",key,i+1,total[clas2indx.get(key)][i],variance[k][i]));
					}
				}
				k++;
			}
			
			
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
